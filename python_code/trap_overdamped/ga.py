"""
Genetic Algorithm Implementation for Neural Network Protocol Optimization

This module implements a genetic algorithm that evolves neural networks to generate
optimal control protocols for minimal work in stochastic thermodynamic systems.

"""

import torch
import numpy as np
import matplotlib.pyplot as plt
import copy
import logging
import random
from datetime import datetime
from pathlib import Path
from tqdm import tqdm

from nets import FeedforwardNet
from engine_trap_overdamped import calculate_order_parameter, number_of_steps, cee_boundary, visualize_protocol

# Disable gradient computation and set precision
torch.set_grad_enabled(False)
torch.set_default_dtype(torch.float64)

# Configure matplotlib for publication-quality plots
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = 'Computer Modern'
plt.rcParams['text.usetex'] = True  # Enable LaTeX rendering
plt.rcParams['font.size'] = 16
plt.rcParams['xtick.labelsize'] = 16
plt.rcParams['ytick.labelsize'] = 16
plt.rcParams['legend.fontsize'] = 16
plt.rcParams['axes.titlesize'] = 16
plt.rcParams['axes.labelsize'] = 16
plt.set_loglevel("error")

class GeneticAlgorithm:
    """
    Genetic Algorithm for optimizing neural network-generated control protocols.

    This class implements a complete evolutionary system that evolves neural networks
    to generate control protocols that minimize work in stochastic thermodynamic systems.
    """

    def __init__(self, config = None):
        """
        Initialize the genetic algorithm with configuration parameters.

        Args:
            config: Configuration dictionary with GA parameters
        """
        # Default configuration
        self.default_config = {
            # Population parameters
            'population_size': 50,
            'num_generations': 100,
            'elite_fraction': 0.1,

            # Mutation parameters
            'mutation_scale': 0.025,

            # Evaluation parameters
            'num_trajectories': 10000,

            # Neural network architecture
            'network_config': {
                'hidden_dim': 64,
                'num_layers': 3,
                'input_dim': 1,
                'output_dim': 1
            },

            # System parameters
            'random_seed': 0,
            'device': 'cuda' if torch.cuda.is_available() else 'cpu',

            # Output parameters
            'save_all_protocols': True,
            'visualization_interval': 10,
            'output_dir': 'ga_protocols'
        }

        # Merge user config with defaults
        self.config = self.default_config.copy()
        if config:
            self._update_config(self.config, config)

        # Set random seeds for reproducibility
        if self.config['random_seed'] is not None:
            torch.manual_seed(self.config['random_seed'])
            np.random.seed(self.config['random_seed'])
            random.seed(self.config['random_seed'])

        # Initialize algorithm state
        self.population = []
        self.fitness_scores = []
        self.current_generation = 0
        self.best_fitness = float('inf')
        self.best_individual = None
        self.fitness_history = []

        # Setup logging
        self._setup_logging()

        # Create output directories
        self._setup_directories()

        self.logger.info(f"Initialized GeneticAlgorithm with config: {self.config}")

    def _update_config(self, base_config, user_config):
        """Recursively update configuration dictionary."""
        for key, value in user_config.items():
            if key in base_config and isinstance(base_config[key], dict) and isinstance(value, dict):
                self._update_config(base_config[key], value)
            else:
                base_config[key] = value

    def _setup_logging(self):
        """Setup logging configuration."""
        log_dir = Path(self.config['output_dir']) / 'logs'
        log_dir.mkdir(parents=True, exist_ok=True)

        log_file = log_dir / f"ga_optimization_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"

        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(log_file),
                logging.StreamHandler()
            ]
        )
        self.logger = logging.getLogger('GeneticAlgorithm')

    def _setup_directories(self):
        """Create necessary output directories."""
        base_dir = Path(self.config['output_dir'])

        # Create main directories
        for subdir in ['generations', 'best_networks', 'results', 'protocols']:
            (base_dir / subdir).mkdir(parents=True, exist_ok=True)

    def initialize_population(self):
        """Initialize population with random neural networks."""
        self.logger.info(f"Initializing population of size {self.config['population_size']}")

        self.population = []
        net_config = self.config['network_config']

        for i in range(self.config['population_size']):
            # Create neural network
            net = FeedforwardNet(
                n_in=net_config['input_dim'],
                hidden_dim=net_config['hidden_dim'],
                num_layers=net_config['num_layers'],
                n_out=net_config['output_dim']
            )

            # Initialize with random weights
            for param in net.parameters():
                param.normal_(std=self.config['mutation_scale'])

            self.population.append(net)

        # Initialize fitness scores
        self.fitness_scores = [float('inf')] * self.config['population_size']

        self.logger.info("Population initialization complete")

    def generate_protocol(self, network: FeedforwardNet):
        """
        Generate a control protocol from a neural network.

        Args:
            network: Neural network to generate protocol

        Returns:
            Protocol tensor compatible with engine
        """
        # Create protocol tensor with boundary conditions
        protocol = torch.zeros(number_of_steps + 2, 1)
        protocol[0, 0] = cee_boundary[0, 0]   # Initial boundary
        protocol[-1, 0] = cee_boundary[1, 0]  # Final boundary

        # Generate time input for neural network
        time_input = torch.linspace(0, 1, number_of_steps).unsqueeze(1)

        # Forward pass through network
        with torch.no_grad():
            protocol_values = network(time_input)
            protocol[1:-1, 0] = protocol_values[:, 0]

        return protocol

    def save_protocol(self, protocol: torch.Tensor, filename: str):
        """
        Save protocol to file in engine-compatible format.

        Args:
            protocol: Protocol tensor to save
            filename: Output filename
        """
        # Extract middle values (excluding boundaries)
        protocol_data = protocol[1:-1, 0].cpu().numpy()
        np.savetxt(filename, protocol_data)

    def evaluate_individual(self, individual_id: int, network: FeedforwardNet):
        """
        Evaluate fitness of a single individual.

        Args:
            individual_id: Index of individual in population
            network: Neural network to evaluate

        Returns:
            Fitness score (lower is better)
        """
        try:
            # Generate protocol
            protocol = self.generate_protocol(network)
            # Run simulation
            op = calculate_order_parameter(protocol, self.config['num_trajectories'])

            # Calculate fitness (mean work)
            fitness = op

            # Save protocol if requested
            if self.config['save_all_protocols']:
                protocol_file = (Path(self.config['output_dir']) / 'protocols' /
                               f"gen_{self.current_generation:03d}_ind_{individual_id:03d}.dat")
                self.save_protocol(protocol, str(protocol_file))

            return fitness

        except Exception as e:
            self.logger.warning(f"Evaluation failed for individual {individual_id}: {e}")
            return float('inf')

    def evaluate_population(self):
        """Evaluate fitness for all individuals in the population."""
        self.logger.info(f"Evaluating population for generation {self.current_generation}")

        # Reset fitness scores
        self.fitness_scores = []

        # Evaluate each individual
        for i, network in enumerate(tqdm(self.population, desc="Evaluating individuals")):
            fitness = self.evaluate_individual(i, network)
            self.fitness_scores.append(fitness)


        # Update best individual
        best_idx = np.argmin(self.fitness_scores)
        if self.fitness_scores[best_idx] < self.best_fitness:
            self.best_fitness = self.fitness_scores[best_idx]
            self.best_individual = copy.deepcopy(self.population[best_idx])

            # Save best network
            best_net_file = (Path(self.config['output_dir']) / 'best_networks' /
                           f"best_gen_{self.current_generation:03d}.pt")
            self.best_individual.save_model(str(best_net_file))

        self.logger.info(f"Generation {self.current_generation} - Best fitness: {self.best_fitness:.6f}")

    def mutate_net(self, network: FeedforwardNet):
        """
        Create a mutated copy of a neural network.

        Args:
            network: Network to mutate

        Returns:
            Mutated network copy
        """
        mutated_net = copy.deepcopy(network)
        for param in mutated_net.parameters():
            param.add_(self.config['mutation_scale'] * torch.randn_like(param))
        return mutated_net

    def create_offspring(self):
        """
        Create offspring through mutation.

        Args:
            parent_indices: Indices of selected parents

        Returns:
            List of offspring networks
        """
        offspring = []

        # Keep elite individuals
        elite_count = int(self.config['elite_fraction'] * self.config['population_size'])
        elite_indices = sorted(range(len(self.fitness_scores)),
                             key=lambda x: self.fitness_scores[x])[:elite_count]

        for idx in elite_indices:
            offspring.append(copy.deepcopy(self.population[idx]))

        # Create mutated offspring to fill remaining population
        while len(offspring) < self.config['population_size']:
            parent_idx = random.choice(elite_indices)
            mutated_child = self.mutate_net(self.population[parent_idx])
            offspring.append(mutated_child)

        return offspring[:self.config['population_size']]

    def evolve_generation(self):
        """Execute one generation of evolution."""
        generation_start_time = datetime.now()

        # Evaluate current population
        self.evaluate_population()

        # Calculate and store statistics
        fitness_stats = {
            'best': min(self.fitness_scores),
            'worst': max(self.fitness_scores),
            'mean': np.mean(self.fitness_scores),
            'std': np.std(self.fitness_scores)
        }
        self.fitness_history.append(fitness_stats)


        # Create next generation
        self.population = self.create_offspring()

        # Update generation counter
        self.current_generation += 1

        # Track timing
        generation_time = (datetime.now() - generation_start_time).total_seconds()

        self.logger.info(f"Generation {self.current_generation-1} completed in {generation_time:.2f}s")
        self.logger.info(f"Fitness - Best: {fitness_stats['best']:.6f}, "
                        f"Mean: {fitness_stats['mean']:.6f}, "
                      )

    def visualize_progress(self):
        """Generate visualization plots of optimization progress."""
        if not self.fitness_history:
            return

        results_dir = Path(self.config['output_dir']) / 'results'

        # Fitness evolution plot
        fig, (ax1) = plt.subplots(1, 1, figsize=(10, 8))

        generations = range(len(self.fitness_history))
        best_fitness = [stats['best'] for stats in self.fitness_history]
        mean_fitness = [stats['mean'] for stats in self.fitness_history]

        ax1.plot(generations, best_fitness, 'b-', label='Best Fitness', linewidth=2)
        ax1.plot(generations, mean_fitness, 'r--', label='Mean Fitness', linewidth=1)
        ax1.set_xlabel('Generation')
        ax1.set_ylabel('Fitness (Mean Work)')
        ax1.set_title('Fitness Evolution')
        ax1.legend()

        plt.tight_layout()
        plt.savefig(results_dir / 'fitness_evolution.png', dpi=300, bbox_inches='tight')
        plt.close()

        # Best protocol visualization
        if self.best_individual is not None:
            self.visualize_best_protocol()

    def visualize_best_protocol(self):
        """Visualize the best protocol found so far using the engine's visualization function."""
        if self.best_individual is None:
            return

        # Generate best protocol
        best_protocol = self.generate_protocol(self.best_individual)

        # Use the comprehensive visualization function from engine.py
        visualize_protocol(best_protocol, 10*self.config['num_trajectories'])

        # Save best protocol to file
        results_dir = Path(self.config['output_dir']) / 'results'
        protocol_file = results_dir / 'best_protocol.dat'
        self.save_protocol(best_protocol, str(protocol_file))

    def run(self):
        """
        Run the complete genetic algorithm optimization.

        Args:
            resume_from_checkpoint: Optional checkpoint filename to resume from

        Returns:
            Tuple of (best_network, best_fitness)
        """
        self.logger.info("Starting genetic algorithm optimization")

        try:
            self.initialize_population()

            # Main evolution loop
            target_generations = self.config['num_generations']

            while self.current_generation < target_generations:
                # Evolve one generation
                self.evolve_generation()

                # Generate visualizations periodically
                if self.current_generation % self.config['visualization_interval'] == 0:
                    self.visualize_progress()

            # Final statistics and visualization
            self.visualize_progress()

            self.logger.info("Optimization completed successfully")
            self.logger.info(f"Best fitness achieved: {self.best_fitness:.6f}")

            return self.best_individual, self.best_fitness

        except KeyboardInterrupt:
            self.logger.info("Optimization interrupted by user")
            raise
        except Exception as e:
            self.logger.error(f"Optimization failed: {e}")
            raise


def main():
    """Main function demonstrating genetic algorithm usage."""
    # Example configuration
    # config = {
    #     'population_size': 30,
    #     'num_generations': 50,
    #     'num_parents': 6,
    #     'mutation_scale': 0.05,
    #     'num_trajectories': 1000,
    #     'network_config': {
    #         'hidden_dim': 32,
    #         'num_layers': 2
    #     },
    #     'checkpoint_interval': 5,
    #     'visualization_interval': 5,
    #     'output_dir': 'ga_optimization_results'
    # }

    # Create and run genetic algorithm
    # ga = GeneticAlgorithm(config)

    ga = GeneticAlgorithm()
    best_network, best_fitness = ga.run()

    print(f"\nOptimization completed!")
    print(f"Best fitness: {best_fitness:.6f}")

if __name__ == "__main__":
    main()
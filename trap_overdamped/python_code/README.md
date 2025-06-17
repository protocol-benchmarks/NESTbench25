## Code Structure

- `engine_trap_overdamped.py`: Main simulation engine that implements the overdamped Langevin dynamics
- `ga.py`: Implements a genetic algorithm for protocol optimization using neural networks
- `input_control_parameters_learned.dat`: A pre-learned protocol
- `requirements.txt`: Python package dependencies

## Usage

The simulation models a particle in a time-dependent harmonic potential U(x,t) = 0.5*(x-λ(t))², where λ(t) is a time-dependent control parameter that can be optimized.

### Running the Simulation Engine

To run the simulation engine with default settings:
```python
python engine_trap_overdamped.py
```

This will:
1. Load the default protocol
2. Visualize the protocol with 10^5 trajectories
3. Calculate the mean work using 10^6 trajectories

### Key Functions

- `load_protocol(filename)`: Loads a control protocol from a file
- `load_default_protocol()`: Generates the optimal protocol with modified ramp behavior
- `run_protocol(protocol, number_of_trajectories, visualize=False)`: Simulates particle trajectories
- `visualize_protocol(protocol, number_of_trajectories)`: Creates visualizations of the protocol execution
- `final_answer(protocol)`: Calculates the order parameter (mean work) with high precision

### Protocol Optimization

The genetic algorithm can be used to optimize the control protocol:
```python
python ga.py
```

This will:
1. Initialize a population of neural networks
2. Evolve the networks to generate optimal protocols
3. Save the best protocol and visualization results

## Dependencies

Install required packages with:
```
pip install -r requirements.txt
```

Required packages:
- matplotlib (>=3.10.3)
- numpy (>=2.3.0)
- torch (>=2.7.1)


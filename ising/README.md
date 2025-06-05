# Ising Model Benchmark for Protocol Design

This code implements a benchmark system for testing methods of protocol optimization in nonequilibrium statistical mechanics. The system consists of a 2D Ising model undergoing state change, with the objective of minimizing entropy production.

## Problem Description

We consider the 2D Ising model on a square lattice of $N=32^2$ sites with periodic boundary conditions. On each site $i$ is a binary spin $S_i = \pm 1$. The lattice energy function is:

$$E = -J\sum_{\langle ij \rangle} S_i S_j -h\sum_{i=1}^N S_i$$

where $J=1$ is the Ising coupling, $h$ is the magnetic field, and the first sum runs over all nearest-neighbor bonds.

The model evolves by Glauber Monte Carlo dynamics. At each time step, a lattice site is chosen and a spin flip is proposed, which is accepted with probability:

$$P_{\text{Glauber}}(\Delta E)=\left( 1+ \exp(\Delta E/T) \right)^{-1}$$

where $\Delta E$ is the energy change under the proposed move.

## Optimization Objective

The control parameters of this problem are temperature and magnetic field, $\mathbf{c}(t) = (T(t), h(t))$. The aim is to change these parameters from the values $\mathbf{c}(0) = (0.7,-1)$ to the values $\mathbf{c}(t_f) = (0.7,1)$, in $t_f=100 N$ Monte Carlo steps (100 Monte Carlo sweeps). 

The order parameter to be minimized is the mean entropy production $\langle\sigma\rangle$, where the entropy produced over the course of a simulation is:

$$\sigma=E_f/T_f - E_0/T_0 - \sum_{k=1}^{t_f} \Delta E_k/T_k$$

Here $E_f$ and $E_0$ are the final and initial energies of the system, and $\Delta E_k$ and $T_k$ are the energy change and temperature at step $k$ of the simulation.

## Default Protocol

The default protocol is based on the work of Rotskoff and Crooks (2015), which is the optimal protocol within the near-equilibrium approximation. This protocol avoids the first-order phase transition line and the critical point of the Ising model.

## Code Structure

The code consists of:
- `engine_ising.c`: Main simulation code
- `engine_ising.h`: Header file with API definitions
- `movie.py`: Visualization script to generate snapshots and animations
- `default_protocol.dat`: Default protocol parameters
- `input_control_parameters_learned.dat`: Example of an optimized protocol

## Usage

### Building the Code

The code can be built either as a standalone executable or as a library for integration with optimization algorithms:

```bash
# Build as standalone executable
make standalone

# Build as library
make library
```

### Running the Simulation

To run the simulation with the default protocol:

```bash
./sim
```

This will:
1. Load the default protocol
2. Run the simulation for 100 Monte Carlo sweeps
3. Calculate the entropy production
4. Generate visualization files (`output_picture.pdf`, `output_movie.mp4`)

### Using a Custom Protocol

To use a custom protocol:
1. Create a file named `input_control_parameters.dat` with dimensions $2 \times 1000$, where each row contains the temperature and magnetic field values at 1000 equally spaced time points.
2. Run the simulation as above.

### API Functions

The following functions are exposed in the library:

- `final_answer()`: Calculates the order parameter (entropy production) over 100 samples of 10^3 trajectories
- `load_protocol()`: Reads protocol from `input_control_parameters.dat`
- `load_default_protocol()`: Loads the default protocol
- `visualize_protocol()`: Generates visualizations of the current protocol
- `calculate_order_parameter(n_traj)`: Computes order parameter over n_traj trajectories

## Visualization

The simulation generates several output files:

- `report_lattice_time_*.dat`: Lattice state snapshots at different times
- `report_entprod_histogram.dat`: Entropy production histogram
- `output_picture.pdf`: Grid of lattice snapshots showing system evolution
- `output_movie.mp4`: Animation of the simulation

## Protocol Learning

This benchmark is designed to test methods for optimizing time-dependent protocols. The file `input_control_parameters_learned.dat` contains an example of a protocol learned via neuroevolution that produces significantly less entropy than the default protocol.

For this system, small perturbations about a simple protocol produce little change in the microscopic state because at $T=0.7$ the energy barrier to flipping single spins is more than 10 (in units of $k_BT$). Effective optimization strategies include starting with large random perturbations or using informed initial guesses.

## Requirements

- C++ compiler with C++11 support
- Python 3.6+ with numpy, matplotlib, and imageio (for visualization)
- GNU Make

## References

This benchmark is based on the work described in:
- Rotskoff, G. M., & Crooks, G. E. (2015). Optimal control in nonequilibrium systems: Dynamic Riemannian geometry of the Ising model. Physical Review E, 92(6), 060102.
- Gingrich, T. R., Rotskoff, G. M., Vaikuntanathan, S., & Geissler, P. L. (2016). Efficiency and large deviations in time-asymmetric stochastic heat engines. New Journal of Physics, 18(2), 023007.
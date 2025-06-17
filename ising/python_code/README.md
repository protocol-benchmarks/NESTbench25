## Code Structure

The code consists of:
- `engine_ising.py`: Main simulation code implementing 2D Ising model dynamics
- `default_protocol.dat`: Default protocol parameters for temperature and magnetic field

## Usage

### Running the Simulation

To run the simulation with the default protocol:

```bash
python engine_ising.py
```

This will:
1. Load the default protocol
2. Run the simulation for 100 Monte Carlo sweeps
3. Calculate the entropy production
4. Generate visualization files (`protocol.png`, `entprod_histogram.png`, `output.mp4`)

### Using a Custom Protocol

To use a custom protocol:
1. Create a file named `input_control_parameters.dat` with dimensions $2 \times N$, where each row contains the temperature and magnetic field values at N equally spaced time points.
2. Modify the code to load your custom protocol using `load_protocol()` instead of `load_default_protocol()`.

### API Functions

The following functions are available in the Python module:

- `load_protocol(filename)`: Reads protocol from a specified file
- `load_default_protocol()`: Loads the default protocol
- `run_protocol(protocol, number_of_trajectories, visualize)`: Simulates the Ising model with given protocol
- `visualize_protocol(protocol, number_of_trajectories)`: Generates visualizations of the current protocol
- `calculate_order_parameter(protocol, number_of_trajectories)`: Computes order parameter over n_traj trajectories
- `final_answer(protocol)`: Calculates the order parameter (entropy production) with high statistical precision

## Visualization

The simulation generates several output files:

- `protocol.png`: Visualization of the temperature and magnetic field protocol
- `plots/output_*.png`: Lattice state snapshots at different times
- `entprod_histogram.png`: Entropy production histogram
- `output.mp4`: Animation of the simulation showing spin configurations over time

## Protocol Learning

This benchmark is designed to test methods for optimizing time-dependent protocols. The default protocol is based on the work of Rotskoff and Crooks (2015), which is the optimal protocol within the near-equilibrium approximation.

For this system, small perturbations about a simple protocol produce little change in the microscopic state because at $T=0.7$ the energy barrier to flipping single spins is more than 10 (in units of $k_BT$). Effective optimization strategies include starting with large random perturbations or using informed initial guesses.

## Implementation Details

The Python implementation uses PyTorch for efficient tensor operations and supports GPU acceleration when available. The simulation employs a checkerboard updating scheme for the Metropolis algorithm to ensure proper simulation of the Ising model dynamics.

## Requirements

- Python 3.6+
- PyTorch
- NumPy
- Matplotlib
- FFmpeg (for movie generation)
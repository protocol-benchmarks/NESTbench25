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


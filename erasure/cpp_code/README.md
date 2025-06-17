## Code Structure

The code consists of:
- `engine_erasure.c`: Main simulation code
- `engine_erasure.h`: Header file with API definitions
- `movie.py`: Visualization script to generate snapshots and animations
- `input_control_parameters_learned.dat`: Example of an optimized protocol
- `Makefile`: Build configuration

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
2. Run the simulation for the specified erasure time
3. Calculate the erasure error rate
4. Generate visualization files (`output_picture.pdf`, `output_movie.mp4`)

### Using a Custom Protocol

To use a custom protocol:
1. Create a file named `input_control_parameters.dat` with dimensions $2 \times 1000$, where each row contains the potential shape parameters $c_0$ and $c_1$ at 1000 equally spaced time points.
2. Run the simulation as above.

### API Functions

The following functions are exposed in the library:

- `final_answer()`: Calculates the order parameter (erasure error rate) over 100 samples of 10^4 trajectories
- `load_protocol()`: Reads protocol from `input_control_parameters.dat`
- `load_default_protocol()`: Loads the default protocol
- `visualize_protocol()`: Generates visualizations of the current protocol
- `calculate_order_parameter(n_traj)`: Computes order parameter over n_traj trajectories
- `calculate_custom_order_parameter(n_traj)`: Alternative order parameter for curriculum learning

## Visualization

The simulation generates several output files:

- `report_work_histogram.dat`: Work distribution data
- `output_picture.pdf`: Grid of snapshots showing system evolution
- `output_movie.mp4`: Animation of the simulation

The visualization shows the potential $U(x,\mathbf{c}(t))$, the associated Boltzmann distribution $\rho_0(x)$, and the distribution of particle positions $\rho(x,t)$.

## Protocol Learning

This benchmark is designed to test methods for optimizing time-dependent protocols. The file `input_control_parameters_learned.dat` contains an example of a protocol learned via neuroevolution that achieves a much higher erasure success rate than the default protocol.

For this system, the challenge is to manipulate the potential in a way that effectively transfers the particle to the left well without imparting excessive kinetic energy that would allow it to escape. The system includes a quiescent period of duration $2t_f$ after the protocol has run to check the stability of the memory.

## Requirements

- C++ compiler with C++11 support
- Python 3.6+ with numpy, matplotlib, and imageio (for visualization)
- GNU Make

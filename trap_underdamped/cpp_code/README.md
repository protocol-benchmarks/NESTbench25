## Code Structure

The code consists of:
- `engine_trap_underdamped.c`: Main simulation code
- `engine_trap_underdamped.h`: Header file with API definitions
- `movie.py`: Visualization script to generate snapshots and animations
- `Makefile`: Build instructions
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
1. Load the default protocol (the optimal protocol)
2. Run the simulation for the specified time
3. Calculate the work required
4. Generate visualization files (`output_figure.pdf`, `output_movie.mp4`)

### Using a Custom Protocol

To use a custom protocol:
1. Create a file named `input_control_parameters.dat` with dimensions $1 \times 1000$, where each row contains the trap position values at 1000 equally spaced time points.
2. Run the simulation as above.

### API Functions

The following functions are exposed in the library:

- `final_answer()`: Calculates the order parameter (mean work) over 100 samples of 10^4 trajectories
- `load_protocol()`: Reads protocol from `input_control_parameters.dat`
- `load_default_protocol()`: Loads the optimal protocol
- `visualize_protocol()`: Generates visualizations of the current protocol
- `calculate_order_parameter(n_traj)`: Computes order parameter over n_traj trajectories

## Visualization

The simulation generates several output files:

- `report_position_time_*.dat`: Position distribution snapshots at different times
- `report_work_histogram.dat`: Work distribution histogram
- `output_figure.pdf`: Grid of snapshots showing system evolution
- `output_movie.mp4`: Animation of the simulation

## Protocol Learning

This benchmark is designed to test methods for optimizing time-dependent protocols. The file `input_control_parameters_learned.dat` contains an example of a protocol learned via neuroevolution that produces work close to the theoretical optimum.

The challenge in this system is to approximate the effect of delta function impulses within the constraints of the protocol representation. The learned protocol uses rapid changes at the beginning and end of the protocol to approximate these impulses.

## Requirements

- C++ compiler with C++11 support
- Python 3.6+ with numpy, matplotlib, and imageio (for visualization)
- GNU Make

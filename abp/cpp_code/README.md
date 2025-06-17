
## Code Structure

The code consists of:
- `engine_abp.c`: Main simulation code
- `engine_abp.h`: Header file with API definitions
- `movie.py`: Visualization script to generate snapshots and animations
- `target_p_r.dat`: Target radial distribution for the active steady state
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
2. Run the simulation for the specified time
3. Calculate the order parameter (transformation accuracy and work)
4. Generate visualization files (`output_figure.pdf`, `output_movie.mp4`)

### Using a Custom Protocol

To use a custom protocol:
1. Create a file named `input_control_parameters.dat` with dimensions $2 \times 1000$, where each row contains the trap stiffness $\kappa$ and swim speed $\lambda$ values at 1000 equally spaced time points.
2. Run the simulation as above.

### API Functions

The following functions are exposed in the library:

- `final_answer()`: Calculates the order parameter over 10^6 trajectories
- `load_protocol()`: Reads protocol from `input_control_parameters.dat`
- `load_default_protocol()`: Loads the default protocol
- `visualize_protocol()`: Generates visualizations of the current protocol
- `calculate_order_parameter(n_traj)`: Computes order parameter over n_traj trajectories

## Visualization

The simulation generates several output files:

- `report_position_time_*.dat`: Radial distribution snapshots at different times
- `report_work_histogram.dat`: Work distribution histogram
- `output_figure.pdf`: Grid of snapshots showing system evolution
- `output_movie.mp4`: Animation of the simulation

## Protocol Learning

This benchmark is designed to test methods for optimizing time-dependent protocols. The file `input_control_parameters_learned.dat` contains an example of a protocol learned via neuroevolution that successfully enacts the state-to-state transformation while extracting approximately 3 $k_B T$ units of work from the system.

The learned protocol typically involves an initial input of work as the trap parameter $\kappa$ is increased, with the extraction of work happening later when $\kappa$ is returned to its resting value.

## Requirements

- C++ compiler with C++11 support
- Python 3.6+ with numpy, matplotlib, and imageio (for visualization)
- GNU Make


## Code Structure

The code consists of:
- `engine_trap_underdamped.py`: Main simulation code with both physics engine and visualization capabilities
- `input_control_parameters_learned.dat`: Example of an optimized protocol

## Usage

### Running the Simulation

To run the simulation with the default protocol:

```bash
python engine_trap_underdamped.py
```

This will:
1. Load the default protocol (the optimal protocol)
2. Run the simulation for the specified time
3. Calculate the work required
4. Generate visualization files (`time_evolution.png`, `work_histogram.png`, `output.mp4`)

### Using a Custom Protocol

To use a custom protocol:
1. Create a file named `input_control_parameters_learned.dat` with dimensions $1 \times 1000$, where each row contains the trap position values at 1000 equally spaced time points.
2. Modify the main block in `engine_trap_underdamped.py` to use `load_protocol()` instead of `load_default_protocol()`.
3. Run the simulation as above.

### API Functions

The following functions are available in the Python module:

- `final_answer(protocol, kick_velocity=False)`: Calculates the order parameter (mean work) over 100 samples of 10^4 trajectories
- `load_protocol(filename="input_control_parameters_learned.dat")`: Reads protocol from file with automatic handling of different protocol lengths
- `load_default_protocol()`: Loads the analytically optimal protocol
- `visualize_protocol(protocol, number_of_trajectories, kick_velocity=False)`: Generates visualizations of the current protocol
- `calculate_order_parameter(protocol, number_of_trajectories)`: Computes order parameter over n_traj trajectories
- `run_protocol(protocol, number_of_trajectories, kick_velocity=False, visualize=False)`: Core simulation function that runs the Langevin dynamics

## Visualization

The simulation generates several output files:

- `protocol.png`: Visualization of the control protocol
- `time_evolution.png`: Grid of snapshots showing system evolution at different time points
- `work_histogram.png`: Work distribution histogram
- `output.mp4`: Animation of the simulation
- Individual frames in the `plots/` directory

## Protocol Learning

This benchmark is designed to test methods for optimizing time-dependent protocols. The file `input_control_parameters_learned.dat` contains an example of a protocol learned via neuroevolution that produces work close to the theoretical optimum.

The challenge in this system is to approximate the effect of delta function impulses within the constraints of the protocol representation. The learned protocol uses rapid changes at the beginning and end of the protocol to approximate these impulses.

## Implementation Details

The simulation implements underdamped Langevin dynamics with a time-dependent harmonic potential. Key features include:

- GPU acceleration via PyTorch when available
- Efficient batch processing of multiple trajectories in parallel
- High-quality visualization with LaTeX rendering
- Statistical analysis of work distribution with error estimation

## Requirements

- Python 3.6+ with PyTorch, NumPy, and Matplotlib
- ffmpeg (for video generation)
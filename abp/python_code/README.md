## Code Structure

The code consists of:
- `engine_abp.py`: Main simulation code with API definitions and visualization functions
- `input_control_parameters_learned.dat`: Example of an optimized protocol
- `target_p_r.dat`: Target radial distribution for the active steady state

## Usage

### Running the Simulation

To run the simulation with the default protocol:

```bash
python engine_abp.py
```

This will:
1. Load the default protocol
2. Run the simulation for the specified time
3. Calculate the order parameter (transformation accuracy and work)
4. Generate visualization files (`time_evolution.png`, `work_histogram.png`, `output.mp4`)

### Using a Custom Protocol

To use a custom protocol:
1. Create a file named `input_control_parameters.dat` with dimensions $2 \times N$, where each row contains the trap stiffness $\kappa$ and swim speed $\lambda$ values. The protocol will be automatically interpolated to match the required number of timesteps.
2. Run the simulation as above.

### API Functions

The following functions are available in the Python module:

- `final_answer(protocol)`: Calculates the order parameter over 10^6 trajectories
- `load_protocol(filename)`: Reads protocol from a specified file
- `load_default_protocol()`: Loads the default protocol
- `visualize_protocol(protocol, number_of_trajectories)`: Generates visualizations of the current protocol
- `calculate_order_parameter(protocol, number_of_trajectories)`: Computes order parameter over n_traj trajectories
- `run_protocol(protocol, number_of_trajectories, visualize)`: Runs the simulation with given parameters

## Visualization

The simulation generates several output files:

- `plots/output_*.png`: Radial distribution snapshots at different times
- `work_histogram.png`: Work distribution histogram
- `time_evolution.png`: Grid of snapshots showing system evolution
- `output.mp4`: Animation of the simulation
- `protocol.png`: Visualization of the control parameters over time

## Protocol Learning

This benchmark is designed to test methods for optimizing time-dependent protocols. The file `input_control_parameters_learned.dat` contains an example of a protocol learned via neuroevolution that successfully enacts the state-to-state transformation while extracting approximately 3 $k_B T$ units of work from the system.

The learned protocol typically involves an initial input of work as the trap parameter $\kappa$ is increased, with the extraction of work happening later when $\kappa$ is returned to its resting value.

## Requirements

- Python 3.6+ with:
  - PyTorch (with CUDA support recommended for faster simulation)
  - NumPy
  - Matplotlib (with LaTeX support for better visualization)
  - FFmpeg (for movie generation)
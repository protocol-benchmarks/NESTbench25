## Code Structure

The code consists of:
- `engine_erasure.py`: Main simulation code for stochastic thermodynamic system
- `input_control_parameters_learned.dat`: Example of an optimized protocol

## Usage

### Running the Simulation

To run the simulation with the default protocol:

```bash
python engine_erasure.py
```

This will:
1. Generate a default protocol
2. Run the simulation for the specified erasure time
3. Calculate the erasure error rate
4. Generate visualization files (`time_evolution.png`, `output.mp4`, `work_histogram.png`, `protocol.png`)

### Using a Custom Protocol

To use a custom protocol:
1. Create a file named `input_control_parameters.dat` with dimensions $2 \times 1000$, where each row contains the potential shape parameters $c_0$ and $c_1$ at 1000 equally spaced time points.
2. Modify the main section of `engine_erasure.py` to load this protocol:
   ```python
   protocol = load_protocol("input_control_parameters.dat")
   ```
3. Run the simulation as above.

### API Functions

The following functions are available in the module:

- `final_answer(protocol)`: Calculates the order parameter (erasure error rate) over 100 samples of 10^4 trajectories
- `load_protocol(filename)`: Reads protocol from the specified file with automatic interpolation
- `load_default_protocol()`: Generates a default protocol with piecewise constant behavior
- `visualize_protocol(protocol, number_of_trajectories)`: Generates visualizations of the current protocol
- `calculate_order_parameter(protocol, number_of_trajectories)`: Computes order parameter over n_traj trajectories
- `run_protocol(protocol, number_of_trajectories, visualize)`: Runs the simulation with the given protocol

## Visualization

The simulation generates several output files:

- `work_histogram.png`: Work distribution data
- `time_evolution.png`: Grid of snapshots showing system evolution
- `output.mp4`: Animation of the simulation
- `protocol.png`: Plot of the control parameters over time

The visualization shows the potential $U(x,\mathbf{c}(t))$, the associated Boltzmann distribution $\rho_0(x)$, and the distribution of particle positions $\rho(x,t)$.

## Protocol Learning

This benchmark is designed to test methods for optimizing time-dependent protocols. The file `input_control_parameters_learned.dat` contains an example of a protocol learned via neuroevolution that achieves a much higher erasure success rate than the default protocol.

For this system, the challenge is to manipulate the potential in a way that effectively transfers the particle to the left well without imparting excessive kinetic energy that would allow it to escape. The system includes a quiescent period of duration $2t_f$ after the protocol has run to check the stability of the memory.

## System Description

The system simulates an underdamped Langevin equation with parameters based on a mechanical oscillator with frequency f₀ = 1090 Hz and quality factor Q = 7. The potential is defined as:

$U(x,t) = 0.5*(x - \text{sign}(x-c₀(t))*c₁(t))² + c₀(t)*c₁(t)*(\text{sign}(x-c₀(t)) + \text{sign}(c₀(t)))$

where $c₀(t)$ and $c₁(t)$ are time-dependent control parameters that can be optimized.

## Requirements

- Python 3.6+ with:
  - PyTorch
  - NumPy
  - Matplotlib
  - FFmpeg (for video generation)
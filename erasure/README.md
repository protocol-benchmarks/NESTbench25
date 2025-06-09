# Logic Erasure Benchmark for Protocol Design

This code implements a benchmark system for testing methods of protocol optimization in nonequilibrium statistical mechanics. The system simulates logic erasure using an underdamped Langevin particle in a double-well potential, with the objective of maximizing the probability of successful erasure.

## Problem Description

We consider an underdamped Langevin particle that can start in either well of a double-well potential, representing a one-bit memory (left well = 0, right well = 1). The particle's motion is governed by:

$$m \ddot{x} + \gamma \dot{x} = -\frac{\partial U(x,\mathbf{c}(t))}{\partial x} + \left( 2 \gamma k_B T\right)^{1/2}\xi(t)$$

where $\xi(t)$ is Gaussian white noise with zero mean and correlation function $\langle\xi(t) \xi(t')\rangle = \delta(t-t')$. The potential is:

$$U(x, \mathbf{c}(t)) = \frac{k}{2} \left( x - S(x - c_0(t))\, c_1(t) \right)^2 + k c_0(t)\, c_1(t) \left( S(x - c_0(t)) + S(c_0(t)) \right)$$

where $S$ is the sign function. The system parameters are based on a micromechanical cantilever with:
- Fundamental frequency: $f_0 = 1090$ Hz
- Quality factor: $Q = 7$

## Optimization Objective

The control parameters of this problem are the potential shape parameters $\mathbf{c}(t) = (c_0(t), c_1(t))$. The aim is to manipulate these parameters to bring the particle to the left-hand well (logic state 0) at time $t_f$, regardless of its initial state, thus erasing the memory.

The protocol must operate within a very short time $t_f = 0.5 \times (2\pi/\omega_0)$, which is half the fundamental oscillation period of the system. This short duration makes perfect erasure impossible, and the objective is to minimize the erasure error rate.

The boundary conditions are $\mathbf{c}(0) = \mathbf{c}(t_f) = (0, 5)$, which correspond to a symmetric double-well potential with an energy barrier of $12.5$ (in units of $k_B T$).

## Default Protocol

The default protocol attempts to merge the double wells, translate the resulting single well to the left, and then reconstitute the double-well potential. While this approach works well for longer erasure times, it is ineffective for the short duration specified in this benchmark, achieving only about 53% successful erasure.

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

## References

This benchmark is based on the work described in:
- Dago, S., et al. (2021). Information erasure under fast non-isothermal operations. Physical Review Research, 3(1), 013161.
- Dago, S., et al. (2024). Approaching the physical limits of information erasure in a micromechanical system. Applied Physics Reviews, 11, 021405.
- Landauer, R. (1961). Irreversibility and heat generation in the computing process. IBM Journal of Research and Development, 5(3), 183-191.
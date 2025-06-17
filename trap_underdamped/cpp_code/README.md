# Underdamped Particle in a Harmonic Trap Benchmark

This code implements a benchmark system for testing methods of protocol optimization in nonequilibrium statistical mechanics. The system consists of an underdamped Langevin particle in a moving harmonic trap, with the objective of minimizing work during trap translation.

## Problem Description

We consider an underdamped particle with position $x$ experiencing the Langevin dynamics:

$$m \ddot{x} + \gamma \dot{x} = -\frac{\partial U(x,\lambda(t))}{\partial x} + \left(2 \gamma k_B T\right)^{1/2}\xi(t)$$

where $m$ is the particle mass, $\gamma$ is the damping coefficient, and $\xi(t)$ is Gaussian white noise with zero mean and correlation function $\langle\xi(t)\xi(t')\rangle=\delta(t-t')$.

The particle is confined in a harmonic potential:

$$U(x,\lambda) = \frac{k}{2}(x-\lambda(t))^2$$

where $k$ is the spring constant and $\lambda(t)$ is the time-dependent position of the trap center.

The natural scales of this problem are:
- Length: $\sigma_0 = \sqrt{k_B T/k} = 1$ nm
- Time: $\omega_0^{-1} = \sqrt{m/k}$ where $\omega_0 = 2\pi \times 1090$ Hz
- Quality factor: $Q = m\omega_0/\gamma = 7$

These values are motivated by micromechanical cantilever experiments, so the "particle" in this problem represents the position of a cantilever tip.

## Optimization Objective

The control parameter of this problem is the trap center position $\lambda(t)$. The aim is to move the trap from position $\lambda(0)=0$ to position $\lambda(t_f)=5$ in finite time $t_f = 0.15 \times (2\pi/\omega_0)$, which is short relative to the fundamental oscillation time.

The order parameter to be minimized is the mean work required to perform many independent realizations of this translation.

## Optimal Protocol

The optimal protocol for this problem is known analytically and is given by:

$$\lambda^*(t) = \lambda(t_f)\frac{Q\omega_0 t + 1}{Q\omega_0 t_f + 2} + \Lambda[\delta(t) - \delta(t-t_f)]$$

where $\Lambda = \lambda(t_f)/(\omega_0^2 t_f + 2\omega_0 Q^{-1})$.

This protocol contains delta-function kicks at the beginning and end, which effect an instantaneous change in the velocity of the particle. This poses a challenge for optimization algorithms since direct control of the particle velocity is not allowed, only control of the trap position $\lambda(t)$.

The mean work done under this optimal protocol is:

$$W^* = \frac{k\lambda(t_f)^2}{2 + Q\omega_0 t_f} \approx 2.90787$$

in units of $k_B T$.

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

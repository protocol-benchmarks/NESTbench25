# Active Brownian Particle Benchmark

This code implements a benchmark system for testing methods of protocol optimization in nonequilibrium statistical mechanics. The system consists of an active Brownian particle (ABP) confined in a harmonic trap, with the dual objective of enacting a state-to-state transformation and maximizing extracted work.

## Problem Description

We consider an active Brownian particle in two-dimensional space with position vector $\bm{r}$ and orientation $\theta$. The particle moves in the direction $\hat{\bm{e}}(\theta)=(\cos\theta, \sin\theta)$ with controlled speed $\lambda$. Its motion is governed by the Langevin equations:

$$\dot{\bm{r}} = \lambda(t) \hat{\bm{e}}(\theta) - \kappa(t) \bm{r} + \sqrt{2}\bm{\xi}_r(t)$$

$$\dot{\theta} = \sqrt{2}\xi_\theta(t)$$

where $\bm{\xi}_r(t)$ and $\xi_\theta(t)$ are Gaussian white noise terms with zero mean and unit variance. The force term $-\kappa(t) \bm{r}$ is derived from a harmonic restoring potential:

$$U(\bm{r},\kappa(t))=\frac{1}{2} \kappa(t) \bm{r}^2$$

The spring constant $\kappa(t)$ and swim speed $\lambda(t)$ are the control parameters of the problem, bounded as $1 \leq \kappa \leq 7$ and $0 \leq \lambda \leq 11$. These bounds are inspired by experiments on spherical Janus particles, whose self-propulsion speed is controlled by light intensity, and which are confined by acoustic traps.

## Optimization Objective

This problem has two objectives:

1. **State-to-state transformation**: Start from a passive equilibrium state with parameters $(\kappa(0),\lambda(0)) = (4,0)$ and reach an active steady state with parameters $(\kappa(t_f),\lambda(t_f)) = (4,5)$ in time $t_f = 0.44$.

2. **Work extraction**: Extract as much work as possible during this transformation.

The order parameter to be minimized is $\Delta + 10^{-3} \langle W \rangle$, where:

- $\Delta$ measures the mean-squared error between the final radial distribution $\rho(r,t_f)$ and the target distribution $\rho_{ss}(r)$
- $\langle W \rangle$ is the mean work done on the system (negative values indicate work extraction)

The prefactor $10^{-3}$ assigns slightly greater weight to the state-to-state transformation, making this the primary objective, with work extraction as a secondary goal.

## Default Protocol

The default protocol immediately sets the control parameter values to their final values. This protocol fails to effect a proper state-to-state transformation and does not extract any work from the system, since $\kappa$ is not changed.

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

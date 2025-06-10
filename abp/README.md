# Active Brownian Particles Benchmark for Protocol Design

This code implements a benchmark system for testing methods of protocol optimization in nonequilibrium statistical mechanics. The system consists of an active Brownian particle (ABP) confined in a harmonic trap, with the objective of enacting a state-to-state transformation while maximizing work extraction.

## Problem Description

We consider an active Brownian particle in two-dimensional space. The particle has position vector $\bm{r}$ and orientation $\theta$, and moves in the direction $\hat{\bm{e}}(\theta)=(\cos\theta, \sin\theta)$ with constant speed $\lambda$. Its motion is governed by the Langevin equations:

$$
\begin{align}
\dot{{\bm r}} &= \lambda(t) \hat{\bm{e}}(\theta) - \kappa(t) \bm{r} + \sqrt{2}\bm{\xi}_r(t),\\ 
\dot{\theta} &= \sqrt{2}\xi_\theta(t),
\end{align}
$$

where $\bm{\xi}_r(t)$ and $\xi_\theta(t)$ are Gaussian white noise terms with zero mean and unit variance. The force term $-\kappa(t) {\bm r}$ is derived from a harmonic restoring potential $U({\bm r},\kappa(t))=\frac{1}{2} \kappa(t) {\bm r}^2$.

This model is inspired by experiments on spherical Janus particles, whose self-propulsion speed is controlled by light intensity, and which are confined by acoustic traps.

## Optimization Objective

The control parameters of this problem are the trap spring constant and swim speed, $\mathbf{c}(t) = (\kappa(t), \lambda(t))$. The aim is to change these parameters from the values $\mathbf{c}(0) = (4,0)$ (passive state) to the values $\mathbf{c}(t_f) = (4,5)$ (active state), in time $t_f=0.44$.

The order parameter to be optimized has two components:
1. The primary objective is to ensure that the final radial distribution $\rho(r,t_f)$ matches the steady-state distribution $\rho_{ss}(r)$ associated with the final parameters.
2. The secondary objective is to extract as much work as possible from the system.

These objectives are combined into a single order parameter:
$$\Delta + 10^{-3} \langle W \rangle$$

where $\Delta$ measures the mean-squared error between the actual profile and target profile, and $\langle W \rangle$ is the mean work.

## Default Protocol

The default protocol immediately sets the control parameters to their final values. This protocol fails to effect a state-to-state transformation, as the system does not have enough time to reach the steady state associated with the final parameters. Additionally, no work is extracted because the trap stiffness $\kappa$ remains constant.

## Code Structure

The code consists of:
- `engine_abp.c`: Main simulation code
- `engine_abp.h`: Header file with API definitions
- `movie.py`: Visualization script to generate snapshots and animations
- `target_p_r.dat`: Target radial distribution function
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
2. Run the simulation
3. Calculate the order parameter (state-to-state transformation error and work)
4. Generate visualization files (`output_figure.pdf`, `output_movie.mp4`)

### Using a Custom Protocol

To use a custom protocol:
1. Create a file named `input_control_parameters.dat` with dimensions $2 \times 1000$, where each row contains the trap stiffness $\kappa$ and swim speed $\lambda$ values at 1000 equally spaced time points.
2. Run the simulation as above.

The protocol parameters are bounded as follows:
- Trap stiffness: $1 \leq \kappa \leq 7$
- Swim speed: $0 \leq \lambda \leq 11$

These bounds are enforced by the code and do not need to be enforced in the input file.

### API Functions

The following functions are exposed in the library:

- `final_answer()`: Calculates the order parameter using 10^6 trajectories
- `load_protocol()`: Reads protocol from `input_control_parameters.dat`
- `load_default_protocol()`: Loads the default protocol
- `visualize_protocol()`: Generates visualizations of the current protocol
- `calculate_order_parameter(n_traj)`: Computes order parameter over n_traj trajectories

## Visualization

The simulation generates several output files:

- `report_position_time_*.dat`: Radial distribution snapshots at different times
- `report_boltz_time_*.dat`: Boltzmann distribution at different times
- `report_potential_time_*.dat`: Potential energy at different times
- `report_work_histogram.dat`: Work distribution
- `output_figure.pdf`: Grid of snapshots showing system evolution
- `output_movie.mp4`: Animation of the simulation

## Protocol Learning

This benchmark is designed to test methods for optimizing time-dependent protocols. The file `input_control_parameters_learned.dat` contains an example of a protocol learned via neuroevolution that successfully enacts a state-to-state transformation while extracting work from the system.

The learned protocol demonstrates that it's possible to extract almost 3 $k_BT$ units of work while ensuring that the final state closely matches the target steady state. The work extraction involves an initial input of work as the trap parameter $\kappa$ is increased, followed by extraction of work when $\kappa$ is returned to its resting value.

## Requirements

- C++ compiler with C++11 support
- Python 3.6+ with numpy, matplotlib, and imageio (for visualization)
- GNU Make

## References

This benchmark is based on the work described in:
- Baldovin, M., Vulpiani, A., Puglisi, A., & Baiesi, M. (2023). Control protocols for active particles. Physical Review E, 107(3), 034603.
- Casert, C., Whitelam, S., & Tamblyn, I. (2024). Learning control protocols for active Brownian particles.
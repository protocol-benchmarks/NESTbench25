# Active Brownian Particle Benchmark for Protocol Design

This code implements a benchmark system for testing methods of protocol optimization in nonequilibrium statistical mechanics. The system consists of an active Brownian particle (ABP) in a harmonic trap, with the dual objective of enacting a state-to-state transformation and maximizing extracted work.

## Problem Description

We consider an active Brownian particle in two-dimensional space. The particle has position vector $\mathbf{r}$ and orientation $\theta$, and moves in the direction $\hat{\mathbf{e}}(\theta)=(\cos\theta, \sin\theta)$ with controllable speed $\lambda$. Its motion is governed by the Langevin equations:

$$\dot{\mathbf{r}} = \lambda(t) \hat{\mathbf{e}}(\theta) - \kappa(t) \mathbf{r} + \sqrt{2}\mathbf{\xi}_r(t)$$
$$\dot{\theta} = \sqrt{2}\xi_\theta(t)$$

where $\mathbf{\xi}_r(t)$ and $\xi_\theta(t)$ are Gaussian white noise terms with zero mean and unit variance. The force term $-\kappa(t) \mathbf{r}$ is derived from a harmonic restoring potential $U(\mathbf{r},\kappa(t))=\frac{1}{2} \kappa(t) \mathbf{r}^2$.

## Optimization Objective

The control parameters of this problem are the spring constant $\kappa(t)$ and swim speed $\lambda(t)$, denoted as $\mathbf{c}(t) = (\kappa(t),\lambda(t))$. These parameters are bounded as $0 \leq \lambda \leq 11$ and $1 \leq \kappa \leq 7$, reflecting typical experimental constraints.

The aim is to:
1. Transform the system from a passive equilibrium state at $(\kappa(0),\lambda(0)) = (4,0)$ to an active state at $(\kappa(t_f),\lambda(t_f)) = (4,5)$ in time $t_f=0.44$
2. Extract as much work as possible during this transformation

The order parameter to be minimized is $\Delta + 10^{-3} \langle W \rangle$, where:

$$\Delta = N^{-1} \sum_{i=1}^{N} \left[ \rho(r_i, t_f) - \rho_{\text{ss}}(r_i) \right]^2$$

measures the mean-squared error between the actual radial distribution $\rho(r,t_f)$ and the target steady-state distribution $\rho_{\text{ss}}(r)$, using $N=100$ values of $r_i$. The term $\langle W \rangle$ is the mean work done on the system (negative values indicate work extraction).

## Default Protocol

The default protocol immediately sets the control parameters to their final values. This protocol fails to effect a proper state-to-state transformation and extracts no work, as the trap stiffness $\kappa$ remains constant throughout.

## Code Structure

The code consists of:
- `engine_abp.c`: Main simulation code
- `engine_abp.h`: Header file with API definitions
- `movie.py`: Visualization script to generate snapshots and animations
- `input_control_parameters_learned.dat`: Example of an optimized protocol
- `target_p_r.dat`: Target radial distribution function

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

- `final_answer()`: Calculates the order parameter over 100 samples of 10^6 trajectories
- `load_protocol()`: Reads protocol from `input_control_parameters.dat`
- `load_default_protocol()`: Loads the default protocol
- `visualize_protocol()`: Generates visualizations of the current protocol
- `calculate_order_parameter(n_traj)`: Computes order parameter over n_traj trajectories

## Visualization

The simulation generates several output files:

- `report_position_time_*.dat`: Particle position distributions at different times
- `report_boltz_time_*.dat`: Boltzmann distributions at different times
- `report_potential_time_*.dat`: Potential energy profiles at different times
- `output_figure.pdf`: Grid of snapshots showing system evolution
- `output_movie.mp4`: Animation of the simulation

## Protocol Learning

This benchmark is designed to test methods for optimizing time-dependent protocols. The file `input_control_parameters_learned.dat` contains an example of a protocol learned via neuroevolution that successfully enacts the state-to-state transformation while extracting approximately 3 $k_BT$ units of work from the system.

The learned protocol involves an initial input of work as the trap parameter $\kappa$ is increased, followed by work extraction when $\kappa$ is returned to its resting value.

## Requirements

- C++ compiler with C++11 support
- Python 3.6+ with numpy, matplotlib, and imageio (for visualization)
- GNU Make

## References

This benchmark is based on the work described in:
- Baldovin, M., Vulpiani, A., Puglisi, A., & Baiesi, M. (2023). Control of active particles: Optimal protocols and speed limits. Physical Review E, 108(1), 014605.
- Casert, C., et al. (2024). Learning efficient protocols for active Brownian particles.
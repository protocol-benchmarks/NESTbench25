from jax import jit, grad, vmap, random, lax, ops, tree_util, value_and_grad
import jax.numpy as jnp
from jax_md import simulate
import functools
from jax.scipy import ndimage

f32 = jnp.float32
f64 = jnp.float64

def Vspring(box_size, k):
  """Return a harmonic spring energy function with periodic boundary conditions.

  Args:
    box_size: length of the periodic simulation box.
    k: spring constant (trap stiffness).

  Returns:
    spring_energy(position, r0): energy of a particle at `position` with trap center `r0`.
  """
  def spring_energy(position, r0=0.0, **kwargs):
    dR = jnp.mod((position - r0) + box_size * f32(0.5), box_size) - f32(0.5) * box_size
    return k/2 * dR ** 2
  return spring_energy

def MakeSchedule_linear(schedule, simulation_steps, r0_init, r0_final):
  """Interpolate a set of points into a full per-time-step protocol schedule.

  Linearly interpolates `schedule` (the optimizable points) across all
  simulation steps, then prepends `r0_init` and appends `r0_final` so that the
  boundary values are always exactly enforced, enabling the optimizer to find discrete jumps at the start and end of the protocol.

  Args:
    schedule: 1-D array of N values to interpolate between.
    simulation_steps: total number of simulation time steps (length of output).
    r0_init: trap position at t=0 (prepended to the array).
    r0_final: trap position at t=tf (appended to the array).

  Returns:
    1-D array of length `simulation_steps` giving the trap position at each step.
  """
  # ideally, len(schedule)+2 should divide evenly into simulation_steps
  sim_steps = jnp.arange(1, simulation_steps - 1)
  N = len(schedule) - 1
  span = sim_steps[-1] - sim_steps[0]
  inds = jnp.array([(N / span) * (sim_steps - 1)])
  # enforce r0[0] = init, r0[tf] = final
  mapped_sched = ndimage.map_coordinates(schedule, inds, 1)
  mapped_sched = jnp.concatenate([mapped_sched, jnp.reshape(r0_final, [1])])
  mapped_sched = jnp.concatenate([jnp.reshape(r0_init, [1]), mapped_sched])
  return mapped_sched

def make_trap_fxn(schedule, sim_steps, r0_init, r0_final):
  """Return a callable that maps a time-step index to the trap center position.

  Args:
    schedule: trap position values for the protocol (see MakeSchedule_linear).
    sim_steps: total number of simulation time steps.
    r0_init: initial trap position.
    r0_final: final trap position.

  Returns:
    Get_r0(step): returns the trap center at integer time step `step`.
  """
  positions = MakeSchedule_linear(schedule, sim_steps, r0_init, r0_final)
  def Get_r0(step):
    return positions[step]
  return Get_r0

def theoretical_minW_schedule(dt, simulation_steps, end_time, init, final):
  """Compute the analytically optimal (minimum-work) trap protocol.

  Implements the Schmiedl-Seifert (EPL 2007) protocol:
    lambda*(t) = (final - init) * (t + 1) / (end_time + 2) + init
  with jump discontinuities at t=0 and t=end_time enforced by the boundary values.

  Args:
    dt: integration time step.
    simulation_steps: total number of time steps.
    end_time: total protocol duration.
    init: initial trap position.
    final: final trap position.

  Returns:
    1-D array of length `simulation_steps` with the optimal trap position at each step.
  """
  times = dt*(jnp.arange(1,simulation_steps-1))
  positions = (final-init)*(times+1)/(end_time+2)+init
  positions = jnp.concatenate([positions, jnp.reshape(final, [1])])
  positions = jnp.concatenate([jnp.reshape(init, [1]), positions])
  return positions

def make_trap_fxn_theoretical(dt, simulation_steps, end_time, init, final):
  """Return a callable for the analytically optimal trap protocol.

  Args:
    dt: integration time step.
    simulation_steps: total number of time steps.
    end_time: total protocol duration.
    init: initial trap position.
    final: final trap position.

  Returns:
    Get_r0(step): returns the optimal trap center at integer time step `step`.
  """
  positions = theoretical_minW_schedule(dt, simulation_steps, end_time, init, final)
  def Get_r0(step):
      return positions[step]
  return Get_r0

def run_brownian_opt(energy_fn, schedule, init_position, r0_init, r0_final, Neq, shift, key, simulation_steps, dt=1e-5, temperature=1e-5):
  """Runs a simulation of an overdamped Brownian particle being dragged by a moving harmonic trap moving according to 'schedule', keeping track of the external work performed.
  Args:
    energy_fn: the function that governs particle energy. Here, an external harmonic potential
    schedule: the parameters of the function that yields the trap position at a given time (i.e. the control protocol parameters)
    r0_init: initial position of the trap, for which equilibration is done
    Neq: number of equilibration steps
    shift: shift_fn governing the Brownian simulation; properly accounts for periodic boundary conditions if used
    key: random key
    simulation_steps: total # simulation steps
    dt: integration time step
    temperature: simulation temperature kT

  Returns:
    the trajectory of particle positions and the cumulative work required to drag the particle, from eq'n 17 in Jarzynski 2008
  """
  trap_fn = make_trap_fxn(schedule, simulation_steps, r0_init, r0_final)
  init, apply = simulate.brownian(energy_fn, shift, dt=dt, kT=temperature, gamma=1.)
  apply = jit(apply)
  
  def equilibrate(state, Neq, apply):
    @jit
    def scan_eq(state, step):
      state = apply(state, t=step, r0=r0_init)
      return state, 0
    state, _ = lax.scan(scan_eq,state,jnp.arange(Neq))
    return state

  def increment_work(state, step):
        return (energy_fn(state.position, r0=trap_fn(step)) - energy_fn(state.position, r0=trap_fn(step-1)))

  @jit
  def scan_fn(state, step):
    dW = increment_work(state, step) #increment based on position BEFORE 'thermal kick' a la Crooks
    # Dynamically pass r0 to apply, which passes it on to energy_fn
    state = apply(state, t=step, r0=trap_fn(step))
    return state, (state.position, dW)

  key, split = random.split(key)
  eq_state = equilibrate(init(split, init_position), Neq, apply)
  state=eq_state
  state, (positions, works) = lax.scan(scan_fn,state,(jnp.arange(simulation_steps-1)+1))
  positions = jnp.concatenate([jnp.reshape(eq_state.position, [1]), positions])
  return positions, works

def single_estimate(energy_fn, init_position, r0_init, r0_final, Neq, shift, simulation_steps, dt, temperature):
  """Return a value-and-grad function for one trajectory's total work.

  Args:
    energy_fn: harmonic trap energy function.
    init_position: initial particle position.
    r0_init: initial trap center (used for equilibration).
    r0_final: final trap center.
    Neq: number of equilibration steps.
    shift: shift function for periodic boundary conditions.
    simulation_steps: total number of simulation time steps.
    dt: integration time step.
    temperature: simulation temperature (kBT).

  Returns:
    _single_estimate(schedule, seed): returns ((total_work, (positions, per_step_works)), grad)
    where grad is the gradient of total_work with respect to schedule.
  """
  @functools.partial(value_and_grad, has_aux=True)
  def _single_estimate(schedule, seed):
      positions, works = run_brownian_opt(energy_fn, schedule, init_position, r0_init, r0_final, Neq, shift, seed, simulation_steps, dt, temperature)
      total_work = works.sum()
      summary = (positions, works)
      return total_work, summary
  return _single_estimate

def estimate_gradient(batch_size, energy_fn, init_position, r0_init, r0_final, Neq, shift, simulation_steps, dt, temperature):
  """Return a JIT-compiled function that estimates the mean gradient over a batch of trajectories.

  Args:
    batch_size: number of independent trajectories to average over.
    energy_fn: harmonic trap energy function.
    init_position: initial particle position.
    r0_init: initial trap center.
    r0_final: final trap center.
    Neq: number of equilibration steps.
    shift: shift function for periodic boundary conditions.
    simulation_steps: total number of simulation time steps.
    dt: integration time step.
    temperature: simulation temperature (kBT).

  Returns:
    _estimate_gradient(schedule, seed): returns (mean_grad, (total_works, summary)) where
    total_works has shape (batch_size,) and summary is (positions, per_step_works) with
    positions of shape (batch_size, simulation_steps).
  """
  mapped_estimate = vmap(single_estimate(energy_fn, init_position, r0_init, r0_final, Neq, shift, simulation_steps, dt, temperature), [None, 0])
  @jit
  def _estimate_gradient(schedule, seed):
    seeds = random.split(seed, batch_size)
    (works, summary), grad = mapped_estimate(schedule, seeds)
    return jnp.mean(grad, axis=0), (works, summary)
  return _estimate_gradient


def plot_with_stddev(x, label=None, n=1, axis=0, ax=None, color='k', alpha=0.):
  """Plot the mean of x with a shaded ±n standard-deviation band.

  Args:
    x: 2-D array; the mean and std are taken along `axis`.
    label: legend label for the mean line.
    n: number of standard deviations for the shaded band.
    axis: axis along which to compute statistics.
    ax: matplotlib Axes to plot on.
    color: line and band color.
    alpha: transparency of the shaded band.
  """
  stddev = jnp.std(x, axis)
  mn = jnp.mean(x, axis)

  ax.fill_between(jnp.arange(mn.shape[0]),
                  mn + n * stddev, mn - n * stddev, alpha=alpha)
  ax.plot(mn, '-o', label=label,color=color)

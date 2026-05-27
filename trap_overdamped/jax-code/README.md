# JAX Implementation: Overdamped Particle in Harmonic Trap

This directory contains a JAX/JAX-MD implementation that learns the minimum-work protocol for dragging an overdamped Brownian particle with a harmonic trap, using automatic differentiation and gradient-based optimization.

## Installation

```bash
pip install -r requirements.txt
```

`ffmpeg`, required to produce the animation, is bundled via `imageio-ffmpeg` and installed automatically with the requirements above — no separate system install needed.

## Files

| File | Description |
|------|-------------|
| `optimize_overdamped_brownian.py` | Main script: optimizes the protocol and plots results |
| `utils.py` | Simulation engine and gradient estimation utilities |
| `make_animation.py` | Creates an animation of the optimal protocol dragging particle distributions and cumulative work over the course of the simulation|
| `requirements.txt` | Python package dependencies |

## Usage

```bash
python optimize_overdamped_brownian.py
```

The key simulation parameters are defined near the top of `optimize_overdamped_brownian.py`:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `end_time` | `1.0` | Protocol duration (s) |
| `dt` | `1e-3` | Integration time step |
| `batch_size` | `5000` | Trajectories per gradient estimate |
| `opt_steps` | `100` | Number of Adam optimization steps |
| `r0_init` | `5.0` | Initial trap center position |
| `r0_final` | `10.0` | Final trap center position |
| `trap_stiffness` | `1.0` | Harmonic trap spring constant |
| `savedir` | `'./'` | Directory for saved output files |

## Output

After running, the script produces:

1. **Terminal output:** the final mean dissipated work averaged over the last batch of trajectories.
2. **`schedules.pkl`** and **`works.pkl`**: the full optimization trajectory (protocol parameters and work values at each step).
3. **A two-panel plot** showing (left) convergence of mean work vs. optimization step and (right) the evolution of the learned protocol against the analytical Schmiedl-Seifert solution.
4. **`particle_animation.mp4`**: two-panel animation showing (left) an ensemble of particle positions and the moving trap potential and (right) mean cumulative work vs. time compared to the theoretical optimum.

## Optimal Protocol

The analytically optimal (minimum-mean-work) protocol is known from Schmiedl & Seifert, EPL **81**, 20003 (2008). It has a linear interior:

$$\lambda^\star(t) = (\lambda_f - \lambda_0)\,\frac{t+1}{t_f+2} + \lambda_0, \quad 0 < t < t_f$$

with jump discontinuities at the boundaries. For the default parameters ($\lambda_0 = 5$, $\lambda_f = 10$, $t_f = 1$), the minimum mean work is:

$$\langle W \rangle^\star = \frac{(\lambda_f - \lambda_0)^2}{t_f + 2} = \frac{25}{3} \approx 8.333\; k_BT$$

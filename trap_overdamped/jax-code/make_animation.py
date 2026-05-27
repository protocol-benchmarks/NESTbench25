import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import imageio_ffmpeg
plt.rcParams['animation.ffmpeg_path'] = imageio_ffmpeg.get_ffmpeg_exe()

import utils


def create_histogram_animation(
    summaries,
    schedules,
    energy_fn,
    box_size,
    simulation_steps,
    batch_size,
    dt,
    r0_init,
    r0_final,
    trap_stiffness,
    savedir,
    fps=30,
    skip_steps=10,
):
    """Create and save a two-panel animation of the final learned protocol.

    Left panel: the moving harmonic trap potential with particle positions
    scattered at their energy coordinates.
    Right panel: mean cumulative work building up over time, with a moving dot,
    compared against the Schmiedl-Seifert theoretical optimum.

    Args:
        summaries: list of (positions, per_step_works) tuples, one per optimization step.
                   positions has shape (batch_size, simulation_steps).
                   per_step_works has shape (batch_size, simulation_steps-1).
        schedules: list of (iteration, schedule_params) tuples from optimization.
        energy_fn: harmonic trap energy function (unused; accepted for API consistency).
        box_size: periodic simulation box size.
        simulation_steps: total number of simulation time steps.
        batch_size: number of trajectory samples per optimization step.
        dt: integration time step.
        r0_init: initial trap center position.
        r0_final: final trap center position.
        trap_stiffness: harmonic trap spring constant.
        savedir: directory path in which to save the animation file.
        fps: frames per second for the saved animation.
        skip_steps: number of simulation steps skipped between animation frames.

    Returns:
        matplotlib FuncAnimation object.
    """
    positions_batch, works_batch = summaries[-1]
    positions_batch = np.array(positions_batch)  # (batch_size, simulation_steps)
    works_batch = np.array(works_batch)           # (batch_size, simulation_steps-1)

    _, final_schedule = schedules[-1]
    trap_fn = utils.make_trap_fxn(final_schedule, simulation_steps, r0_init, r0_final)

    times = dt * np.arange(simulation_steps)
    frame_steps = np.arange(0, simulation_steps, skip_steps)

    # Mean cumulative work at each simulation step; prepend 0 for t=0
    mean_cumwork = np.concatenate(
        [[0.], np.cumsum(works_batch, axis=1).mean(axis=0)]
    )  # (simulation_steps,)

    # Subsample particles for scatter plot (5000 dots is too dense to render well)
    plot_n = min(batch_size, 500)
    rng = np.random.default_rng(0)
    plot_idx = rng.choice(batch_size, plot_n, replace=False)

    x_lo = r0_init - 5
    x_hi = r0_final + 5
    xs = np.linspace(x_lo, x_hi, 300)
    max_energy = 15  # kT
    theory_work = (r0_final - r0_init) ** 2 / (times[-1] + 2)

    fig, (pax, qax) = plt.subplots(1, 2, figsize=[13, 6])
    plt.rc('xtick', labelsize=20)
    plt.rc('ytick', labelsize=20)
    plt.subplots_adjust(left=0.12, right=0.97, wspace=0.35)

    # --- Left panel ---
    trap_curve, = pax.plot([], [], 'r--', lw=2, alpha=0.5, label='Trap potential')
    particles = pax.scatter([], [], c='g', marker='.', s=60, alpha=0.4, label='Particles')
    pax.set_xlim(x_lo, x_hi)
    pax.set_ylim(-0.5, max_energy)
    pax.set_xlabel("Position", fontsize=20)
    pax.set_ylabel("Energy (kT)", fontsize=20)
    pax.legend(loc='upper right', fontsize=14)
    time_text = pax.text(0.03, 0.93, '', transform=pax.transAxes, fontsize=14)

    # --- Right panel ---
    work_line, = qax.plot([], [], 'c-', lw=2, label='Mean work')
    work_dot,  = qax.plot([], [], 'co', markersize=8)
    qax.axhline(theory_work, color='tab:orange', linestyle='--', lw=2,
                label=f'Theory: {theory_work:.2f} kT')
    qax.set_xlim(times[0], times[-1])
    qax.set_ylim(-0.5, max(float(mean_cumwork[-1]), theory_work) * 1.3)
    qax.set_xlabel("Time (s)", fontsize=20)
    qax.set_ylabel("Mean cumulative work (kT)", fontsize=20)
    qax.legend(fontsize=14)

    def update(frame_idx):
        t = int(frame_steps[frame_idx])
        trap_pos = float(trap_fn(t))

        # Left panel: trap parabola + particle scatter
        trap_curve.set_data(xs, trap_stiffness / 2 * (xs - trap_pos) ** 2)
        pos_t = positions_batch[plot_idx, t]
        e_t = trap_stiffness / 2 * (pos_t - trap_pos) ** 2
        particles.set_offsets(np.column_stack([pos_t, e_t]))
        time_text.set_text(f't = {times[t]:.3f} s')

        # Right panel: cumulative mean work trace + moving dot
        work_line.set_data(times[:t + 1], mean_cumwork[:t + 1])
        work_dot.set_data([times[t]], [mean_cumwork[t]])

        return trap_curve, particles, time_text, work_line, work_dot

    anim = animation.FuncAnimation(
        fig, update, frames=len(frame_steps), interval=1000 // fps, blit=True
    )

    if savedir:
        outpath_mp4 = savedir + 'particle_animation.mp4'
        outpath_gif = savedir + 'particle_animation.gif'
        try:
            anim.save(outpath_mp4, fps=fps, extra_args=['-vcodec', 'libx264'])
            print(f'Animation saved to {outpath_mp4}')
        except Exception:
            try:
                anim.save(outpath_gif, fps=fps, writer='pillow')
                print(f'ffmpeg not found; animation saved as GIF to {outpath_gif}')
            except Exception as e:
                print(f'Could not save animation: {e}')

    plt.close(fig)
    return anim

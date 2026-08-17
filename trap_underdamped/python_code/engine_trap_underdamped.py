"""
Simulation Engine

This module implements a stochastic system with a controllable
harmonic potential. It simulates underdamped Langevin dynamics of a particle in a time-dependent
potential and calculates thermodynamic quantities like work and heat.

The system consists of a particle moving in a quadratic potential U(x,t) = 0.5*(x-λ(t))²,
where λ(t) is a time-dependent control parameter that can be optimized.

Key components:
- Protocol loading and generation
- Langevin dynamics simulation
- Work and heat calculation
- Visualization of particle distributions and potential
"""

import torch
import math
import numpy as np
import matplotlib.pyplot as plt
import os
import subprocess

# Configure matplotlib for nicer plots
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = 'Computer Modern'
plt.rcParams['text.usetex'] = True  # Enable LaTeX rendering
plt.rcParams['font.size'] = 16
plt.rcParams['xtick.labelsize'] = 16
plt.rcParams['ytick.labelsize'] = 16
plt.rcParams['legend.fontsize'] = 16
plt.rcParams['axes.titlesize'] = 16
plt.rcParams['axes.labelsize'] = 16

# Set up device and torch settings
device = "cuda" if torch.cuda.is_available() else "cpu"
torch.set_grad_enabled(False)  # Disable gradient tracking for efficiency
torch.set_default_dtype(torch.float64)  # Use double precision for numerical stability

# Simulation parameters
f_0 = 1090 #in Hz
quality = 7
omega_0 = 2*np.pi*f_0/1e6 #in recicropoal microseconds
cycle_time = 1e6/f_0 #in microseconds

trajectory_time = .15*cycle_time  # Total simulation time
timestep = 0.01  # Integration time step
number_of_steps = int(trajectory_time/timestep)  # Total number of simulation steps

# Control parameter configuration
number_of_control_parameters = 1  # Dimension of the control parameter space
cee_boundary = torch.zeros(2,number_of_control_parameters)  # Boundary conditions for control parameters
cee_boundary[0,:] = 0.0  # Initial value (at t=0)
cee_boundary[1,:] = 5.0  # Final value (at t=T)

vee_boundary = torch.zeros(2)
vee_boundary[0] = cee_boundary[1][0]/(2.0/(quality*omega_0)+trajectory_time)
vee_boundary[1] = -cee_boundary[1][0]/(2.0/(quality*omega_0)+trajectory_time)

# Visualization parameters
number_of_pictures = 100  # Number of snapshots to save during visualization
n_bins = 50  # Number of histogram bins for distribution visualization

def load_protocol(filename="input_control_parameters_learned.dat"):
    """
    Load a control protocol from a file with automatic repetition handling.

    The protocol defines how the control parameter λ(t) changes over time.
    If the protocol has fewer points than number_of_steps, each protocol value
    is repeated for multiple consecutive timesteps. If the protocol has more
    points than number_of_steps, intermediate steps are skipped periodically.

    Parameters
    ----------
    filename : str
        Path to the file containing the protocol data

    Returns
    -------
    torch.Tensor
        A tensor of shape (number_of_steps+2, number_of_control_parameters)
        containing the protocol values at each timestep.
        The first and last entries are set to the boundary conditions.

    Raises
    ------
    FileNotFoundError
        If the specified file does not exist
    ValueError
        If the protocol generation fails or file is empty
    """
    # Initialize protocol tensor with zeros
    cee = torch.zeros(number_of_steps+2, number_of_control_parameters, device=device)

    # Load protocol values from file with automatic repetition handling
    try:
        protocol_data = np.loadtxt(filename)

        if len(protocol_data) == number_of_steps:
            # Perfect match - use as is
            protocol_values = protocol_data
        elif len(protocol_data) < number_of_steps:
            # Fewer protocol points than steps - distribute evenly with minimal remaining steps
            # Use linear interpolation to create evenly spaced protocol values
            indices = np.linspace(0, len(protocol_data) - 1, number_of_steps)
            protocol_values = np.interp(indices, np.arange(len(protocol_data)), protocol_data)
        else:
            # More protocol points than steps - skip intermediate steps periodically
            # Use linear interpolation indices to select evenly spaced points
            indices = np.linspace(0, len(protocol_data) - 1, number_of_steps, dtype=int)
            protocol_values = protocol_data[indices]

        # Verify final length matches requirements
        if len(protocol_values) != number_of_steps:
            raise ValueError(f"Failed to generate protocol with {number_of_steps} steps from {len(protocol_data)} protocol timepoints")

        cee[1:-1, 0] = torch.from_numpy(protocol_values)

    except FileNotFoundError:
        raise FileNotFoundError(f"Protocol file '{filename}' not found")
    except Exception as e:
        raise ValueError(f"Error loading protocol from '{filename}': {str(e)}")

    # Set boundary conditions
    cee[0, :] = cee_boundary[0, :]  # Initial value
    cee[-1, :] = cee_boundary[1, :]  # Final value

    return cee

def load_default_protocol():
    """
    Generate a default protocol with modified ramp behavior.

    Creates a protocol where the control parameter λ(t) follows a modified ramp. The formula creates
    discontinuities (jumps) at the beginning and end of the protocol, which is optimal for this physical system.


    Returns
    -------
    torch.Tensor
        A tensor of shape (number_of_steps+2, number_of_control_parameters)
        containing the protocol values at each timestep.
    """
    # Initialize protocol tensor with zeros
    cee = torch.zeros(number_of_steps+2, number_of_control_parameters, device=device)

    # Create modified ramp with jumps at the beginning and end
    t_values = np.linspace(0, trajectory_time, number_of_steps)
    cee[1:-1, 0] = cee_boundary[1, 0] * (quality*omega_0*t_values + 1) / (quality*omega_0*trajectory_time + 2)

    # Set boundary conditions
    cee[0, :] = cee_boundary[0, :]  # Initial value
    cee[-1, :] = cee_boundary[1, :]  # Final value

    return cee



def run_protocol(protocol, number_of_trajectories, kick_velocity= False, visualize=False):
    """
    Run a stochastic thermodynamic simulation with the given protocol.

    Simulates the Langevin dynamics of a particle in a time-dependent
    harmonic potential U(x,t) = 0.5*(x-λ(t))², where λ(t) is defined
    by the protocol. Calculates work and heat for each trajectory.

    Parameters
    ----------
    protocol : torch.Tensor
        Control parameter values at each timestep,
        shape (number_of_steps+2, number_of_control_parameters)
    number_of_trajectories : int
        Number of independent particle trajectories to simulate
    visualize : bool, optional
        Whether to collect data for visualization (default: False)

    Returns
    -------
    work : torch.Tensor
        Work done on each particle, shape (number_of_trajectories,)
    heat : torch.Tensor
        Heat exchanged for each particle, shape (number_of_trajectories,)
    pos_to_visualize : torch.Tensor, optional
        Particle positions at visualization timesteps,
        shape (number_of_pictures, number_of_trajectories)
    parameters_to_visualize : torch.Tensor, optional
        Protocol values at visualization timesteps,
        shape (number_of_pictures, number_of_control_parameters)
    """
    # Initialize visualization arrays if needed

    # Precompute constants for efficiency
    lambd = np.exp(-timestep*omega_0/quality)

    # Initialize state variables
    # Initial positions drawn from standard normal distribution (equilibrium at λ=0)
    positions = torch.randn(number_of_trajectories, device=device)
    velocities = omega_0*torch.randn(number_of_trajectories, device=device)

    heat = torch.zeros(number_of_trajectories, device=device)
    work = torch.zeros(number_of_trajectories, device=device)

    if kick_velocity:
        work+=vee_boundary[0]*(velocities+0.5*vee_boundary[0])/(omega_0*omega_0)
        velocities+=vee_boundary[0]

    # Move protocol to the correct device (GPU if available)
    protocol = protocol.to(device)

    # Preallocate noise tensor for efficiency
    noise = torch.zeros_like(positions)

    if visualize:
        pos_to_visualize = torch.zeros(number_of_pictures+2, number_of_trajectories, device=device)
        parameters_to_visualize = torch.zeros(number_of_pictures+2, number_of_control_parameters, device=device)
        pos_to_visualize[0] = positions
        parameters_to_visualize[0] = protocol[0]

    # Main simulation loop
    for i in range(1, number_of_steps - 1):
        # Calculate work done by changing the protocol
        # Work is the change in potential energy due to protocol change
        # dW = U(x,λ_new) - U(x,λ_old) = 0.5*[(x-λ_new)² - (x-λ_old)²]
        work += 0.5 * ((positions - protocol[i])**2 - (positions - protocol[i-1])**2)

        # Current energy before position update
        old_energy = 0.5 * (positions - protocol[i])**2

        # Update positions using Langevin dynamics
        # dx = -∇U(x,λ)dt + √(2dt)·η, where η is Gaussian white noise
        noise.normal_()  # Generate standard normal noise
        grad = (positions-protocol[i])*quality*omega_0
        positions += timestep * velocities
        velocities *= lambd
        velocities -= (1-lambd)*grad -np.sqrt(1-lambd*lambd)*omega_0*noise

        # Calculate heat exchange (energy change due to position update)
        # dQ = U(x_new,λ) - U(x_old,λ)
        new_energy = 0.5 * (positions - protocol[i])**2
        heat += new_energy - old_energy

        # Store data for visualization at regular intervals
        if visualize and (i-1) % (number_of_steps // number_of_pictures) == 0:
            idx = 1 + (i-1) // (number_of_steps // number_of_pictures)
            pos_to_visualize[idx, :] = positions
            parameters_to_visualize[idx, :] = protocol[i]

    # Final work calculation for the last step

    if kick_velocity:
        work+=vee_boundary[1]*(velocities+0.5*vee_boundary[1])/(omega_0*omega_0)
        velocities+=vee_boundary[1]

    work += 0.5 * ((positions - protocol[-1])**2 - (positions - protocol[-2])**2)

    # Return appropriate results based on visualization flag
    if visualize:
        pos_to_visualize[-1] = positions
        parameters_to_visualize[-1] = protocol[-1]
        return work, heat, pos_to_visualize, parameters_to_visualize
    else:
        return work, heat


def calculate_order_parameter(protocol, number_of_trajectories):
    work, heat = run_protocol(protocol, number_of_trajectories)
    return torch.mean(work).item()

def visualize_protocol(protocol, number_of_trajectories, kick_velocity= False):
    """
    Visualize the evolution of particle distributions during protocol execution.

    Creates a series of plots showing:
    1. The potential energy landscape at different timesteps
    2. The theoretical equilibrium distribution (Boltzmann)
    3. The actual particle distribution from simulation
    4. A histogram of work values at the end

    Parameters
    ----------
    protocol : torch.Tensor
        Control parameter values at each timestep
    number_of_trajectories : int
        Number of particle trajectories to simulate

    Returns
    -------
    None
        Plots are saved to the 'plots' directory
    """
    # Create plots directory if it doesn't exist
    os.makedirs('plots', exist_ok=True)
    fig, ax = plt.subplots(figsize=(6, 6))
    plt.plot(np.linspace(0,trajectory_time, number_of_steps+2), protocol[:,0], color = "g")
    plt.xlabel(r"$t$")
    plt.ylabel(r"$c_0$")
    fig.tight_layout()
    fig.savefig(f'protocol.png', pad_inches=0.05, dpi=300)


    # Run the simulation with visualization enabled
    work, heat, pos_to_visualize, parameters_to_visualize = run_protocol(
        protocol, number_of_trajectories, visualize=True, kick_velocity = kick_velocity
    )

    # Define x-axis range for plotting potential and distributions
    # Range is set to be twice the maximum protocol value in both directions
    x_min = -2.0 * cee_boundary[1, 0]
    x_max = 2.0 * cee_boundary[1, 0]
    n_x_points = 1000  # Number of points for smooth curve plotting
    x_points = torch.linspace(x_min, x_max, n_x_points)

    # Select 6 time points for the 2x3 panel figure
    selected_indices = [0, (number_of_pictures+2)//5, 2*(number_of_pictures+2)//5,
                       3*(number_of_pictures+2)//5, 4*(number_of_pictures+2)//5, number_of_pictures+1]

    # Create 2x3 panel figure
    panel_fig, panel_axes = plt.subplots(2, 3, figsize=(12, 9))
    panel_axes = panel_axes.flatten()
    panel_counter = 0

    # Generate a series of plots showing the evolution of the system
    for i in range(number_of_pictures+2):
        fig, ax = plt.subplots(figsize=(6, 6))

        # Create histogram of particle positions
        hist, bin_edges = np.histogram(pos_to_visualize[i], bins=n_bins, density=True)

        # Compute bin centers for plotting histogram
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

        # Calculate potential energy landscape: U(x,t) = 0.5*(x-λ(t))²
        potential = 0.5 * (x_points - parameters_to_visualize[i, 0])**2

        # Calculate equilibrium Boltzmann distribution: ρ₀(x) ∝ exp(-U(x,t))
        # First compute the partition function (normalization factor)
        partition_function = torch.trapz(torch.exp(-potential), x=x_points)
        # Then the normalized Boltzmann distribution
        boltzmann = torch.exp(-potential) / partition_function

        # Plot the potential energy, equilibrium distribution, and actual distribution
        ax.plot(x_points, potential, 'k', label=r"$U(x,t)$")
        ax.plot(x_points, boltzmann*20, 'k--', label=r"$\rho_0(x)$")  # Scaled for visibility
        ax.plot(bin_centers, hist*20, 'g-', label=r"$\rho(x,t)$")  # Scaled for visibility

        # Set axis limits for consistent visualization
        ax.set_ylim(0, 10)
        ax.set_xlim(x_min, x_max)

        # Customize plot appearance
        ax.set_xticks([-15, 0, 15])
        ax.set_xlabel(r"$x$", fontsize=16)
        ax.set_yticks([])  # Remove vertical axis values for cleaner look
        plt.legend(loc='lower center', bbox_to_anchor=(0.5, 1.0), ncol=3, frameon=False)

        # Save the figure as a high-resolution PNG
        fig.tight_layout()
        fig.savefig(f'plots/output_{i}.png', pad_inches=0.05, dpi=300)
        plt.close(fig)

        # Add to 2x3 panel figure if this is a selected time point
        if i in selected_indices and panel_counter < 6:
            ax_panel = panel_axes[panel_counter]

            # Plot the same data on the panel
            ax_panel.plot(x_points, potential, 'k', label=r"$U(x,t)$")
            ax_panel.plot(x_points, boltzmann*20, 'k--', label=r"$\rho_0(x)$")
            ax_panel.plot(bin_centers, hist*20, 'g-', label=r"$\rho(x,t)$")

            # Set axis limits and appearance
            ax_panel.set_ylim(0, 10)
            ax_panel.set_xlim(x_min, x_max)
            ax_panel.set_xticks([-15, 0, 15])
            ax_panel.set_xlabel(r"$x$", fontsize=14)
            ax_panel.set_yticks([])

            # Add title showing time point
            time_fraction = i / (number_of_pictures + 1)
            ax_panel.set_title(r"$t=$ "+ f"{time_fraction:.2f}"+ r"$t_{\rm{f}}$", fontsize=14)

            if panel_counter == 0:  # Add legend only to first panel
                ax_panel.legend(loc='lower right', fontsize=12, frameon=False)

            panel_counter += 1

    # Save the 2x3 panel figure
    panel_fig.tight_layout()
    panel_fig.savefig("time_evolution.png", pad_inches=0.1, dpi=300)
    plt.close(panel_fig)

    # Create movie from output images using ffmpeg
    try:
        cmd = [
            'ffmpeg', '-y',  # -y to overwrite existing file
            '-framerate', '10',  # 10 frames per second
            '-i', 'plots/output_%d.png',  # Input pattern
            '-c:v', 'libx264',  # Video codec
            '-pix_fmt', 'yuv420p',  # Pixel format for compatibility
            '-crf', '18',  # High quality
            'output.mp4'  # Output file
        ]

        result = subprocess.run(cmd, capture_output=True, text=True, cwd='.')

        if result.returncode == 0:
            print("Movie created successfully: plots/output.mp4")
        else:
            print(f"Error creating movie: {result.stderr}")

    except FileNotFoundError:
        print("ffmpeg not found. Please install ffmpeg to create movies.")

    # Create a histogram of work values
    fig, ax = plt.subplots(figsize=(6, 6))

    # Calculate histogram of work values
    hist, bin_edges = np.histogram(work, bins=n_bins, density=True)

    # Compute bin centers for plotting
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

    # Plot work distribution
    plt.plot(bin_centers, hist, color="g")
    plt.axvline(x = torch.mean(work), ls = "--", color = "g")
    plt.text(s = r"$\langle W \rangle = $" + f"{torch.mean(work).item():.3f}", x = 1.1*torch.mean(work), y = np.max(hist), ha = "left")
    plt.xlabel(r"$W$")
    plt.ylabel(r"$P(W)$")
    plt.ylim([0, 1.1*np.max(hist)])
    plt.xlim([np.min(bin_centers), np.max(bin_centers)])
    # Save work histogram
    fig.tight_layout()
    fig.savefig("work_histogram.png", pad_inches=0.05, dpi=300)
    plt.close(fig)



def final_answer(protocol, kick_velocity= False):
    """
    Calculate the statistical properties of work with high precision.

    Runs a large number of simulations (10^6 trajectories) to accurately
    estimate the mean work and its standard error. Uses batch processing
    to improve statistical analysis.

    Parameters
    ----------
    protocol : torch.Tensor
        Control parameter values at each timestep

    Returns
    -------
    None
        Prints the mean work value and standard error to the console
    """
    # Run a large number of trajectories for statistical significance
    work, heat = run_protocol(protocol, int(1e6), kick_velocity = kick_velocity)
    # Reshape work array into 100 batches of 10^4 trajectories each
    # This allows for better statistical analysis through batch means
    work = torch.reshape(work, (int(1e2), int(1e4)))

    # Calculate mean work for each batch
    sample_means = torch.mean(work, axis=1)

    # Calculate overall mean (mean of the batch means)
    overall_mean = torch.mean(sample_means)

    # Calculate standard deviation of the batch means
    # unbiased=True uses Bessel's correction (n-1 denominator)
    std_of_means = torch.std(sample_means, unbiased=True)

    # Calculate standard error of the mean
    # SE = σ/√n where σ is the standard deviation and n is the number of batches
    standard_error = std_of_means / np.sqrt(100)

    # Print results with 6 decimal places of precision
    print(f"order parameter = {overall_mean:.6f} ± {standard_error:.6f}, n_samples = {100}")

if __name__ == "__main__":
    # protocol = load_protocol()
    # visualize_protocol(protocol, int(1e5),  kick_velocity= False)
    # final_answer(protocol, kick_velocity= False)
    protocol = load_default_protocol()
    visualize_protocol(protocol, int(1e5),  kick_velocity= True)
    final_answer(protocol, kick_velocity= True)

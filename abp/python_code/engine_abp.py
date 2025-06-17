"""
2D Active Brownian Particle Simulation Engine

This module implements a 2D active Brownian particle system with controllable
activity and confinement. It simulates the dynamics of self-propelled particles
transitioning from passive to active phases while minimizing work and matching
a target radial distribution.

The system consists of particles with positions (x,y) and orientations θ evolving according to:
- dx/dt = -κx + v₀cos(θ) + √(2D)ξₓ(t)
- dy/dt = -κy + v₀sin(θ) + √(2D)ξᵧ(t)
- dθ/dt = √(2Dᵣ)ξθ(t)

where κ is the confinement strength, v₀ is the self-propulsion velocity,
D is the translational diffusion coefficient, and Dᵣ is the rotational diffusion coefficient.

Key components:
- Protocol loading and generation for activity control
- Active Brownian particle dynamics simulation
- Work calculation during activity transitions
- Radial distribution analysis and target matching
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
trajectory_time = .44  # Total simulation time
timestep = 0.001  # Integration time step
number_of_steps = int(trajectory_time/timestep)  # Total number of simulation steps
# protocol_update_time=trajectory_time/(1.0*number_of_steps)

# Control parameter configuration
number_of_control_parameters = 2  # Dimension of the control parameter space
cee_boundary = torch.zeros(2,number_of_control_parameters)  # Boundary conditions for control parameters
cee_boundary[:,0] = 4.0
cee_boundary[:,1] = torch.Tensor([0.0, 5.0])

# Visualization parameters
number_of_pictures = 100  # Number of snapshots to save during visualization
n_bins = 51  # Number of histogram bins for distribution visualization

def load_protocol(filename="input_control_parameters.dat"):
    """
    Load a control protocol from a file for active Brownian particle transitions.

    The protocol defines how the control parameters [κ(t), v₀(t)] change over time,
    controlling the transition from passive (v₀=0) to active (v₀>0) phases.
    If the protocol has fewer points than number_of_steps, linear interpolation
    is used to generate the required number of timesteps.

    Parameters
    ----------
    filename : str
        Path to the file containing the protocol data.
        Each row should contain [κ(t), v₀(t)] values.

    Returns
    -------
    torch.Tensor
        A tensor of shape (number_of_steps+2, number_of_control_parameters)
        containing the protocol values at each timestep.
        Column 0: confinement strength κ(t)
        Column 1: self-propulsion velocity v₀(t)

    Raises
    ------
    FileNotFoundError
        If the specified file does not exist
    ValueError
        If the protocol generation fails or file format is incorrect
    """
    # Initialize protocol tensor with zeros
    cee = torch.zeros(number_of_steps+2, number_of_control_parameters, device=device)

    # Load protocol values from file with automatic repetition handling
    try:
        protocol_data = np.loadtxt(filename)

        # Ensure protocol_data is 2D even if it's a single column
        if protocol_data.ndim == 1:
            protocol_data = protocol_data.reshape(-1, 1)

        # Check if the number of columns matches the expected control parameters
        if protocol_data.shape[1] != number_of_control_parameters:
            raise ValueError(f"Protocol file has {protocol_data.shape[1]} columns but expected {number_of_control_parameters} control parameters")

        if len(protocol_data) == number_of_steps:
            # Perfect match - use as is
            protocol_values = protocol_data
        elif len(protocol_data) < number_of_steps:
            # Fewer protocol points than steps - distribute evenly with minimal remaining steps
            # Use linear interpolation to create evenly spaced protocol values
            indices = np.linspace(0, len(protocol_data) - 1, number_of_steps)
            protocol_values = np.zeros((number_of_steps, number_of_control_parameters))
            for param_idx in range(number_of_control_parameters):
                protocol_values[:, param_idx] = np.interp(indices, np.arange(len(protocol_data)), protocol_data[:, param_idx])
        else:
            # More protocol points than steps - skip intermediate steps periodically
            # Use linear interpolation indices to select evenly spaced points
            indices = np.linspace(0, len(protocol_data) - 1, number_of_steps, dtype=int)
            protocol_values = protocol_data[indices]

        # Verify final shape matches requirements
        if protocol_values.shape != (number_of_steps, number_of_control_parameters):
            raise ValueError(f"Failed to generate protocol with shape ({number_of_steps}, {number_of_control_parameters}) from protocol data with shape {protocol_data.shape}")

        cee[1:-1] = torch.from_numpy(protocol_values)

    except FileNotFoundError:
        raise FileNotFoundError(f"Protocol file '{filename}' not found")
    except Exception as e:
        raise ValueError(f"Error loading protocol from '{filename}': {str(e)}")

    # Set boundary conditions
    cee[0, :] = cee_boundary[0, :]  # Initial value
    cee[-1, :] = cee_boundary[-1,:]  # Final value

    return cee

def load_default_protocol():
    """
    Generate a default protocol for passive-to-active transition.

    Creates a simple protocol where particles transition from passive (v₀=0)
    to active (v₀=5) behavior while maintaining constant confinement strength.

    Returns
    -------
    torch.Tensor
        A tensor of shape (number_of_steps+2, number_of_control_parameters)
        containing the protocol values at each timestep.
    """
    # Initialize protocol tensor with zeros
    cee = torch.zeros(number_of_steps+2, number_of_control_parameters, device=device)


    cee[1:-1] = cee_boundary[-1]
    # Set boundary conditions
    cee[0, :] = cee_boundary[0, :]  # Initial value
    cee[-1, :] = cee_boundary[-1,:]  # Final value

    return cee



def run_protocol(protocol, number_of_trajectories, visualize=False):
    """
    Run a 2D active Brownian particle simulation with the given protocol.

    Simulates the dynamics of self-propelled particles with time-dependent
    confinement κ(t) and self-propulsion velocity v₀(t). Particles evolve
    according to overdamped Langevin equations with orientational dynamics.

    Parameters
    ----------
    protocol : torch.Tensor
        Control parameter values at each timestep,
        shape (number_of_steps+2, number_of_control_parameters)
        Column 0: confinement strength κ(t)
        Column 1: self-propulsion velocity v₀(t)
    number_of_trajectories : int
        Number of independent particle trajectories to simulate
    visualize : bool, optional
        Whether to collect data for visualization (default: False)

    Returns
    -------
    work : torch.Tensor
        Work done on each particle during the protocol, shape (number_of_trajectories,)
    heat : torch.Tensor
        Heat exchanged for each particle, shape (number_of_trajectories,)
    positions : torch.Tensor
        Final particle positions, shape (number_of_trajectories, 2)
    pos_to_visualize : torch.Tensor, optional (if visualize=True)
        Particle positions at visualization timesteps,
        shape (number_of_pictures+2, number_of_trajectories, 2)
    parameters_to_visualize : torch.Tensor, optional (if visualize=True)
        Protocol values at visualization timesteps,
        shape (number_of_pictures+2, number_of_control_parameters)
    """
    # Initialize visualization arrays if needed

    # Precompute constants for efficiency
    sqrt2dt = math.sqrt(2 * timestep)  # Noise scaling factor for Langevin dynamics

    # Initialize state variables
    # Initial positions drawn from standard normal distribution (equilibrium at λ=0)
    positions = torch.randn(number_of_trajectories, 2,device=device)/torch.sqrt(cee_boundary[0,0])
    angles = 2*np.pi*torch.rand(number_of_trajectories, device = device)
    heat = torch.zeros(number_of_trajectories, device=device)
    work = torch.zeros(number_of_trajectories, device=device)

    # Move protocol to the correct device (GPU if available)
    protocol = protocol.to(device)

    # Preallocate noise tensor for efficiency
    noise = torch.zeros(number_of_trajectories,device=device)

    if visualize:
        pos_to_visualize = torch.zeros(number_of_pictures+2, number_of_trajectories, 2, device=device)
        parameters_to_visualize = torch.zeros(number_of_pictures+2, number_of_control_parameters, device=device)
        pos_to_visualize[0] = positions
        parameters_to_visualize[0] = protocol[0]

    # Main simulation loop
    for i in range(1, number_of_steps - 1):
        # Calculate work done by changing the protocol
        # Work is the change in potential energy due to protocol change
        # dW = U(x,λ_new) - U(x,λ_old) = 0.5*[(x-λ_new)² - (x-λ_old)²]
        work += .5*(protocol[i,0] - protocol[i-1,0])*(positions**2).sum(dim = 1)

        # Current energy before position update
        old_energy =0.5 * protocol[i,0]*(positions**2).sum(dim=1)

        # Update positions using Langevin dynamics
        # dx = -∇U(x,λ)dt + √(2dt)·η, where η is Gaussian white noise
        noise.normal_(std = sqrt2dt)  # Generate standard normal noise
        positions[:,0] += (protocol[i,1]*torch.cos(angles) - protocol[i,0] * positions[:,0])*timestep + noise
        noise.normal_(std = sqrt2dt)
        positions[:,1] += (protocol[i,1]*torch.sin(angles) - protocol[i,0] * positions[:,1])*timestep + noise
        noise.normal_(std = sqrt2dt)
        angles += noise


        # Calculate heat exchange (energy change due to position update)
        # dQ = U(x_new,λ) - U(x_old,λ)
        new_energy = 0.5 * protocol[i,0]*(positions**2).sum(dim=1)
        heat += new_energy - old_energy

        # Store data for visualization at regular intervals
        if visualize and (i-1) % (number_of_steps // number_of_pictures) == 0:
            idx = 1 + (i-1) // ((number_of_steps)  // number_of_pictures)
            idx = min(idx, number_of_pictures)  # Clamp idx to valid range
            pos_to_visualize[idx, :] = positions
            parameters_to_visualize[idx, :] = protocol[i]

    # Final work calculation for the last step
    work += .5*(protocol[-1,0] - protocol[-2,0])*(positions**2).sum(dim = 1)
    # Return appropriate results based on visualization flag
    if visualize:
        pos_to_visualize[-1] = positions
        parameters_to_visualize[-1] = protocol[-1]
        return work, heat, positions, pos_to_visualize, parameters_to_visualize
    else:
        return work, heat, positions

def visualize_protocol(protocol, number_of_trajectories):
    """
    Visualize the evolution of active Brownian particle distributions during protocol execution.

    Creates a series of plots showing:
    1. The confinement potential energy landscape U(r) = 0.5*κ*r² at different timesteps
    2. The theoretical equilibrium distribution for passive particles (Boltzmann)
    3. The actual particle radial distribution from active Brownian dynamics
    4. The target active phase distribution
    5. A histogram of work values at the end
    6. Time evolution of control parameters [κ(t), v₀(t)]

    Parameters
    ----------
    protocol : torch.Tensor
        Control parameter values at each timestep
        Column 0: confinement strength κ(t)
        Column 1: self-propulsion velocity v₀(t)
    number_of_trajectories : int
        Number of particle trajectories to simulate

    Returns
    -------
    None
        Plots are saved to the 'plots' directory and a movie is generated
    """
    # Create plots directory if it doesn't exist
    os.makedirs('plots', exist_ok=True)
    fig, ax = plt.subplots(figsize=(6, 6))
    plt.plot(np.linspace(0,trajectory_time, number_of_steps+2), protocol[:,0], color = "g")
    plt.plot(np.linspace(0,trajectory_time, number_of_steps+2), protocol[:,1], color = "b")
    plt.xlabel(r"$t$")
    plt.ylabel(r"$c$")
    fig.tight_layout()
    fig.savefig(f'protocol.png', pad_inches=0.05, dpi=300)


    # Run the simulation with visualization enabled
    work, heat,positions, pos_to_visualize, parameters_to_visualize = run_protocol(
        protocol, number_of_trajectories, visualize=True
    )
    # Define x-axis range for plotting potential and distributions
    # Range is set to be twice the maximum protocol value in both directions
    x_min = 1e-3
    x_max = 5/np.sqrt(2)
    n_x_points = 1000  # Number of points for smooth curve plotting
    x_points = torch.linspace(x_min, x_max, n_x_points)

    # Select 6 time points for the 2x3 panel figure
    selected_indices = [0, (number_of_pictures+2)//5, 2*(number_of_pictures+2)//5,
                       3*(number_of_pictures+2)//5, 4*(number_of_pictures+2)//5, number_of_pictures+1]

    # Create 2x3 panel figure
    panel_fig, panel_axes = plt.subplots(2, 3, figsize=(12, 9))
    panel_axes = panel_axes.flatten()
    panel_counter = 0

    target = np.loadtxt("target_p_r.dat")

    # Generate a series of plots showing the evolution of the system
    for i in range(number_of_pictures+2):
        fig, ax = plt.subplots(figsize=(7.5, 6))

        r = torch.sqrt((pos_to_visualize[i]**2).sum(dim=1))
        # Create histogram of particle positions
        hist, bin_edges = np.histogram(r, bins=n_bins, density=True)

        # Compute bin centers for plotting histogram
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

        # Calculate potential energy landscape: U(x,t) = 0.5*(x-λ(t))²
        pot = .5 * parameters_to_visualize[i,0] *x_points**2
        # Calculate equilibrium Boltzmann distribution: ρ₀(x) ∝ exp(-U(x,t))
        # First compute the partition function (normalization factor)
        partition_function = torch.trapz(x_points*torch.exp(-pot), x=x_points)
        # Then the normalized Boltzmann distribution
        boltzmann = x_points*torch.exp(-pot) / partition_function

        # Plot the potential energy, equilibrium distribution, and actual distribution
        ax.plot(x_points, pot, 'k', label=r"$U(r,t)$")
        ax.plot(x_points, boltzmann*20, 'k--', label=r"$r\rho_0(r)$")  # Scaled for visibility
        ax.plot(bin_centers, hist*20, 'g-', label=r"$r\rho(r,t)$")  # Scaled for visibility
        ax.plot(target[:,0], target[:,1]*20, color = "k", ls = ":", label = "Target")
        # Set axis limits for consistent visualization
        ax.set_ylim(0, 35)
        ax.set_xlim(0, x_max)

        # Customize plot appearance
        ax.set_xticks([0, 3])
        ax.set_xlabel(r"$r$", fontsize=16)
        ax.set_yticks([])  # Remove vertical axis values for cleaner look
        plt.legend(loc='lower center', bbox_to_anchor=(0.5, 1.0), ncol=4, frameon=False)

        # Save the figure as a high-resolution PNG
        fig.tight_layout()
        fig.savefig(f'plots/output_{i}.png', pad_inches=0.05, dpi=300)
        plt.close(fig)

        # Add to 2x3 panel figure if this is a selected time point
        if i in selected_indices and panel_counter < 6:
            ax_panel = panel_axes[panel_counter]

            # Plot the same data on the panel
            ax_panel.plot(x_points, pot, 'k', label=r"$U(r,t)$")
            ax_panel.plot(x_points, boltzmann*20, 'k--', label=r"$r\rho_0(r)$")
            ax_panel.plot(bin_centers, hist*20, 'g-', label=r"$r\rho(r,t)$")
            ax_panel.plot(target[:,0], target[:,1]*20, color = "k", ls = ":", label = "Target")

            # Set axis limits and appearance
            ax_panel.set_ylim(0, 35)
            ax_panel.set_xlim(0, x_max)
            ax_panel.set_xticks([ 0, 3])
            ax_panel.set_xlabel(r"$r$", fontsize=14)
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

def calculate_order_parameter(protocol, number_of_trajectories=int(1e5)):
    """
    Calculate the order parameter for active Brownian particle transition.

    The order parameter combines work minimization with radial distribution matching.
    It evaluates how well the protocol achieves the desired active phase distribution
    while minimizing the work cost of the transition.

    Parameters
    ----------
    protocol : torch.Tensor
        Control parameter values at each timestep
    number_of_trajectories : int, optional
        Number of trajectories to simulate (default: 1e5)

    Returns
    -------
    float
        The calculated order parameter combining:
        - Work penalty: scaled mean work done during transition
        - Distribution error: MSE between final radial distribution and target
    """
    work, heat, positions = run_protocol(protocol, number_of_trajectories)
    rs = np.sqrt(positions[:,0]**2 + positions[:,1]**2)
    target = np.loadtxt("target_p_r.dat")
    bin_centers = target[:,0]

    bin_width = bin_centers[1] - bin_centers[0]
    bin_edges = np.concatenate([
        [bin_centers[0] - bin_width/2],  # left edge of first bin
        bin_centers[:-1] + bin_width/2,  # edges between bins
        [bin_centers[-1] + bin_width/2]  # right edge of last bin
    ])

    # Create histogram
    hist_counts, _ = np.histogram(rs, bins=bin_edges, density = True)

    return 1e-3*torch.mean(work).item() + np.mean((hist_counts - target[:,1])**2)


def final_answer(protocol):
    """
    Calculate the final performance metrics for the active Brownian particle protocol.

    Runs a large-scale simulation to accurately estimate the order parameter,
    which quantifies how well the protocol achieves the passive-to-active
    transition while minimizing work and matching the target distribution.

    Parameters
    ----------
    protocol : torch.Tensor
        Control parameter values at each timestep

    Returns
    -------
    None
        Prints the order parameter value to the console
    """
    # Run a large number of trajectories for statistical significance
    # Print results with 6 decimal places of precision
    print(f"order parameter = {calculate_order_parameter(protocol, number_of_trajectories=int(1e6))}")

if __name__ == "__main__":
    protocol = load_default_protocol()
    visualize_protocol(protocol, int(1e5))
    final_answer(protocol)

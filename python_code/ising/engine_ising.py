"""
2D Ising Model Stochastic Thermodynamics Simulation Engine

This module implements a 2D Ising model with controllable temperature and magnetic field
for studying nonequilibrium thermodynamics. It simulates Monte Carlo dynamics of magnetic
spins during field reversal protocols while calculating entropy production.

The system consists of spins σᵢ ∈ {-1, +1} on a square lattice evolving under the Hamiltonian:
H = -J∑⟨i,j⟩σᵢσⱼ - h∑ᵢσᵢ

where J is the coupling strength (set to 1), h(t) is the time-dependent magnetic field,
and the dynamics occur at temperature T(t).

Key components:
- Protocol loading and generation for temperature and magnetic field control
- Monte Carlo dynamics simulation using Metropolis algorithm
- Entropy production calculation during field reversal
- Visualization of spin configurations and thermodynamic trajectories
- Statistical analysis of entropy production distributions
- Checkerboard updating for efficient parallel simulation
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
L = 32  # Linear lattice size (L×L spins total)
number_of_steps = 100  # Total number of Monte Carlo sweep steps

# Control parameter configuration
number_of_control_parameters = 2  # [Temperature T, Magnetic field h]
cee_boundary = torch.zeros(2, number_of_control_parameters)  # Boundary conditions
cee_boundary[:, 0] = 0.7  # Temperature (kept constant at T ≈ 0.7 > Tc)
cee_boundary[:, 1] = torch.Tensor([-1., 1.])  # Magnetic field (flip from -1 to +1)

# Visualization parameters
number_of_pictures = 100  # Number of snapshots to save during visualization
n_bins = 51  # Number of histogram bins for distribution visualization

def load_protocol(filename="input_control_parameters.dat"):
    """
    Load a control protocol from a file for Ising model field reversal.

    The protocol defines how the control parameters [T(t), h(t)] change over time,
    typically keeping temperature constant while reversing the magnetic field.
    Linear interpolation is used if the protocol has different length than number_of_steps.

    Parameters
    ----------
    filename : str
        Path to the file containing the protocol data.
        Each row should contain [T(t), h(t)] values.

    Returns
    -------
    torch.Tensor
        A tensor of shape (number_of_steps+2, number_of_control_parameters)
        containing the protocol values at each timestep.
        Column 0: Temperature T(t)
        Column 1: Magnetic field h(t)

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
    Generate a default protocol for magnetic field reversal.

    Loads a pre-computed optimal protocol from file that minimizes entropy
    production during the field flip from h = -1 to h = +1 at constant temperature.

    Returns
    -------
    torch.Tensor
        A tensor of shape (number_of_steps+2, number_of_control_parameters)
        containing the protocol values at each timestep.
    """
    # Initialize protocol tensor with zeros
    cee = torch.zeros(number_of_steps+2, number_of_control_parameters, device=device)

    cee[1:-1] = torch.from_numpy(np.loadtxt("default_protocol.dat")[::10])
    # Set boundary conditions
    cee[0, :] = cee_boundary[0, :]  # Initial value
    cee[-1, :] = cee_boundary[-1,:]  # Final value

    return cee



def run_protocol(protocol, number_of_trajectories, visualize=False):
    """
    Run a 2D Ising model simulation with the given temperature and field protocol.

    Simulates Monte Carlo dynamics of spins on a square lattice using the Metropolis
    algorithm with checkerboard updating. Calculates entropy production during the
    magnetic field reversal process.

    Parameters
    ----------
    protocol : torch.Tensor
        Control parameter values at each timestep,
        shape (number_of_steps+2, number_of_control_parameters)
        Column 0: Temperature T(t)
        Column 1: Magnetic field h(t)
    number_of_trajectories : int
        Number of independent spin configurations to simulate
    visualize : bool, optional
        Whether to collect spin configurations for visualization (default: False)

    Returns
    -------
    entprod : torch.Tensor
        Entropy production for each trajectory, shape (number_of_trajectories,)
    lattices_to_visualize : torch.Tensor, optional (if visualize=True)
        Spin configurations at visualization timesteps,
        shape (number_of_pictures+2, n_lattices_plot, L, L)
    parameters_to_visualize : torch.Tensor, optional (if visualize=True)
        Protocol values at visualization timesteps,
        shape (number_of_pictures+2, number_of_control_parameters)
    """

    spins = -1*torch.ones((number_of_trajectories, L, L), device = device)

    # Initialize state variables
    # Initial positions drawn from standard normal distribution (equilibrium at λ=0)
    entprod = torch.zeros(number_of_trajectories, device=device)

    # Move protocol to the correct device (GPU if available)
    protocol = protocol.to(device)

    # Preallocate noise tensor for efficiency
    probs = torch.zeros(number_of_trajectories, L, L, device=device)

    if visualize:
        n_lattices_plot = 10
        lattices_to_visualize = torch.zeros(number_of_pictures+2, n_lattices_plot, L,L, device=device)
        parameters_to_visualize = torch.zeros(number_of_pictures+2, number_of_control_parameters, device=device)
        lattices_to_visualize[0] = spins[:n_lattices_plot]
        parameters_to_visualize[0] = protocol[0]

    # Main simulation loop

    i_indices, j_indices = torch.meshgrid(torch.arange(L, device=device), torch.arange(L, device=device), indexing='ij')
    even_sites_mask = (i_indices + j_indices) % 2 == 0
    odd_sites_mask = ~even_sites_mask

    entprod -= (-spins * (torch.roll(spins, 1 , dims = 1) + torch.roll(spins, 1 , dims = 2)) - protocol[0,1] * spins).sum(dim = (1,2))/protocol[0,0]

    for i in range(1,number_of_steps+1):
        for j in range(2):
            if j == 0:
                mask = even_sites_mask
            else:
                mask = odd_sites_mask

            probs.uniform_()
            delta_E = 2*spins * (torch.roll(spins, 1 , dims = 1) + torch.roll(spins, -1 , dims = 1) + torch.roll(spins, 1 , dims = 2) + torch.roll(spins, -1 , dims = 2)) + 2 * protocol[i, 1] * spins
            update_condition = (1/(1+torch.exp(delta_E/protocol[i,0])) > probs) & mask.unsqueeze(0)
            spins = torch.where(update_condition, -spins, spins)
            entprod -= (delta_E * update_condition).sum(dim=(1,2))/protocol[i,0]

        # Store data for visualization at regular intervals
        if visualize and (i-1) % (number_of_steps // number_of_pictures) == 0:
            idx = 1 + (i-1) // ((number_of_steps)  // number_of_pictures)
            lattices_to_visualize[idx, :] = spins[:n_lattices_plot]
            parameters_to_visualize[idx, :] = protocol[i]

    entprod += (-spins * (torch.roll(spins, 1 , dims = 1) +  torch.roll(spins, 1 , dims = 2)) - protocol[-1,1] * spins).sum(dim = (1,2))/protocol[-1,0]


    # Return appropriate results based on visualization flag
    if visualize:
        lattices_to_visualize[-1] = spins[:n_lattices_plot]
        parameters_to_visualize[-1] = protocol[-1]
        return entprod, lattices_to_visualize, parameters_to_visualize
    else:
        return entprod

def visualize_protocol(protocol, number_of_trajectories):
    """
    Visualize the evolution of spin configurations during magnetic field reversal.

    Creates a series of plots showing:
    1. The protocol trajectory in (T,h) parameter space
    2. Spin configurations at different time points during the field flip
    3. A histogram of entropy production values
    4. Time evolution movie of the magnetic phase transition

    Parameters
    ----------
    protocol : torch.Tensor
        Control parameter values at each timestep
    number_of_trajectories : int
        Number of spin configurations to simulate

    Returns
    -------
    None
        Plots are saved to the 'plots' directory and a movie is generated
    """
    # Create plots directory if it doesn't exist
    os.makedirs('plots', exist_ok=True)
    fig, ax = plt.subplots(figsize=(6, 6))
    plt.plot(protocol[:,0], protocol[:,1], color = "g")
    # plt.plot(np.linspace(0,trajectory_time, number_of_steps+2), protocol[:,1], color = "b")
    plt.xlabel(r"$T$")
    plt.ylabel(r"$h$")
    fig.tight_layout()
    fig.savefig(f'protocol.png', pad_inches=0.05, dpi=300)


    # Run the simulation with visualization enabled
    entprod, lattices_to_visualize, parameters_to_visualize = run_protocol(
        protocol, number_of_trajectories, visualize=True
    )

    # Create 2x3 panel figure

    # Generate a series of plots showing the evolution of the system
    for i in range(number_of_pictures+2):
        panel_fig, panel_axes = plt.subplots(2, 3, figsize=(12, 9))
        panel_axes = panel_axes.flatten()

        # Add to 2x3 panel figure if this is a selected time point
        panel_axes[3].plot(protocol[:,0], protocol[:,1], color = "k")
        panel_axes[3].scatter(protocol[i,0], protocol[i,1], color = "k")

        panel_axes[3].set_xlabel(r"$T$")
        panel_axes[3].set_ylabel(r"$h$", rotation='horizontal', ha='right')
        for j in range(6):
            if j == 3:
                pass
            else:
                ax_panel = panel_axes[j]

                # Plot the same data on the panel
                ax_panel.matshow(lattices_to_visualize[i,j], cmap='gray', vmin=-1, vmax=1)
                ax_panel.set_xticks([])
                ax_panel.set_yticks([])

        panel_fig.tight_layout()
        panel_fig.savefig(f'plots/output_{i}.png', pad_inches=0.05, dpi=300)
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
    hist, bin_edges = np.histogram(entprod, bins=n_bins, density=True)

    # Compute bin centers for plotting
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    print(f"{torch.mean(entprod).item():.3f}")
    # Plot work distribution
    plt.plot(bin_centers, hist, color="g")
    plt.axvline(x = torch.mean(entprod), ls = "--", color = "g")
    plt.text(s = r"$\langle \sigma \rangle = $" + f"{torch.mean(entprod).item():.3f}", x = 1.1*torch.mean(entprod), y = np.max(hist), ha = "left")
    plt.xlabel(r"$\sigma$")
    plt.ylabel(r"$P(\sigma)$")
    plt.ylim([0, 1.1*np.max(hist)])
    plt.xlim([np.min(bin_centers), np.max(bin_centers)])
    # Save work histogram
    fig.tight_layout()
    fig.savefig("entprod_histogram.png", pad_inches=0.05, dpi=300)
    plt.close(fig)

def calculate_order_parameter(protocol, number_of_trajectories=int(1e5)):
    """
    Calculate the order parameter for the magnetic field reversal protocol.

    The order parameter is the mean entropy production during the field flip,
    which quantifies the thermodynamic cost of the nonequilibrium process.
    Optimal protocols minimize this quantity.

    Parameters
    ----------
    protocol : torch.Tensor
        Control parameter values at each timestep
    number_of_trajectories : int, optional
        Number of trajectories to simulate (default: 1e5)

    Returns
    -------
    float
        The mean entropy production (order parameter)
    """
    entprod = run_protocol(protocol, number_of_trajectories)
    return torch.mean(entprod).item()

def final_answer(protocol):
    """
    Calculate the entropy production statistics with high precision.

    Runs a large number of simulations (10^6 trajectories) to accurately
    estimate the mean entropy production and its standard error using
    batch analysis for improved statistical reliability.

    Parameters
    ----------
    protocol : torch.Tensor
        Control parameter values at each timestep

    Returns
    -------
    None
        Prints the mean entropy production value and standard error to the console
    """
    # Run a large number of trajectories for statistical significance
    entprod = run_protocol(protocol, int(1e6))

    # Reshape entropy production array into 100 batches of 10^4 trajectories each
    # This allows for better statistical analysis through batch means
    entprod = torch.reshape(entprod, (int(1e2), int(1e4)))

    # Calculate mean entropy production for each batch
    sample_means = torch.mean(entprod, axis=1)

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
    protocol = load_default_protocol()
    visualize_protocol(protocol, int(1e4))
    final_answer(protocol)

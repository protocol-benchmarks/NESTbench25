import numpy as np
import matplotlib.pyplot as plt
import glob
import imageio
import os

# Function to read data from a file
def read_data(file_name):
    try:
        data = np.loadtxt(file_name)
        if data.ndim == 1 and data.shape[0] != 2:
            raise ValueError(f"File '{file_name}' does not contain two columns.")
        return (data[0:1], data[1:2]) if data.ndim == 1 else (data[:, 0], data[:, 1])
    except Exception as e:
        raise RuntimeError(f"Failed to read '{file_name}': {e}")

# Find all relevant files and determine N
pos_files = glob.glob('report_position_time_*.dat')
boltz_files = glob.glob('report_boltz_time_*.dat')
potential_files = glob.glob('report_potential_time_*.dat')

# Determine N based on the number of pos files (assuming all three sets of files are consistent)
N = len(pos_files)

# Set the font to Computer Modern (LaTeX's default font)
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = 'Computer Modern'
plt.rcParams['text.usetex'] = True
plt.rcParams['font.size'] = 16
plt.rcParams['xtick.labelsize'] = 16
plt.rcParams['ytick.labelsize'] = 16
plt.rcParams['legend.fontsize'] = 16
plt.rcParams['axes.titlesize'] = 16
plt.rcParams['axes.labelsize'] = 16

for n in range(N):
    # Construct file names based on the current value of n
    pos_file = f'report_position_time_{n}.dat'
    boltz_file = f'report_boltz_time_{n}.dat'
    potential_file = f'report_potential_time_{n}.dat'

    # Check if all files exist for this value of n
    if not (pos_file in pos_files and boltz_file in boltz_files and potential_file in potential_files):
        print(f"Missing files for n={n}, skipping...")
        continue

    # Read data from files
    try:
        x_pos, y_pos = read_data(pos_file)
        x_boltz, y_boltz = read_data(boltz_file)
        x_potential, y_potential = read_data(potential_file)
    except RuntimeError as e:
        print(e)
        continue

    # Scale pos and boltz vertical values by 20
    y_pos_scaled = y_pos * 20
    y_boltz_scaled = y_boltz * 20

    # Create a figure and axis with higher resolution
    fig, ax = plt.subplots(figsize=(10.0267, 6.0267), dpi=300)  # Increase resolution with dpi and larger size

    # Plot the data
    ax.plot(x_pos, y_pos_scaled, 'g-', label='pos')
    ax.plot(x_boltz, y_boltz_scaled, 'k--', label='boltz')
    ax.plot(x_potential, y_potential, 'k', label='potential')

    # Set axis limits
    ax.set_ylim(0, 10)
    ax.set_xlim(-10, 10)

    # Customize x-axis ticks
    ax.set_xticks([-5, 0, 5])
    ax.set_xlabel(r"$x$", fontsize=16)

    # Add legend
    ax.set_yticks([])  # Remove vertical axis values

    # Save the figure as a PNG
    output_file = f'output_{n}.png'
    fig.savefig(output_file, format='png', dpi=300)  # Higher resolution
    plt.close(fig)  # Close the figure to free up memory

    print(f"PNG saved as '{output_file}'")

# Create a list of the output file names
output_files = [f'output_{n}.png' for n in range(N) if f'report_position_time_{n}.dat' in pos_files]

# Generate a figure with 6 equally-spaced frames
num_frames = 6
selected_indices = np.linspace(0, N - 1, num_frames, dtype=int)

fig, axs = plt.subplots(2, 3, figsize=(12.0267, 6.0267), dpi=300)
axs = axs.flatten()

for i, n in enumerate(selected_indices):
    pos_file = f'report_position_time_{n}.dat'
    boltz_file = f'report_boltz_time_{n}.dat'
    potential_file = f'report_potential_time_{n}.dat'

    try:
        x_pos, y_pos = read_data(pos_file)
        x_boltz, y_boltz = read_data(boltz_file)
        x_potential, y_potential = read_data(potential_file)
    except RuntimeError as e:
        print(e)
        continue

    y_pos_scaled = y_pos * 20
    y_boltz_scaled = y_boltz * 20

    ax = axs[i]
    ax.plot(x_pos, y_pos_scaled, 'g-', label='pos')
    ax.plot(x_boltz, y_boltz_scaled, 'k--', label='boltz')
    ax.plot(x_potential, y_potential, 'k', label='potential')
    ax.set_xlim(-10, 10)
    ax.set_ylim(0, 10)
    ax.set_xticks([-5, 0, 5])
    ax.set_xlabel(r"$x$", fontsize=16)
    ax.set_yticks([])
    ax.set_title(r"$t = {:.2f}$".format(n / (N - 1)), fontsize=16)

for ax in axs:
    ax.label_outer()

fig.tight_layout()
fig.savefig("output_figure.pdf", format='pdf', dpi=300)
print("figure saved as 'output_figure.pdf'")

# Stitch the images into a movie
with imageio.get_writer('output_movie.mp4', fps=7) as writer:
    for filename in output_files:
        image = imageio.imread(filename)
        writer.append_data(image)

print("movie saved as 'output_movie.mp4'")

# Remove intermediate PNG files
for filename in output_files:
    os.remove(filename)

# Remove input data files
for n in range(N):
    for prefix in ['report_position_time_', 'report_boltz_time_', 'report_potential_time_']:
        try:
            os.remove(f'{prefix}{n}.dat')
        except FileNotFoundError:
            pass

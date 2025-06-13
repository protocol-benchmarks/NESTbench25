import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import imageio.v2 as imageio
from matplotlib.colors import ListedColormap

# Find all lattice files
lattice_files = sorted(glob.glob('report_lattice_time_*.dat'), key=lambda x: int(x.split('_')[-1].split('.')[0]))
N = len(lattice_files)

# Custom colormap: white for -1, light blue for +1
cmap = ListedColormap(["white", "#add8e6"])  # light blue hex

# Load lattice data
def load_lattice(filename):
    return np.loadtxt(filename)

# Select 9 evenly spaced indices for grid
num_frames = 9
selected_indices = np.linspace(0, N - 1, num_frames, dtype=int)

# Plot grid of selected snapshots
fig, axs = plt.subplots(3, 3, figsize=(9, 9), dpi=300)
axs = axs.flatten()

for i, idx in enumerate(selected_indices):
    lattice = load_lattice(lattice_files[idx])
    ax = axs[i]
    ax.imshow((lattice + 1) // 2, cmap=cmap, vmin=0, vmax=1, origin='lower')
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(f"$t = {idx / (N - 1):.2f}$", fontsize=12)
    for spine in ax.spines.values():
        spine.set_edgecolor('black')
        spine.set_linewidth(0.5)

fig.tight_layout()
fig.savefig("output_picture.pdf", format='pdf', dpi=300)
plt.close(fig)

# Create animated movie from all frames
frame_files = []
for i, fname in enumerate(lattice_files):
    lattice = load_lattice(fname)
    fig, ax = plt.subplots(figsize=(4.267, 4.267), dpi=150)
    ax.imshow((lattice + 1) // 2, cmap=cmap, vmin=0, vmax=1, origin='lower')
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_edgecolor('black')
        spine.set_linewidth(0.5)
    frame_path = f"_frame_{i:04d}.png"
    fig.savefig(frame_path)
    plt.close(fig)
    frame_files.append(frame_path)

with imageio.get_writer("output_movie.mp4", fps=6) as writer:
    for fname in frame_files:
        writer.append_data(imageio.imread(fname))

# Clean up temporary frame files
for fname in frame_files:
    os.remove(fname)

# Clean up report_lattice input files
for fname in lattice_files:
    os.remove(fname)

print("Lattice figure saved as 'output_picture.pdf'")
print("Animation saved as 'output_movie.mp4'")

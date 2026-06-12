import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import re
import argparse

def parse_profiling_data(filename):
    """Parse the profiling results file and extract data."""
    sign_assignment_data = []
    degree_data = []

    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            # Parse sign assignment times (no degree specified)
            sign_match = re.search(r'Time taken for sign assignment at length (\d+): ([\d.]+) seconds, Max memory: ([\d.]+) MiB', line)
            if sign_match:
                length = int(sign_match.group(1))
                time = float(sign_match.group(2))
                memory = float(sign_match.group(3))
                sign_assignment_data.append((length, time, memory))

            # Parse degree-specific data
            degree_match = re.search(r'Braid length: (\d+), Degree: (\d+), Time taken: ([\d.]+) seconds, Max memory: ([\d.]+) MiB', line)
            if degree_match:
                length = int(degree_match.group(1))
                degree = int(degree_match.group(2))
                time = float(degree_match.group(3))
                memory = float(degree_match.group(4))
                degree_data.append((length, degree, time, memory))

    return sign_assignment_data, degree_data

def create_plots(filename):
    """Create the two requested plots."""
    # Parse the data
    sign_data, degree_data = parse_profiling_data(filename)

    # Create figure with subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))

    # Plot 1: 1D curves - Sign assignment times and memory vs braid length
    if sign_data:
        lengths = [item[0] for item in sign_data]
        times = [item[1] for item in sign_data]
        memories = [item[2] for item in sign_data]

        ax1_twin = ax1.twinx()

        line1 = ax1.plot(lengths, times, color='#2E8B57', marker='o', label='Sign Assignment Time', linewidth=2, markersize=6)
        line2 = ax1_twin.plot(lengths, memories, color='#FF6347', marker='s', label='Sign Assignment Memory', linewidth=2, markersize=6)

        ax1.set_xlabel('Braid Length')
        ax1.set_ylabel('Time (seconds)', color='#2E8B57')
        ax1_twin.set_ylabel('Memory (MiB)', color='#FF6347')
        ax1.set_title('Sign Assignment Performance vs Braid Length')
        ax1.tick_params(axis='y', labelcolor='#2E8B57')
        ax1_twin.tick_params(axis='y', labelcolor='#FF6347')
        ax1.grid(True, alpha=0.3)

        # Combine legends
        lines = line1 + line2
        labels = [l.get_label() for l in lines]
        ax1.legend(lines, labels, loc='upper left')

    # Plot 2: 3D surface plots - (degree, braid length) vs time and memory
    if degree_data:
        # Extract unique lengths and degrees
        lengths_deg = sorted(list(set([item[0] for item in degree_data])))
        degrees = sorted(list(set([item[1] for item in degree_data])))

        # Create meshgrid
        L, D = np.meshgrid(lengths_deg, degrees)

        # Initialize time and memory arrays
        Time = np.full(L.shape, np.nan)
        Memory = np.full(L.shape, np.nan)

        # Fill arrays with data
        for length, degree, time, memory in degree_data:
            i = degrees.index(degree)
            j = lengths_deg.index(length)
            Time[i, j] = time
            Memory[i, j] = memory

        # Create 3D plot
        ax2.remove()
        ax2 = fig.add_subplot(122, projection='3d')

        # Plot time surface
        surf1 = ax2.plot_surface(L, D, Time, alpha=0.8, cmap='coolwarm', label='Time')

        # Plot memory surface (offset slightly for visibility)
        Memory_scaled = Memory / np.nanmax(Memory) * np.nanmax(Time)  # Scale memory to time range
        surf2 = ax2.plot_surface(L, D, Memory_scaled, alpha=0.6, cmap='summer', label='Memory (scaled)')

        ax2.set_xlabel('Braid Length')
        ax2.set_ylabel('Degree')
        ax2.set_zlabel('Time (seconds)')
        ax2.set_title('Numerical Algebra Performance vs (Degree, Braid Length)')

        # Add color bars
        cbar1 = plt.colorbar(surf1, ax=ax2, shrink=0.5, aspect=20, pad=0.1)
        cbar1.set_label('Time (seconds)')

        cbar2 = plt.colorbar(surf2, ax=ax2, shrink=0.5, aspect=20, pad=0.2)
        cbar2.set_label('Memory (scaled to time range)')

        # Add text annotation for memory scaling
        ax2.text2D(0.02, 0.98, 'Green-yellow surface: Memory (scaled to time range)',
                   transform=ax2.transAxes, fontsize=8, verticalalignment='top')

    plt.tight_layout()
    plt.savefig('profiling_plots.png', dpi=300, bbox_inches='tight')
    plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plot profiling results from a data file')
    parser.add_argument('filename', help='Path to the profiling results file')
    args = parser.parse_args()
    
    create_plots(args.filename)

import numpy as np
import matplotlib.pyplot as plt

# List of file paths
file_paths = [
    "/home/eros/workspace-lsm/wildturkey/dbbench-result/lifetime/64output_keys_level_0.txt",
    "/home/eros/workspace-lsm/wildturkey/dbbench-result/lifetime/64output_keys_level_1.txt",
    "/home/eros/workspace-lsm/wildturkey/dbbench-result/lifetime/64output_keys_level_2.txt"
]

# Function to read timestamps from a file
def read_timestamps(file_path):
    with open(file_path, 'r') as file:
        timestamps = [int(line.strip()) for line in file if line.strip().isdigit()]
    return timestamps

# Initialize plot
plt.figure(figsize=(8, 5))

# Colors and labels for each file
colors = ['b', 'g', 'r']
labels = ['Level 0', 'Level 1', 'Level 2']

# Reading timestamps and plotting CDF for each file
for file_path, color, label in zip(file_paths, colors, labels):
    timestamps = read_timestamps(file_path)
    timestamps_sorted = np.sort(timestamps)
    cdf = np.arange(1, len(timestamps_sorted) + 1) / len(timestamps_sorted)
    plt.plot(timestamps_sorted, cdf, marker='o', linestyle='-', color=color, label=label)

# Set x-axis to logarithmic scale
plt.xscale('log')
plt.title('64 CDF of Timestamps (Log Scale) for Each Level')
plt.xlabel('Timestamps (Log Scale)')
plt.ylabel('CDF')
plt.grid(True)
plt.legend()

# Save the plot
output_image_path = "/home/eros/workspace-lsm/wildturkey/dbbench-result/lifetime/64cdf_plot_all_levels_log.png"
plt.savefig(output_image_path, dpi=300)
plt.close()

print(f"CDF plot for all levels saved as {output_image_path}.")

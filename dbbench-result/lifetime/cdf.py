# Python code to read timestamps from a text file and plot their CDF

import numpy as np
import matplotlib.pyplot as plt

# Path to the text file
file_path = "/home/eros/workspace-lsm/wildturkey/dbbench-result/lifetime/64output_keys_level_1.txt"  # Update this path to your actual file

# Function to read timestamps from the file
def read_timestamps(file_path):
    with open(file_path, 'r') as file:
        timestamps = [int(line.strip()) for line in file if line.strip().isdigit()]
    return timestamps

# Reading the timestamps
timestamps = read_timestamps(file_path)

# Sorting timestamps and calculating CDF
timestamps_sorted = np.sort(timestamps)
cdf = np.arange(1, len(timestamps_sorted) + 1) / len(timestamps_sorted)

# Plotting the CDF
plt.figure(figsize=(8, 5))
plt.plot(timestamps_sorted, cdf, marker='o', linestyle='-', color='b')
plt.title('CDF of Timestamps')
plt.xlabel('Timestamps')
plt.ylabel('CDF')
plt.xscale('log')
plt.grid(True)
output_image_path = "/home/eros/workspace-lsm/wildturkey/dbbench-result/lifetime/64cdf_plot1.png"  # You can customize the filename
plt.savefig(output_image_path, dpi=300)
plt.close()


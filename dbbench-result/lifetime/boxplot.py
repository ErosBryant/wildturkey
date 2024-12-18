import numpy as np
import matplotlib.pyplot as plt

# List of file paths
file_paths = [
    "/home/eros/workspace-lsm/wildturkey/dbbench-result/lifetime/output_keys_level_0.txt",
    "/home/eros/workspace-lsm/wildturkey/dbbench-result/lifetime/output_keys_level_1.txt",
    "/home/eros/workspace-lsm/wildturkey/dbbench-result/lifetime/output_keys_level_2.txt",
    "/home/eros/workspace-lsm/wildturkey/dbbench-result/lifetime/output_keys_level_3.txt"
]

# Function to read timestamps from a file
def read_timestamps(file_path):
    with open(file_path, 'r') as file:
        timestamps = [int(line.strip()) for line in file if line.strip().isdigit()]
    return timestamps

# Reading timestamps from all files
data = [read_timestamps(file_path) for file_path in file_paths]

# Creating the box plot
plt.figure(figsize=(8, 5))
plt.boxplot(data, patch_artist=True, showmeans=True, 
            boxprops=dict(facecolor='lightgreen', color='black'), 
            meanprops=dict(marker='o', markerfacecolor='red', markersize=8), 
            medianprops=dict(color='red'))

# Customizing the plot
plt.xticks([1, 2, 3, 4], ['Level 0', 'Level 1', 'Level 2', 'Level 3'])
plt.xlabel('Level')
plt.ylabel('life time')
plt.title('2')
plt.grid(axis='y')

# Save the plot
output_image_path = "/home/eros/workspace-lsm/wildturkey/dbbench-result/lifetime/2box_plot.png"
plt.savefig(output_image_path, dpi=300)
plt.close()

print(f"Box plot saved as {output_image_path}.")

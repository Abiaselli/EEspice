import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Load the data from CSV files
threads_time_df = pd.read_csv('Threads time.csv')
simulation_end_time_df = pd.read_csv('simulation end time.csv')

# Strip whitespace from column names
threads_time_df.columns = threads_time_df.columns.str.strip()

# Calculate the time gap between the longest and shortest task per timestep
tgaps = []

for i in range(len(threads_time_df)):
    durations = [
        threads_time_df.at[i, 'End Mid'] - threads_time_df.at[i, 'Begin Mid'],
        threads_time_df.at[i, 'End Up'] - threads_time_df.at[i, 'Begin Up'],
        threads_time_df.at[i, 'End Down'] - threads_time_df.at[i, 'Begin Down']
    ]
    tgap = max(durations) - min(durations)
    tgaps.append(tgap)

# Convert tgaps to a numpy array for easier manipulation
tgaps = np.array(tgaps)

# Calculate the mean and standard deviation
mean_tgap = np.mean(tgaps)
std_tgap = np.std(tgaps)

# Print the mean and standard deviation
print(f"Mean of tgaps: {mean_tgap}")
print(f"Standard Deviation of tgaps: {std_tgap}")

# Create a histogram of the time gaps
plt.figure(figsize=(10, 6))
plt.hist(tgaps, bins=20, log=True)


# Customize the scales to be in log base 10
plt.xscale('log', base=10)
plt.yscale('log', base=10)


# Set labels and title
plt.xlabel('Time Gap (ms) - Log10 Scale')
plt.ylabel('Frequency - Log10 Scale')
plt.title('Histogram of Time Gaps (tgap) Between Longest and Shortest Task per Iteration')

# plt.xlabel('Time Gap (ms)')
# plt.ylabel('Frequency')
# plt.title('Histogram of Time Gaps (tgap) Between Longest and Shortest Task per Iteration wwithout log scale')

# Display the plot
plt.show()
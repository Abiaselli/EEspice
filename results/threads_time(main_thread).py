import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Load the data from CSV files
threads_time_df = pd.read_csv('Threads time.csv')
simulation_end_time_df = pd.read_csv('simulation end time.csv')

# Strip whitespace from column names
threads_time_df.columns = threads_time_df.columns.str.strip()

# Extract simulation end time
simulation_end_time = simulation_end_time_df.iloc[0, 0]

# Define the threads and colors
threads = ['Main', 'Mid', 'Up', 'Down']
colors = {'Mid': 'blue', 'Up': 'green', 'Down': 'red', 'Main': 'purple'}

# Calculate the main thread's active periods
main_active_periods = []
current_time = 0

for i in range(len(threads_time_df)):
    # Find the minimum start time and maximum end time for the current set of threads
    min_start_time = min(threads_time_df.at[i, 'Begin Mid'], threads_time_df.at[i, 'Begin Up'], threads_time_df.at[i, 'Begin Down'])
    max_end_time = max(threads_time_df.at[i, 'End Mid'], threads_time_df.at[i, 'End Up'], threads_time_df.at[i, 'End Down'])
    
    # Append the active period before the current threads start
    main_active_periods.append((current_time, min_start_time))
    
    # Update the current time to the maximum end time
    current_time = max_end_time

# Append the final active period after the last threads finish
main_active_periods.append((current_time, simulation_end_time))

# Create the plot
plt.figure(figsize=(10, 6))

# Plot the main thread's active periods first
for start, end in main_active_periods:
    plt.barh('Main', end - start, left=start, color=colors['Main'])

# Plot the other threads
for thread in ['Mid', 'Up', 'Down']:
    begin_col = f'Begin {thread}'
    end_col = f'End {thread}'
    for i in range(len(threads_time_df)):
        plt.barh(thread, threads_time_df[end_col][i] - threads_time_df[begin_col][i], left=threads_time_df[begin_col][i], color=colors[thread])

# Set the x-axis limit from 0 to the simulation end time
plt.xlim(0, simulation_end_time)

# Set labels and title
plt.xlabel('Time (ms)')
plt.ylabel('Threads')
plt.yticks(ticks=[0, 1, 2, 3], labels=['Main', 'Mid', 'Up', 'Down'])  # Explicitly set the order of threads on y-axis
plt.gca().invert_yaxis()  # Invert the y-axis to have Main at the top
plt.title('Threads Execution Time')

# Display the plot
plt.show()

# First 3 rows:

# Filter the first three rows
filtered_threads_time_df = threads_time_df.head(3)

# Determine the x-axis limits based on the filtered data
min_time = 0
max_time = filtered_threads_time_df[['End Mid', 'End Up', 'End Down']].max().max() + 1

# Calculate the main thread's active periods for the filtered data
main_active_periods_filtered = []
current_time_filtered = 0

for i in range(len(filtered_threads_time_df)):
    # Find the minimum start time and maximum end time for the current set of threads
    min_start_time_filtered = min(filtered_threads_time_df.at[i, 'Begin Mid'], filtered_threads_time_df.at[i, 'Begin Up'], filtered_threads_time_df.at[i, 'Begin Down'])
    max_end_time_filtered = max(filtered_threads_time_df.at[i, 'End Mid'], filtered_threads_time_df.at[i, 'End Up'], filtered_threads_time_df.at[i, 'End Down'])
    
    # Append the active period before the current threads start
    main_active_periods_filtered.append((current_time_filtered, min_start_time_filtered))
    
    # Update the current time to the maximum end time
    current_time_filtered = max_end_time_filtered

# Append the final active period after the last threads finish
main_active_periods_filtered.append((current_time_filtered, max_time))

# Create the plot for the first 3 rows
plt.figure(figsize=(10, 6))

# Plot the main thread's active periods for the filtered data first
for start, end in main_active_periods_filtered:
    plt.barh('Main', end - start, left=start, color=colors['Main'])

# Plot the other threads for the filtered data
for thread in ['Mid', 'Up', 'Down']:
    begin_col = f'Begin {thread}'
    end_col = f'End {thread}'
    for i in range(len(filtered_threads_time_df)):
        plt.barh(thread, filtered_threads_time_df[end_col][i] - filtered_threads_time_df[begin_col][i], left=filtered_threads_time_df[begin_col][i], color=colors[thread])

# Set the x-axis limit based on the relevant data range
plt.xlim(min_time, max_time)

# Set labels and title
plt.xlabel('Time (ms)')
plt.ylabel('Threads')
plt.yticks(ticks=[0, 1, 2, 3], labels=['Main', 'Mid', 'Up', 'Down'])  # Explicitly set the order of threads on y-axis
plt.gca().invert_yaxis()  # Invert the y-axis to have Main at the top
plt.title('Threads Execution Time (First 3 Rows)')

# Display the plot
plt.show()
# -------------------------------------------------------------------------------

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

# Sum of tgaps
sum_tgaps = np.sum(tgaps)

# Calculate the total active time of the main thread
sum_main_active = sum([end - start for start, end in main_active_periods])

# Calculate the total time
total_time = simulation_end_time

# Calculate util1
util1 = 1 - sum_tgaps / (total_time - sum_main_active)

# Calculate util2
util2 = 1 - sum_tgaps / total_time

# Print util2
print(f"Utilization (util2): {util2}")

# Print util1
print(f"Utilization (util1): {util1}")
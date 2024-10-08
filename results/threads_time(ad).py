import pandas as pd
import matplotlib.pyplot as plt

# Load the data from CSV files
mid_df = pd.read_csv('Threads time mid.csv')
up_df = pd.read_csv('Threads time up.csv')
down_df = pd.read_csv('Threads time down.csv')
simulation_end_time_df = pd.read_csv('simulation end time.csv')

# Strip whitespace from column names
mid_df.columns = mid_df.columns.str.strip()
up_df.columns = up_df.columns.str.strip()
down_df.columns = down_df.columns.str.strip()

# Extract simulation end time
simulation_end_time = simulation_end_time_df.iloc[0, 0]

# Define the threads and colors
threads = ['Main', 'Mid', 'Up', 'Down']
colors = {'Mid': 'blue', 'Up': 'green', 'Down': 'red', 'Main': 'purple'}

# Create the plot
plt.figure(figsize=(10, 6))

# Function to plot thread intervals
def plot_intervals(df, thread, color):
    for i in range(len(df)):
        start_time = df.iloc[i, 0]
        end_time = df.iloc[i, 1]
        plt.barh(thread, end_time - start_time, left=start_time, color=color)


# Combine all intervals from Mid, Up, and Down and sort them
all_intervals = pd.concat([
    mid_df.rename(columns={'Begin Mid': 'Begin', 'End Mid': 'End'}),
    up_df.rename(columns={'Begin Up': 'Begin', 'End Up': 'End'}),
    down_df.rename(columns={'Begin Down': 'Begin', 'End Down': 'End'})
])
all_intervals_sorted = all_intervals.sort_values(by='Begin').reset_index(drop=True)

# Calculate the main thread's active periods
main_active_periods = []
current_time = 0

for i in range(len(all_intervals_sorted)):
    start_time = all_intervals_sorted.iloc[i, 0]
    end_time = all_intervals_sorted.iloc[i, 1]

    # Append the active period before the current threads start
    if current_time < start_time:
        main_active_periods.append((current_time, start_time))
    
    # Update the current time to the maximum end time
    if current_time < end_time:
        current_time = end_time

# Append the final active period after the last threads finish
if current_time < simulation_end_time:
    main_active_periods.append((current_time, simulation_end_time))

# Calculate total main thread time
total_main_thread_time = sum(end - start for start, end in main_active_periods)

# Plot the main thread's active periods
for start, end in main_active_periods:
    plt.barh('Main', end - start, left=start, color=colors['Main'])

# Plot intervals for each thread
plot_intervals(mid_df, 'Mid', colors['Mid'])
plot_intervals(up_df, 'Up', colors['Up'])
plot_intervals(down_df, 'Down', colors['Down'])

# Set the x-axis limit from 0 to the simulation end time
plt.xlim(0, simulation_end_time)

# Set labels and title
plt.xlabel('Time (ms)')
plt.ylabel('Threads')
plt.yticks(ticks=[0, 1, 2, 3], labels=['Main', 'Mid', 'Up', 'Down'])  # Explicitly set the order of threads on y-axis
plt.gca().invert_yaxis()  # Invert the y-axis to have 'Main' at the top
plt.title('Threads Execution Time')

# Display the plot
plt.show()

# Calculate total utilization
total_thread_time = (total_main_thread_time +
                    mid_df.apply(lambda row: row['End Mid'] - row['Begin Mid'], axis=1).sum() +
                    up_df.apply(lambda row: row['End Up'] - row['Begin Up'], axis=1).sum() +
                    down_df.apply(lambda row: row['End Down'] - row['Begin Down'], axis=1).sum())

total_utilization = total_thread_time / simulation_end_time
print("Total Utilization:", total_utilization)
import pandas as pd
import matplotlib.pyplot as plt

# Load the data from CSV files
threads_time_df = pd.read_csv('Threads time.csv')
simulation_end_time_df = pd.read_csv('simulation end time.csv')

# Strip whitespace from column names
threads_time_df.columns = threads_time_df.columns.str.strip()

# Print column names to debug
print("Columns in Threads time.csv after stripping whitespace:", threads_time_df.columns)

# Extract simulation end time
simulation_end_time = simulation_end_time_df.iloc[0, 0]

# Define the threads and colors
threads = ['Mid', 'Up', 'Down']
colors = {'Mid': 'blue', 'Up': 'green', 'Down': 'red'}

# Create the plot
plt.figure(figsize=(10, 6))

for thread in threads:
    begin_col = f'Begin {thread}'
    end_col = f'End {thread}'
    for i in range(len(threads_time_df)):
        plt.barh(thread, threads_time_df[end_col][i] - threads_time_df[begin_col][i], left=threads_time_df[begin_col][i], color=colors[thread])

# Set the x-axis limit from 0 to the simulation end time
plt.xlim(0, simulation_end_time)

# Set labels and title
plt.xlabel('Time (ms)')
plt.ylabel('Threads')
plt.title('Threads Execution Time')

# Display the plot
plt.show()

# First 3 rows:

# Load the data from CSV files
threads_time_df = pd.read_csv('Threads time.csv')
simulation_end_time_df = pd.read_csv('simulation end time.csv')

# Strip whitespace from column names
threads_time_df.columns = threads_time_df.columns.str.strip()

# Extract simulation end time
simulation_end_time = simulation_end_time_df.iloc[0, 0]

# Filter the first three rows
filtered_threads_time_df = threads_time_df.head(3)

# Define the threads and colors
threads = ['Mid', 'Up', 'Down']
colors = {'Mid': 'blue', 'Up': 'green', 'Down': 'red'}

# Determine the x-axis limits based on the filtered data
# min_time = filtered_threads_time_df[['Begin Mid', 'Begin Up', 'Begin Down']].min().min()
min_time = 0
max_time = filtered_threads_time_df[['End Mid', 'End Up', 'End Down']].max().max()+1

# Create the plot
plt.figure(figsize=(10, 6))

for thread in threads:
    begin_col = f'Begin {thread}'
    end_col = f'End {thread}'
    for i in range(len(filtered_threads_time_df)):
        plt.barh(thread, filtered_threads_time_df[end_col][i] - filtered_threads_time_df[begin_col][i], left=filtered_threads_time_df[begin_col][i], color=colors[thread])

# Set the x-axis limit based on the relevant data range
plt.xlim(min_time, max_time)

# Set labels and title
plt.xlabel('Time (ms)')
plt.ylabel('Threads')
plt.title('Threads Execution Time (First 3 Rows)')

# Display the plot
plt.show()
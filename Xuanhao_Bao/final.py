import matplotlib.pyplot as plt
import pandas as pd

# Load the data from the CSV file
df = pd.read_csv('final_solution.csv', header=None)
# Correctly name the columns based on your data structure
df.columns = ['Time', 'Time Step', 'Voltage 1', 'Voltage 2', 'Voltage 3', 'Current', 'Current 2']

# Ensure the CSV structure matches the expected format
# It's crucial that your CSV data is correctly formatted to match the expected six columns

# Plotting
plt.figure(figsize=(12, 8))

# Voltage 1, Voltage 2, and Voltage 3
plt.subplot(2, 1, 1)  # 2 rows, 1 column, 1st subplot
plt.plot(df['Time'], df['Voltage 1'], label='Voltage 1', marker='o')  # Time is the first column, Voltage 1 is the third
plt.plot(df['Time'], df['Voltage 2'], label='Voltage 2', marker='x')  # Voltage 2 is the fourth column
plt.plot(df['Time'], df['Voltage 3'], label='Voltage 3', color='red', marker='^')  # Voltage 3 is the fifth column
plt.xlabel('Time')
plt.ylabel('Voltage')
plt.title('Voltage vs. Time')
plt.legend()

# Current
plt.subplot(2, 1, 2)  # 2 rows, 1 column, 2nd subplot
plt.plot(df['Time'], df['Current 2'], label='Current', color='green', marker='^')  # Current is the sixth column
plt.xlabel('Time')
plt.ylabel('Current')
plt.title('Current vs. Time')
plt.legend()

plt.tight_layout()
plt.show()

# Plotting Current 2 vs. Voltage 2
plt.figure(figsize=(12, 8))
plt.plot(df['Voltage 2'], -df['Current 2'], label='Current vs. Voltage 2', marker='^', color='blue')
plt.xlabel('Voltage 2')
plt.ylabel('Current 2')
plt.title('Current 2 vs. Voltage 2')
plt.legend()
plt.grid(True)  # Optional: Adds a grid for better readability
plt.show()


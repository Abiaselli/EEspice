import matplotlib.pyplot as plt
import pandas as pd

# Load the data from the CSV file
df = pd.read_csv('final_solution.csv', header=None)

# If your CSV file includes headers as the first row, remove ', header=None'
# and ensure your columns are named appropriately in the file, or rename them after loading:
df.columns = ['Time', 'Time Step', 'Voltage 1', 'Voltage 2', 'Voltage 3', 'Current']

# Plotting
plt.figure(figsize=(12, 8))

# Voltage 1 and Voltage 2
plt.subplot(2, 1, 1) # 2 rows, 1 column, 1st subplot
plt.plot(df.iloc[:, 0], df.iloc[:, 2], label='Voltage 1', marker='o') # Assuming Time is the first column, Voltage 1 is the third
plt.plot(df.iloc[:, 0], df.iloc[:, 3], label='Voltage 2', marker='x') # Assuming Voltage 2 is the fourth column
plt.plot(df.iloc[:, 0], df.iloc[:, 4], label='Voltage 3', color='red', marker='^') # Assuming Current is the fifth column
plt.xlabel('Time')
plt.ylabel('Voltage')
plt.title('Voltage vs. Time')
plt.legend()

# Current
plt.subplot(2, 1, 2) # 2 rows, 1 column, 2nd subplot
plt.plot(df.iloc['Time'], df.iloc['Current'], label='Current', color='green', marker='^') # Assuming Current is the fifth column
plt.xlabel('Time')
plt.ylabel('Current')
plt.title('Current vs. Time')
plt.legend()

plt.tight_layout()
plt.show()

import matplotlib.pyplot as plt
import pandas as pd

# Load the data
df = pd.read_csv('final_solution.csv')
df.columns = df.columns.str.strip()  # Ensure column names are trimmed of any whitespace

# Determine which columns are for voltage and current
voltage_columns = [col for col in df.columns if 'Voltage' in col]
current_columns = [col for col in df.columns if 'Current' in col]

# Set up figures
fig1, ax1 = plt.subplots(figsize=(9, 6))
fig2, axes = plt.subplots(3, 1, figsize=(9, 18), sharex=True)  # 3 rows, 1 column for Figure 2
fig3, ax3 = plt.subplots(figsize=(9, 6))

# Figure 1: Voltages vs. Time
for voltage in voltage_columns:
    ax1.plot(df['Time'], df[voltage], label=voltage, marker='o')
ax1.set_xlabel('Time')
ax1.set_ylabel('Voltage')
ax1.set_title('Voltage vs. Time')
ax1.legend()
ax1.grid(True)

# Figure 2: Create 3-row subplot for voltages, currents, and time step
# Plot Voltages vs. Time
for voltage in voltage_columns:
    axes[0].plot(df['Time'], df[voltage], label=voltage + ' vs Time', marker='o')
axes[0].set_ylabel('Voltage')
axes[0].set_title('Voltages vs. Time')
axes[0].legend()
axes[0].grid(True)

# Plot Currents vs. Time
for current in current_columns:
    axes[1].plot(df['Time'], df[current], label=current + ' vs Time', marker='x')
axes[1].set_ylabel('Current')
axes[1].set_title('Currents vs. Time')
axes[1].legend()
axes[1].grid(True)

# Plot Time Step vs. Time
axes[2].plot(df['Time'], df['Time Step'], label='Time Step vs Time', marker='o', color='green')
axes[2].set_xlabel('Time')
axes[2].set_ylabel('Time Step')
axes[2].set_title('Time Step vs. Time')
axes[2].legend()
axes[2].grid(True)

# Figure 3: Current vs. Voltage
# if voltage_columns and current_columns:
#     selected_voltage = voltage_columns[0]
#     selected_current = current_columns[0]
#     df_sorted = df.sort_values(by=selected_voltage)
#     ax3.plot(df_sorted[selected_voltage], df_sorted[selected_current], label=f'Current vs. {selected_voltage}', marker='^', color='purple')
# ax3.set_xlabel(selected_voltage)
# ax3.set_ylabel(selected_current)
# ax3.set_title('Current vs. Voltage')
# ax3.legend()
# ax3.grid(True)

# Figure 3: Selected Voltages vs. Time
#Specify the selected voltages by index or name here
selected_voltages = [voltage_columns[8]]  # Adjust indices according to your needs

for voltage in selected_voltages:
    ax3.plot(df['Time'], df[voltage], label=f'{voltage} vs Time', marker='^')
ax3.set_xlabel('Time')
ax3.set_ylabel('Voltage')
ax3.set_title('Selected Voltages vs. Time')
ax3.legend()
ax3.grid(True)

# Show all figures
plt.show()



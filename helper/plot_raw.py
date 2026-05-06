#!/usr/bin/env python3
"""Module to plot EESpice binary .raw files using matplotlib."""

import sys
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

# Add the current directory to sys.path to ensure we can import decode_raw if run from helper/
sys.path.append(str(Path(__file__).parent))

try:
    from decode_raw import read_raw_file
except ImportError:
    # If invoked from project root, we might need this
    from helper.decode_raw import read_raw_file

def plot_raw_file(raw_file_path, output_png_path):
    """
    Reads a .raw file and saves a plot of all variables against the scale variable to a PNG.
    """
    plots = read_raw_file(raw_file_path)
    if not plots:
        print(f"Warning: No plots found in {raw_file_path}")
        return

    # For simplicity, we plot the first plot block. 
    # Batch simulations might have multiple, but we'll focus on the first for this tool.
    header, variables, rows = plots[0]
    
    if not rows:
        print(f"Warning: No data points in {raw_file_path}")
        return

    # Convert to numpy array for easier slicing: (points, variables)
    data = np.array(rows)
    
    # Scale variable is typically at index 0
    scale_var = variables[0]
    scale_idx, scale_name, scale_type = scale_var
    
    x_data = data[:, scale_idx]
    
    plt.figure(figsize=(10, 6))
    
    # Plot each other variable
    for i in range(len(variables)):
        idx, name, typ = variables[i]
        if idx == scale_idx:
            continue
        
        # If complex, plot the magnitude
        y_data = data[:, idx]
        if np.iscomplexobj(y_data):
            y_data = np.abs(y_data)
            label = f"{name} (magnitude)"
        else:
            label = name
            
        plt.plot(x_data, y_data, label=label)

    plt.xlabel(f"{scale_name} ({scale_type})")
    plt.ylabel("Value")
    plt.title(f"EESpice Plot: {header.get('Plotname', 'Unknown')}")
    plt.grid(True)
    plt.legend()
    
    plt.savefig(output_png_path)
    plt.close()
    print(f"Plot saved to {output_png_path}")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python plot_raw.py <input.raw> <output.png>")
        sys.exit(1)
    
    plot_raw_file(sys.argv[1], sys.argv[2])

#!/bin/bash
set -e

# Setup Python virtual environment
echo "=== Setting up Python venv ==="
python3 -m venv venv
source venv/bin/activate

echo "=== Installing dependencies ==="
pip install -r requirements.txt

# Create Results directory if it doesn't exist
mkdir -p Results

echo "=== Running EEspice Simulation (via Docker) ==="
# We use binary format to generate the .raw file
./eespice.sh run -o Results/inverter.raw -f binary Netlist/Inverter.cir

echo "=== Plotting Results ==="
python helper/plot_raw.py Results/inverter.raw Results/inverter_plot.png

echo "=== Verifying Output ==="
if [ -f "Results/inverter_plot.png" ]; then
    echo "SUCCESS: Results/inverter_plot.png was created."
else
    echo "FAILURE: Results/inverter_plot.png was not created."
    exit 1
fi

# Test importing the module
echo "=== Testing module import ==="
python -c "from helper.plot_raw import plot_raw_file; print('Import successful')"

echo "=== Test Completed Successfully ==="
deactivate

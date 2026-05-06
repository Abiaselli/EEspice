#!/bin/bash
set -e

# EEspice Docker Test Script
# This script runs a test simulation using the eespice.sh helper
# to verify that the Docker image is working correctly and producing output.

NETLIST="Netlist/Inverter.cir"
OUTPUT_FILE="Netlist/Results/inverter.txt"

echo "=== Running EEspice Docker Test ==="

# 1. Clean up any existing output
if [ -f "$OUTPUT_FILE" ]; then
    echo "Removing existing $OUTPUT_FILE..."
    rm "$OUTPUT_FILE"
fi

# 2. Run the simulation
echo "Executing: ./eespice.sh run -o $OUTPUT_FILE -f ascii $NETLIST"
./eespice.sh run -o "$OUTPUT_FILE" -f ascii "$NETLIST"

# 3. Verify output
if [ -f "$OUTPUT_FILE" ]; then
    echo "SUCCESS: $OUTPUT_FILE was created."
    
    # Check ownership
    OWNER=$(ls -l "$OUTPUT_FILE" | awk '{print $3}')
    echo "File owner: $OWNER"
    
    if [ "$OWNER" == "root" ]; then
        echo "WARNING: File is owned by root. Check --user mapping in eespice.sh."
    else
        echo "SUCCESS: File is owned by current user ($OWNER)."
    fi

    # Clean up
    echo "Cleaning up test output..."
    rm "$OUTPUT_FILE"
else
    echo "FAILURE: $OUTPUT_FILE was not created."
    exit 1
fi

echo "=== Test Completed Successfully ==="

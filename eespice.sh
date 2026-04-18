#!/bin/bash

# EEspice Docker Helper Script
# This script simplifies building and running EEspice using Docker.

IMAGE_NAME="eespice"

show_help() {
    echo "Usage: ./eespice.sh [command] [arguments]"
    echo ""
    echo "Commands:"
    echo "  build         Build the EEspice Docker image"
    echo "  run [file]    Run a simulation (e.g., ./eespice.sh run Netlist/Inverter.cir)"
    echo "  shell         Open an interactive shell inside the container"
    echo "  help          Show this help message"
    echo ""
}

case "$1" in
    build)
        echo "Building EEspice Docker image..."
        docker build -t $IMAGE_NAME .
        ;;
    run)
        if [ -z "$2" ]; then
            echo "Error: No netlist file specified."
            show_help
            exit 1
        fi

        # Get absolute path of the file to handle mounting correctly
        FILE_PATH=$(realpath "$2")
        DIR_PATH=$(dirname "$FILE_PATH")
        FILE_NAME=$(basename "$FILE_PATH")

        echo "Running simulation: $FILE_NAME"
        # Mount the directory containing the netlist to /sim in the container
        # This ensures the output files are written back to the host directory
        docker run --rm -v "$DIR_PATH":/sim $IMAGE_NAME "$FILE_NAME"
        ;;
    shell)
        echo "Opening interactive shell..."
        docker run --rm -it --entrypoint /bin/bash -v "$(pwd)":/sim $IMAGE_NAME
        ;;
    help|*)
        show_help
        ;;
esac

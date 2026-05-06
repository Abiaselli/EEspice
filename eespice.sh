#!/bin/bash

# EEspice Docker Helper Script
# This script simplifies building and running EEspice using Docker.

IMAGE_NAME="eespice"

show_help() {
    echo "Usage: ./eespice.sh [command] [arguments]"
    echo ""
    echo "Commands:"
    echo "  build [options] Build the EEspice Docker image"
    echo "      --native            Compile with -march=native (best performance for this machine)"
    echo "      --native arm64      Compile for generic modern ARM64 (armv8-a)"
    echo "      --native x86        Compile for generic modern x86-64 (x86-64-v3)"
    echo "      (default)           Compile for generic baseline (x86-64 or armv8-a)"
    echo "  run [options] [file]  Run a simulation (e.g., ./eespice.sh run Netlist/Inverter.cir)"
    echo "      -o, --output <path>    Specify output file path"
    echo "      -f, --format <format>  Enforce output format (ascii or binary)"
    echo "      -p, --plot [path]      Generate a PNG plot from the results (binary format only)"
    echo "  shell         Open an interactive shell inside the container"

    echo "  help          Show this help message"
    echo ""
}

case "$1" in
    build)
        # Default to generic baseline
        if [[ "$(uname -m)" == "x86_64" ]]; then
            MARCH="x86-64"
        else
            MARCH="armv8-a"
        fi

        shift
        while [[ $# -gt 0 ]]; do
            case $1 in
                --native)
                    if [[ "$2" == "arm64" ]]; then
                        MARCH="armv8-a"
                        shift 2
                    elif [[ "$2" == "x86" ]]; then
                        MARCH="x86-64-v3"
                        shift 2
                    else
                        MARCH="native"
                        shift
                    fi
                    ;;
                *)
                    echo "Unknown build option: $1"
                    exit 1
                    ;;
            esac
        done

        echo "Building EEspice Docker image with MARCH=$MARCH..."
        docker build --build-arg MARCH=$MARCH -t $IMAGE_NAME .
        ;;
    run)
        shift
        EESPICE_ARGS=()
        DOCKER_MOUNTS=()
        NETLIST_FILE=""
        PLOT_PATH=""
        GENERATE_PLOT=false

        while [[ $# -gt 0 ]]; do
            case $1 in
                -o|--output)
                    OUT_PATH=$(realpath -m "$2")
                    OUT_DIR=$(dirname "$OUT_PATH")
                    DOCKER_MOUNTS+=("-v" "$OUT_DIR:$OUT_DIR")
                    EESPICE_ARGS+=("-o" "$OUT_PATH")
                    shift 2
                    ;;
                -f|--format)
                    EESPICE_ARGS+=("-f" "$2")
                    shift 2
                    ;;
                -p|--plot)
                    GENERATE_PLOT=true
                    if [[ $# -gt 1 && ! "$2" == -* ]]; then
                        PLOT_PATH="$2"
                        shift 2
                    else
                        shift
                    fi
                    ;;
                *)
                    if [ -z "$NETLIST_FILE" ]; then
                        NETLIST_FILE="$1"
                    else
                        echo "Error: Multiple netlist files specified or unknown argument: $1"
                        exit 1
                    fi
                    shift
                    ;;
            esac
        done

        if [ -z "$NETLIST_FILE" ]; then
            echo "Error: No netlist file specified."
            show_help
            exit 1
        fi

        # If plotting is requested without an explicit output path, ensure binary format is used
        # because the plotter requires binary .raw files.
        if [ "$GENERATE_PLOT" = true ]; then
            # Find if format was already specified
            FORMAT_SPECIFIED=false
            OUTPUT_SPECIFIED=false
            for arg in "${EESPICE_ARGS[@]}"; do
                if [[ "$arg" == "ascii" ]]; then
                    echo "Error: Plotting requires binary format, but ascii was specified."
                    exit 1
                fi
                if [[ "$arg" == "binary" ]]; then
                    FORMAT_SPECIFIED=true
                fi
                if [[ "$arg" == "-o" ]]; then
                    OUTPUT_SPECIFIED=true
                fi
            done
            if [ "$FORMAT_SPECIFIED" = false ]; then
                EESPICE_ARGS+=("-f" "binary")
            fi
            
            # If -o is missing, provide a default output path so the plotter knows where to look
            if [ "$OUTPUT_SPECIFIED" = false ]; then
                DEFAULT_OUT="Results/output.raw"
                mkdir -p Results
                OUT_PATH=$(realpath -m "$DEFAULT_OUT")
                OUT_DIR=$(dirname "$OUT_PATH")
                DOCKER_MOUNTS+=("-v" "$OUT_DIR:$OUT_DIR")
                EESPICE_ARGS+=("-o" "$OUT_PATH")
                echo "Notice: No output path specified. Defaulting to $DEFAULT_OUT for plotting."
            fi
        fi

        # Get absolute path of the file to handle mounting correctly
        FILE_PATH=$(realpath "$NETLIST_FILE")
        DIR_PATH=$(dirname "$FILE_PATH")
        FILE_NAME=$(basename "$FILE_PATH")

        echo "Running simulation: $FILE_NAME"
        # Mount the netlist's directory to /netlist so eespice can read it.
        # Mount the current working directory to /sim (the default WORKDIR in the container).
        # We also mount any custom output directories.
        docker run --rm --user $(id -u):$(id -g) \
            -v "$DIR_PATH":/netlist \
            -v "$(pwd)":/sim \
            "${DOCKER_MOUNTS[@]}" \
            $IMAGE_NAME "${EESPICE_ARGS[@]}" "/netlist/$FILE_NAME"

        # Handle plotting after simulation
        if [ "$GENERATE_PLOT" = true ]; then
            # We need to find where the .raw file was saved.
            # 1. Check if -o was provided
            RAW_FILE=""
            for ((i=0; i<${#EESPICE_ARGS[@]}; i++)); do
                if [[ "${EESPICE_ARGS[$i]}" == "-o" ]]; then
                    RAW_FILE="${EESPICE_ARGS[$((i+1))]}"
                fi
            done

            # 2. If no -o, we'd need to parse the netlist for .output, which is complex for a bash script.
            # For now, we expect -o to be used with --plot for reliability, or we can warn.
            if [ -z "$RAW_FILE" ]; then
                 echo "Warning: No output path specified with -o. Plotting might fail if .raw file location is unknown."
                 # Try to guess or parse? Let's look for Results/ directory as a fallback if it exists.
            fi

            if [ -n "$RAW_FILE" ] && [ -f "$RAW_FILE" ]; then
                if [ -z "$PLOT_PATH" ]; then
                    PLOT_PATH="${RAW_FILE%.*}.png"
                fi
                
                # Check if venv exists and use it, otherwise use system python
                PYTHON_CMD="python3"
                if [ -d "venv" ]; then
                    PYTHON_CMD="./venv/bin/python3"
                fi
                
                echo "Generating plot: $PLOT_PATH"
                $PYTHON_CMD helper/plot_raw.py "$RAW_FILE" "$PLOT_PATH"
            else
                echo "Error: Could not find output file for plotting. Use -o to specify the output path."
            fi
        fi
        ;;
    shell)
        echo "Opening interactive shell..."
        docker run --rm -it --user $(id -u):$(id -g) --entrypoint /bin/bash -v "$(pwd)":/sim $IMAGE_NAME
        ;;
    help|*)
        show_help
        ;;
esac

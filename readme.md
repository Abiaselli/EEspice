# EEspice

EEspice is a SPICE-compatible circuit simulator optimized for performance. It features BSIM4 transistor model support, multithreaded simulation, and advanced matrix solver integration.

## Docker Usage (Recommended)

The easiest way to build and run EEspice is using Docker. This ensures all dependencies (SuiteSparse, Armadillo, OpenBLAS) are correctly configured without modifying your host system.

### 1. Build the Image
Build the persistent Docker image once (or whenever you update the source code). By default, it builds for a generic baseline architecture for maximum compatibility.

```bash
# Build with generic baseline (default)
./eespice.sh build

# Build optimized for the current machine
./eespice.sh build --native

# Build for generic modern x86 (e.g., AVX2 support)
./eespice.sh build --native x86

# Build for generic modern ARM64
./eespice.sh build --native arm64
```

### 2. Run a Simulation
Run a simulation by providing the path to a netlist file. You can control the output location and format using command-line flags.

```bash
# Basic run (defaults to binary tran_solution.raw)
./eespice.sh run Netlist/Inverter.cir

# Specify custom output path and format
./eespice.sh run -o results.txt -f ascii Netlist/Inverter.cir

# Specify custom output path (auto-detects format from extension)
./eespice.sh run -o results.raw Netlist/Inverter.cir
```

**Options:**
- `-o, --output <path>`: Specify the destination file path.
- `-f, --format <ascii|binary>`: Enforce the output format.
  - `binary`: Produces an `ngspice`-compatible binary RAW file (default).
  - `ascii`: Produces a human-readable CSV file (or a directory of CSVs in batch mode).

**Format Discovery:**
If `-f` is not provided, the format is determined by the file extension of the `-o` path:
- `.txt` or `.csv` &rarr; **ASCII**
- Anything else (or no extension) &rarr; **Binary**

*Note: The legacy `.output` directive inside netlists is now ignored in favor of these command-line flags.*

### 3. Batch Simulation
EEspice supports batch simulations by specifying ranges or lists for component values in the netlist (e.g., `R1 1 0 [1k:1k:10k]` or `V1 1 0 (1.2 1.5 1.8)`).

**Batch Output Behavior:**
- **Binary (Default)**: All simulation runs are saved into a **single multi-plot `.raw` file**. This is the recommended mode for large batches.
- **ASCII**: If you enforce `-f ascii` or use a `.csv` extension, results are saved as a **directory of CSV files** (one file per run).

```bash
# Run batch simulation and save to a single binary file
./eespice.sh run -o sweep_results.raw Netlist/batch.cir

# Run batch simulation and save to a directory of CSVs
./eespice.sh run -o sweep_folder -f ascii Netlist/batch.cir
```

### 4. Interactive Shell
If you need to explore the environment inside the container:
```bash
./eespice.sh shell
```

---

## Native Installation (Manual)

If you prefer to build EEspice directly on your host, ensure you have the following dependencies installed.

### Dependencies
On Ubuntu/Debian:
```bash
sudo apt-get update
sudo apt-get install -y \
    build-essential \
    cmake \
    pkg-config \
    libarmadillo-dev \
    libsuitesparse-dev \
    libsuperlu-dev \
    libopenblas-dev
```

### Build Instructions
```bash
cmake -S . -B build -DEESPICE_BLAS_BACKEND=OpenBLAS -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)
```

The resulting binary will be located at `./build/eespice`.

## Project Structure
- `src/`: Core simulator source code.
- `Netlist/`: Sample circuit netlists for testing.
- `Dockerfile`: Multi-stage build configuration.
- `eespice.sh`: Helper script for Docker operations.

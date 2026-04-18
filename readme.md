# EEspice

EEspice is a SPICE-compatible circuit simulator optimized for performance. It features BSIM4 transistor model support, multithreaded simulation, and advanced matrix solver integration.

## Docker Usage (Recommended)

The easiest way to build and run EEspice is using Docker. This ensures all dependencies (SuiteSparse, Armadillo, OpenBLAS) are correctly configured without modifying your host system.

### 1. Build the Image
Build the persistent Docker image once (or whenever you update the source code):
```bash
./eespice.sh build
```

### 2. Run a Simulation
Run a simulation by providing the path to a netlist file. The script automatically mounts the necessary directories so that output files (e.g., `.csv` results) are saved back to your host machine.
```bash
./eespice.sh run Netlist/Inverter.cir
```

### 3. Interactive Shell
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

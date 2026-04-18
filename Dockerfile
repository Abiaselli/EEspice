# Stage 1: Build environment
FROM ubuntu:24.04 AS builder

ENV DEBIAN_FRONTEND=noninteractive

# Install build-time dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    pkg-config \
    libarmadillo-dev \
    libsuitesparse-dev \
    libsuperlu-dev \
    libopenblas-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app
COPY . .

# Build the executable
RUN cmake -S . -B build \
    -DEESPICE_BLAS_BACKEND=OpenBLAS \
    -DCMAKE_BUILD_TYPE=Release

RUN cmake --build build -j$(nproc)

# Stage 2: Runtime environment
FROM ubuntu:24.04 AS runtime

ENV DEBIAN_FRONTEND=noninteractive

# Install only the runtime libraries (much smaller than build-time)
# Ubuntu 24.04 uses these specific versions
RUN apt-get update && apt-get install -y \
    libarmadillo12 \
    libsuitesparse-dev \
    libsuperlu6 \
    libopenblas0 \
    libgomp1 \
    && rm -rf /var/lib/apt/lists/*

# Copy the binary from the builder stage
COPY --from=builder /app/build/eespice /usr/local/bin/eespice

# Set the working directory to where the user will mount their netlists
WORKDIR /sim

# Default command
ENTRYPOINT ["eespice"]

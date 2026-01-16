/*
 * Benchmark: Dense vs Sparse Matrix Copy Performance
 *
 * This program compares the copy speed of dense and sparse matrices using:
 * 1. Assignment operator (=)
 * 2. std::memcpy
 *
 * Dense matrices: Compares full matrix copy (n*n elements)
 * Sparse matrices: Compares:
 *   - Assignment (=): Copies entire CSC structure (values + row_indices + col_ptrs)
 *   - memcpy: Copies only the values array (assumes same sparsity pattern)
 *
 * Configuration:
 * - Matrix sizes: 100, 500, 1000, 5000
 * - Sparsity: 1% (typical for circuit matrices)
 * - Iterations: 1000 per test
 *
 * Expected results:
 * - Dense: memcpy ~= assignment (both optimized)
 * - Sparse: memcpy 2-7x faster (only copies nnz values vs entire structure)
 */

#include <iostream>
#include <armadillo>
#include <chrono>
#include <cstring>
#include <iomanip>
#include <vector>

template<typename Func>
double measure_time_ms(Func&& func, int iterations) {
    auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < iterations; i++) {
        func();
    }
    auto end = std::chrono::high_resolution_clock::now();
    return std::chrono::duration<double, std::milli>(end - start).count() / iterations;
}

void benchmark_dense(size_t n, int iterations) {
    arma::mat src(n, n, arma::fill::randu);
    arma::mat dst(n, n);

    // Warmup
    dst = src;
    std::memcpy(dst.memptr(), src.memptr(), n * n * sizeof(double));

    // Method 1: Assignment operator
    double time_assign = measure_time_ms([&]() { dst = src; }, iterations);

    // Method 2: memcpy
    double time_memcpy = measure_time_ms([&]() {
        std::memcpy(dst.memptr(), src.memptr(), n * n * sizeof(double));
    }, iterations);

    double speedup = time_assign / time_memcpy;

    std::cout << "Size: " << n << "x" << n << std::endl;
    std::cout << "  Assignment (=):  " << std::fixed << std::setprecision(4) << time_assign << " ms" << std::endl;
    std::cout << "  memcpy:          " << std::fixed << std::setprecision(4) << time_memcpy << " ms" << std::endl;
    std::cout << "  Speedup:         " << std::fixed << std::setprecision(2) << speedup << "x" << std::endl;
    std::cout << std::endl;
}

void benchmark_sparse(size_t n, double sparsity, int iterations) {
    arma::sp_mat src = arma::sprandu<arma::sp_mat>(n, n, sparsity);

    // Create destination with same sparsity pattern
    arma::sp_mat dst = src;

    // Warmup
    dst = src;
    src.sync();
    dst.sync();
    std::memcpy(const_cast<double*>(dst.values), src.values, src.n_nonzero * sizeof(double));

    // Method 1: Assignment operator
    double time_assign = measure_time_ms([&]() { dst = src; }, iterations);

    // Method 2: memcpy (values only)
    double time_memcpy = measure_time_ms([&]() {
        src.sync();
        dst.sync();
        std::memcpy(const_cast<double*>(dst.values), src.values, src.n_nonzero * sizeof(double));
    }, iterations);

    double speedup = time_assign / time_memcpy;

    std::cout << "Size: " << n << "x" << n << " (nnz: " << src.n_nonzero << ")" << std::endl;
    std::cout << "  Assignment (=):  " << std::fixed << std::setprecision(4) << time_assign << " ms" << std::endl;
    std::cout << "  memcpy (values): " << std::fixed << std::setprecision(4) << time_memcpy << " ms" << std::endl;
    std::cout << "  Speedup:         " << std::fixed << std::setprecision(2) << speedup << "x" << std::endl;
    std::cout << std::endl;
}

int main() {
    std::vector<size_t> sizes = {100, 500, 1000, 5000};
    double sparsity = 0.01;  // 1%
    int iterations = 1000;

    std::cout << "=== Dense Matrix Copy Benchmark ===" << std::endl;
    std::cout << "Iterations: " << iterations << std::endl << std::endl;
    for (size_t n : sizes) {
        benchmark_dense(n, iterations);
    }

    std::cout << "=== Sparse Matrix Copy Benchmark ===" << std::endl;
    std::cout << "Iterations: " << iterations << ", Sparsity: " << (sparsity * 100) << "%" << std::endl << std::endl;
    for (size_t n : sizes) {
        benchmark_sparse(n, sparsity, iterations);
    }

    return 0;
}

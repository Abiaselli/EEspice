// file: arrayVSvec.cpp
// Benchmark comparing std::array (stack) vs std::vector (heap) performance
// Build: g++ -O3 -march=native -std=c++20 arrayVSvec.cpp -o arrayVSvec
// Run:   ./arrayVSvec

#include <iostream>
#include <vector>
#include <array>
#include <chrono>
#include <random>
#include <iomanip>
#include <numeric>
#include <algorithm>
#include <cstdint>

using namespace std::chrono;

constexpr double MIN_TEST_DURATION = 1.0; // seconds

// Prevent compiler from optimizing away computations
template<typename T>
void do_not_optimize(T const& value) {
    asm volatile("" : : "r,m"(value) : "memory");
}

// Benchmark with auto-calibrated iterations to run for at least MIN_TEST_DURATION seconds
template<typename F>
std::pair<double, uint64_t> benchmark_auto_iterations(F&& func) {
    uint64_t iterations = 0;
    auto start = high_resolution_clock::now();
    auto elapsed = 0.0;

    // Run until we reach minimum duration
    while (elapsed < MIN_TEST_DURATION) {
        func();
        iterations++;

        auto now = high_resolution_clock::now();
        elapsed = duration<double>(now - start).count();
    }

    return {elapsed, iterations};
}

// Sequential write benchmark
template<typename Container>
std::pair<double, uint64_t> benchmark_sequential_write(Container& container) {
    return benchmark_auto_iterations([&]() {
        for (size_t i = 0; i < container.size(); ++i) {
            container[i] = static_cast<double>(i * 3.14159);
        }
        do_not_optimize(container.data());
    });
}

// Sequential read benchmark
template<typename Container>
std::pair<double, uint64_t> benchmark_sequential_read(const Container& container) {
    return benchmark_auto_iterations([&]() {
        double sum = 0.0;
        for (size_t i = 0; i < container.size(); ++i) {
            sum += container[i];
        }
        do_not_optimize(sum);
    });
}

// Random access write benchmark
template<typename Container>
std::pair<double, uint64_t> benchmark_random_write(Container& container, const std::vector<size_t>& indices) {
    return benchmark_auto_iterations([&]() {
        for (size_t idx : indices) {
            container[idx] = static_cast<double>(idx * 2.71828);
        }
        do_not_optimize(container.data());
    });
}

// Random access read benchmark
template<typename Container>
std::pair<double, uint64_t> benchmark_random_read(const Container& container, const std::vector<size_t>& indices) {
    return benchmark_auto_iterations([&]() {
        double sum = 0.0;
        for (size_t idx : indices) {
            sum += container[idx];
        }
        do_not_optimize(sum);
    });
}

// Generate random indices for testing
std::vector<size_t> generate_random_indices(size_t container_size, size_t num_accesses) {
    std::vector<size_t> indices(num_accesses);
    std::mt19937_64 rng(42); // Fixed seed for reproducibility
    std::uniform_int_distribution<size_t> dist(0, container_size - 1);

    for (auto& idx : indices) {
        idx = dist(rng);
    }
    return indices;
}

void print_result(const char* test_name, size_t N,
                  double arr_time, uint64_t arr_iters,
                  double vec_time, uint64_t vec_iters) {

    double arr_ops_per_sec = arr_iters / arr_time;
    double vec_ops_per_sec = vec_iters / vec_time;

    double arr_bandwidth = (arr_iters * N * sizeof(double)) / arr_time / 1e9;
    double vec_bandwidth = (vec_iters * N * sizeof(double)) / vec_time / 1e9;

    double speedup = vec_time / arr_time;

    std::cout << "\n" << test_name << ":\n";
    std::cout << "  std::array:  " << std::setw(10) << arr_iters << " iterations in "
              << std::setw(6) << std::setprecision(3) << arr_time << " s  ->  "
              << std::setw(8) << std::setprecision(2) << arr_bandwidth << " GB/s\n";
    std::cout << "  std::vector: " << std::setw(10) << vec_iters << " iterations in "
              << std::setw(6) << std::setprecision(3) << vec_time << " s  ->  "
              << std::setw(8) << std::setprecision(2) << vec_bandwidth << " GB/s\n";
    std::cout << "  Speedup: " << std::setw(6) << std::setprecision(4) << speedup << "x";

    if (speedup > 1.02) {
        std::cout << "  [array FASTER]";
    } else if (speedup < 0.98) {
        std::cout << "  [vector FASTER]";
    } else {
        std::cout << "  [virtually EQUAL]";
    }
    std::cout << "\n";
}

void print_result_no_bandwidth(const char* test_name,
                               double arr_time, uint64_t arr_iters,
                               double vec_time, uint64_t vec_iters) {

    double arr_ops_per_sec = arr_iters / arr_time;
    double vec_ops_per_sec = vec_iters / vec_time;

    double speedup = vec_time / arr_time;

    std::cout << "\n" << test_name << ":\n";
    std::cout << "  std::array:  " << std::setw(10) << arr_iters << " iterations in "
              << std::setw(6) << std::setprecision(3) << arr_time << " s  ->  "
              << std::setw(10) << std::setprecision(0) << arr_ops_per_sec << " ops/s\n";
    std::cout << "  std::vector: " << std::setw(10) << vec_iters << " iterations in "
              << std::setw(6) << std::setprecision(3) << vec_time << " s  ->  "
              << std::setw(10) << std::setprecision(0) << vec_ops_per_sec << " ops/s\n";
    std::cout << "  Speedup: " << std::setw(6) << std::setprecision(4) << speedup << "x";

    if (speedup > 1.02) {
        std::cout << "  [array FASTER]";
    } else if (speedup < 0.98) {
        std::cout << "  [vector FASTER]";
    } else {
        std::cout << "  [virtually EQUAL]";
    }
    std::cout << "\n";
}

template<size_t N>
void run_benchmark(const char* size_label) {
    constexpr size_t num_random_accesses = std::min(N, size_t(10000));

    std::cout << "\n" << std::string(80, '=') << "\n";
    std::cout << "Benchmark: " << size_label << " (" << N << " elements, "
              << std::setprecision(2) << std::fixed << (N * sizeof(double) / 1024.0) << " KB)\n";
    std::cout << std::string(80, '=') << "\n";

    // Prepare containers
    std::array<double, N> arr;
    std::vector<double> vec(N);

    // Initialize both containers
    arr.fill(0.0);
    std::fill(vec.begin(), vec.end(), 0.0);

    // Generate random indices
    auto random_indices = generate_random_indices(N, num_random_accesses);

    std::cout << "\nRunning tests (each test runs for " << MIN_TEST_DURATION << " second minimum)...\n";

    // Sequential Write
    std::cout << "\n[1/4] Sequential Write...";
    std::cout.flush();
    auto [arr_seq_write_time, arr_seq_write_iters] = benchmark_sequential_write(arr);
    auto [vec_seq_write_time, vec_seq_write_iters] = benchmark_sequential_write(vec);
    print_result("Sequential Write", N, arr_seq_write_time, arr_seq_write_iters, vec_seq_write_time, vec_seq_write_iters);

    // Sequential Read
    std::cout << "\n[2/4] Sequential Read...";
    std::cout.flush();
    auto [arr_seq_read_time, arr_seq_read_iters] = benchmark_sequential_read(arr);
    auto [vec_seq_read_time, vec_seq_read_iters] = benchmark_sequential_read(vec);
    print_result("Sequential Read", N, arr_seq_read_time, arr_seq_read_iters, vec_seq_read_time, vec_seq_read_iters);

    // Random Write
    std::cout << "\n[3/4] Random Write (" << num_random_accesses << " accesses per iteration)...";
    std::cout.flush();
    auto [arr_rand_write_time, arr_rand_write_iters] = benchmark_random_write(arr, random_indices);
    auto [vec_rand_write_time, vec_rand_write_iters] = benchmark_random_write(vec, random_indices);
    print_result_no_bandwidth("Random Write", arr_rand_write_time, arr_rand_write_iters, vec_rand_write_time, vec_rand_write_iters);

    // Random Read
    std::cout << "\n[4/4] Random Read (" << num_random_accesses << " accesses per iteration)...";
    std::cout.flush();
    auto [arr_rand_read_time, arr_rand_read_iters] = benchmark_random_read(arr, random_indices);
    auto [vec_rand_read_time, vec_rand_read_iters] = benchmark_random_read(vec, random_indices);
    print_result_no_bandwidth("Random Read", arr_rand_read_time, arr_rand_read_iters, vec_rand_read_time, vec_rand_read_iters);
}

int main() {
    std::cout << "\n";
    std::cout << "╔════════════════════════════════════════════════════════════════════════════╗\n";
    std::cout << "║   std::array (stack) vs std::vector (heap) Performance Benchmark          ║\n";
    std::cout << "╚════════════════════════════════════════════════════════════════════════════╝\n";
    std::cout << "\nHypothesis: Stack-allocated arrays should be faster than heap-allocated vectors\n";
    std::cout << "Test configuration: double precision (8 bytes per element)\n";
    std::cout << "Each test runs for minimum " << MIN_TEST_DURATION << " second(s) for statistical accuracy\n";

    // Small size - fits in L1 cache (typically 32-64 KB)
    run_benchmark<100>("Tiny - L1 Cache");

    // Medium-small - still in L1/L2
    run_benchmark<1000>("Small - L1/L2 Cache");

    // Medium size - fits in L2/L3 cache
    run_benchmark<10000>("Medium - L2/L3 Cache");

    // Large size - L3 cache boundary
    run_benchmark<100000>("Large - L3/Memory");

    // Very large - definitely main memory bandwidth limited
    // Note: 1M elements = 8MB, will likely cause stack overflow for std::array
    // Uncomment to test (will likely crash):
    // run_benchmark<1000000>("Very Large - Main Memory");

    std::cout << "\n" << std::string(80, '=') << "\n";
    std::cout << "SUMMARY\n";
    std::cout << std::string(80, '=') << "\n";
    std::cout << "Speedup > 1.0 means std::array is faster\n";
    std::cout << "Speedup < 1.0 means std::vector is faster\n";
    std::cout << "\nKey findings:\n";
    std::cout << "- Small arrays (L1 cache): Expected to show advantage for std::array\n";
    std::cout << "- Large arrays: Performance should converge (both memory-bandwidth limited)\n";
    std::cout << "- Very large std::array sizes WILL cause stack overflow (practical limitation)\n";
    std::cout << "\nConclusion:\n";
    std::cout << "For small, fixed-size data, std::array may have slight advantages.\n";
    std::cout << "For large or variable-size data, std::vector is safer and equally fast.\n";
    std::cout << std::string(80, '=') << "\n\n";

    return 0;
}

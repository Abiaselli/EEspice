#include <armadillo>
#include <array>
#include <chrono>
#include <iostream>
#include <iomanip>
#include <random>
#include <vector>
#include <cmath>

inline void Stamp(arma::mat &mat, int row, int col, double val)
{
    mat(row, col) += val;
}

inline void Stamp_index(arma::mat &mat, int row, int col, double val)
{
    double* p = mat.memptr();
    p[row + col * mat.n_rows] += val;
}

// Verification test: check if both functions produce identical results
bool verify_correctness(int size = 10) {
    std::cout << "\n=== Correctness Verification (size=" << size << ") ===" << std::endl;

    arma::mat mat1(size, size, arma::fill::zeros);
    arma::mat mat2(size, size, arma::fill::zeros);

    std::mt19937 rng(42);
    std::uniform_int_distribution<int> dist(0, size - 1);
    std::uniform_real_distribution<double> val_dist(-10.0, 10.0);

    // Perform 1000 random stamps
    std::vector<std::tuple<int, int, double>> operations;
    for (int i = 0; i < 1000; ++i) {
        int row = dist(rng);
        int col = dist(rng);
        double val = val_dist(rng);
        operations.push_back({row, col, val});
    }

    // Apply using both methods
    for (const auto& [row, col, val] : operations) {
        Stamp(mat1, row, col, val);
        Stamp_index(mat2, row, col, val);
    }

    // Compare results
    double max_diff = 0.0;
    for (arma::uword i = 0; i < mat1.n_elem; ++i) {
        double diff = std::abs(mat1(i) - mat2(i));
        max_diff = std::max(max_diff, diff);
    }

    std::cout << "Maximum difference: " << std::scientific << max_diff << std::endl;

    bool passed = (max_diff < 1e-10);
    std::cout << "Verification: " << (passed ? "PASSED ✓" : "FAILED ✗") << std::endl;

    return passed;
}

// Statistics helper
struct Stats {
    double mean;
    double stddev;
    double min;
    double max;

    Stats(const std::vector<double>& times) {
        mean = std::accumulate(times.begin(), times.end(), 0.0) / times.size();

        double sq_sum = 0.0;
        min = times[0];
        max = times[0];
        for (double t : times) {
            sq_sum += (t - mean) * (t - mean);
            min = std::min(min, t);
            max = std::max(max, t);
        }
        stddev = std::sqrt(sq_sum / times.size());
    }
};

// Benchmark template
template<typename StampFunc>
Stats benchmark_stamp(StampFunc stamp_func, int matrix_size, int num_operations, const std::string& name) {
    const int num_runs = 20;
    std::vector<double> times;

    std::mt19937 rng(12345);
    std::uniform_int_distribution<int> dist(0, matrix_size - 1);
    std::uniform_real_distribution<double> val_dist(-1.0, 1.0);

    // Pre-generate random operations
    std::vector<std::tuple<int, int, double>> operations;
    for (int i = 0; i < num_operations; ++i) {
        operations.push_back({dist(rng), dist(rng), val_dist(rng)});
    }

    // Warm-up run
    arma::mat mat(matrix_size, matrix_size, arma::fill::zeros);
    for (const auto& [row, col, val] : operations) {
        stamp_func(mat, row, col, val);
    }

    // Timed runs
    for (int run = 0; run < num_runs; ++run) {
        arma::mat mat(matrix_size, matrix_size, arma::fill::zeros);

        auto start = std::chrono::high_resolution_clock::now();
        for (const auto& [row, col, val] : operations) {
            stamp_func(mat, row, col, val);
        }
        auto end = std::chrono::high_resolution_clock::now();

        std::chrono::duration<double, std::micro> duration = end - start;
        times.push_back(duration.count());
    }

    return Stats(times);
}

void print_comparison(const std::string& test_name, int size, int ops,
                     const Stats& stats_regular, const Stats& stats_index) {
    std::cout << "\n--- " << test_name << " (size=" << size << ", ops=" << ops << ") ---" << std::endl;
    std::cout << std::fixed << std::setprecision(2);

    std::cout << "\nStamp (regular):" << std::endl;
    std::cout << "  Mean: " << stats_regular.mean << " μs" << std::endl;
    std::cout << "  Std:  " << stats_regular.stddev << " μs" << std::endl;
    std::cout << "  Min:  " << stats_regular.min << " μs" << std::endl;
    std::cout << "  Max:  " << stats_regular.max << " μs" << std::endl;

    std::cout << "\nStamp_index:" << std::endl;
    std::cout << "  Mean: " << stats_index.mean << " μs" << std::endl;
    std::cout << "  Std:  " << stats_index.stddev << " μs" << std::endl;
    std::cout << "  Min:  " << stats_index.min << " μs" << std::endl;
    std::cout << "  Max:  " << stats_index.max << " μs" << std::endl;

    double speedup = stats_regular.mean / stats_index.mean;
    double percent_change = ((stats_index.mean - stats_regular.mean) / stats_regular.mean) * 100.0;

    std::cout << "\nSpeedup: " << std::setprecision(3) << speedup << "x ";
    if (speedup > 1.0) {
        std::cout << "(Stamp_index is " << std::setprecision(1) << (speedup - 1.0) * 100.0 << "% faster)";
    } else {
        std::cout << "(Stamp_index is " << std::setprecision(1) << std::abs(percent_change) << "% slower)";
    }
    std::cout << std::endl;

    double ns_per_op_regular = (stats_regular.mean * 1000.0) / ops;
    double ns_per_op_index = (stats_index.mean * 1000.0) / ops;
    std::cout << "Time per operation: " << std::setprecision(2)
              << ns_per_op_regular << " ns (Stamp) vs "
              << ns_per_op_index << " ns (Stamp_index)" << std::endl;
}

// Benchmark with sequential access pattern (similar to BSIM4 stamping)
void benchmark_sequential_pattern(int size, int num_ops) {
    auto stamp_lambda = [](arma::mat& mat, int row, int col, double val) {
        Stamp(mat, row, col, val);
    };
    auto stamp_index_lambda = [](arma::mat& mat, int row, int col, double val) {
        Stamp_index(mat, row, col, val);
    };

    Stats stats_regular = benchmark_stamp(stamp_lambda, size, num_ops, "Sequential");
    Stats stats_index = benchmark_stamp(stamp_index_lambda, size, num_ops, "Sequential");

    print_comparison("Sequential Access Pattern", size, num_ops, stats_regular, stats_index);
}

// Benchmark realistic BSIM4-like stamping pattern
void benchmark_bsim4_pattern(int size, int num_devices) {
    const int num_runs = 20;
    std::vector<double> times_regular, times_index;

    std::mt19937 rng(54321);
    std::uniform_int_distribution<int> node_dist(0, size - 1);

    // Simulate BSIM4 node structure for each device
    struct DeviceNodes {
        int d, g, s, b, dp, gp, sp, bp;
    };

    std::vector<DeviceNodes> devices;
    for (int i = 0; i < num_devices; ++i) {
        DeviceNodes nodes;
        nodes.d = node_dist(rng);
        nodes.g = node_dist(rng);
        nodes.s = node_dist(rng);
        nodes.b = node_dist(rng);
        nodes.dp = node_dist(rng);
        nodes.gp = node_dist(rng);
        nodes.sp = node_dist(rng);
        nodes.bp = node_dist(rng);
        devices.push_back(nodes);
    }

    std::uniform_real_distribution<double> val_dist(-1.0, 1.0);

    // Benchmark Stamp
    for (int run = 0; run < num_runs; ++run) {
        arma::mat mat(size, size, arma::fill::zeros);

        auto start = std::chrono::high_resolution_clock::now();
        for (const auto& dev : devices) {
            // Simulate typical BSIM4 stamping pattern (8 stamps per device)
            Stamp(mat, dev.dp, dev.d, val_dist(rng));
            Stamp(mat, dev.dp, dev.dp, val_dist(rng));
            Stamp(mat, dev.dp, dev.gp, val_dist(rng));
            Stamp(mat, dev.dp, dev.sp, val_dist(rng));
            Stamp(mat, dev.gp, dev.gp, val_dist(rng));
            Stamp(mat, dev.sp, dev.sp, val_dist(rng));
            Stamp(mat, dev.bp, dev.bp, val_dist(rng));
            Stamp(mat, dev.d, dev.dp, val_dist(rng));
        }
        auto end = std::chrono::high_resolution_clock::now();

        std::chrono::duration<double, std::micro> duration = end - start;
        times_regular.push_back(duration.count());
    }

    // Benchmark Stamp_index
    for (int run = 0; run < num_runs; ++run) {
        arma::mat mat(size, size, arma::fill::zeros);

        auto start = std::chrono::high_resolution_clock::now();
        for (const auto& dev : devices) {
            Stamp_index(mat, dev.dp, dev.d, val_dist(rng));
            Stamp_index(mat, dev.dp, dev.dp, val_dist(rng));
            Stamp_index(mat, dev.dp, dev.gp, val_dist(rng));
            Stamp_index(mat, dev.dp, dev.sp, val_dist(rng));
            Stamp_index(mat, dev.gp, dev.gp, val_dist(rng));
            Stamp_index(mat, dev.sp, dev.sp, val_dist(rng));
            Stamp_index(mat, dev.bp, dev.bp, val_dist(rng));
            Stamp_index(mat, dev.d, dev.dp, val_dist(rng));
        }
        auto end = std::chrono::high_resolution_clock::now();

        std::chrono::duration<double, std::micro> duration = end - start;
        times_index.push_back(duration.count());
    }

    Stats stats_regular(times_regular);
    Stats stats_index(times_index);

    int total_ops = num_devices * 8;
    print_comparison("BSIM4-like Pattern", size, total_ops, stats_regular, stats_index);
}

int main(){
    std::cout << "========================================" << std::endl;
    std::cout << "  Matrix Stamp Function Benchmark" << std::endl;
    std::cout << "========================================" << std::endl;

    // Correctness verification
    if (!verify_correctness(10)) {
        std::cerr << "ERROR: Correctness verification failed!" << std::endl;
        return 1;
    }
    if (!verify_correctness(100)) {
        std::cerr << "ERROR: Correctness verification failed!" << std::endl;
        return 1;
    }

    std::cout << "\n========================================" << std::endl;
    std::cout << "  Performance Benchmarks" << std::endl;
    std::cout << "========================================" << std::endl;

    // Test different matrix sizes
    std::cout << "\n### Random Access Pattern ###" << std::endl;
    benchmark_sequential_pattern(50, 50000);
    benchmark_sequential_pattern(100, 50000);
    benchmark_sequential_pattern(500, 50000);
    benchmark_sequential_pattern(1000, 50000);

    // BSIM4-like pattern (more realistic)
    std::cout << "\n\n### BSIM4-like Stamping Pattern ###" << std::endl;
    benchmark_bsim4_pattern(100, 1000);   // 100x100 matrix, 1000 devices
    benchmark_bsim4_pattern(500, 5000);   // 500x500 matrix, 5000 devices
    benchmark_bsim4_pattern(1000, 10000); // 1000x1000 matrix, 10000 devices

    std::cout << "\n========================================" << std::endl;
    std::cout << "  Benchmark Complete" << std::endl;
    std::cout << "========================================" << std::endl;

    return 0;
}
#define ARMA_DONT_USE_WRAPPER
#define ARMA_USE_MKL_ALLOC
#include <iostream>
#include <chrono>
#include <vector>
#include <iomanip>
#include <memory>
#include <algorithm>
#include <numeric>
#include <random>
#include <fstream>
#include <omp.h>

// EEspice includes
#include "CKT.hpp"
#include "global.hpp"
#include "hybrid_matrix.hpp"
#include "bsim4v82/bsim4v82setup.hpp"
#include "bsim4v82/bsim4v82temp.hpp"
#include "bsim4v82/bsim4v82load/bsim4v82applyStamps.hpp"

// Benchmarking includes
#include "loadomp.hpp"
#include "single.hpp"


std::shared_ptr<bsim4::BSIM4model> CreateBSIM4Model(){
    auto model = std::make_shared<bsim4::BSIM4model>();
    // Setup the model using the default temperature
    bsim4::modelSetup(*model, nomTemp);
    // Setup the model using the actual temperature
    bsim4::modelTemp(*model, nomTemp);
    return model;
}

void CreateNodes(std::vector<BSIM4> &bsim4){
    // Create nodes for the circuit - each instance gets unique node
    for (int i = 0; i < bsim4.size(); ++i) {
        const int node = i + 1;
        bsim4[i].bsim4v82Instance.BSIM4dNode = node;
        bsim4[i].bsim4v82Instance.BSIM4gNodeExt = node;
        bsim4[i].bsim4v82Instance.BSIM4sNode = node;
        bsim4[i].bsim4v82Instance.BSIM4bNode = node;
        bsim4[i].bsim4v82Instance.BSIM4dNodePrime = node;
        bsim4[i].bsim4v82Instance.BSIM4gNodePrime = node;
        bsim4[i].bsim4v82Instance.BSIM4gNodeMid = node;
        bsim4[i].bsim4v82Instance.BSIM4sNodePrime = node;
        bsim4[i].bsim4v82Instance.BSIM4bNodePrime = node;
        bsim4[i].bsim4v82Instance.BSIM4dbNode = node;
        bsim4[i].bsim4v82Instance.BSIM4sbNode = node;
        bsim4[i].bsim4v82Instance.BSIM4qNode = node;
    }
}

std::vector<BSIM4> CreateBSIM4Instances(std::shared_ptr<bsim4::BSIM4model> &model, size_t num_instances){
    std::vector<BSIM4> bsim4(num_instances);
    CKTcircuit ckt; // Dummy circuit for instance setup
    CreateNodes(bsim4); // Create nodes for each instance
    for(auto &b4 : bsim4){
        // Assign the model pointer
        b4.bsim4v82Instance.BSIM4modPtr = model;
        // Setup the instance using the default temperature
        bsim4::instanceSetup(*model, b4.bsim4v82Instance, ckt);
        // Setup the instance using the actual temperature
        bsim4::instanceTemp(b4.bsim4v82Instance, *model);
    }
    return bsim4;
}


int main(){
    // Open log file for terminal output
    std::ofstream log_file("loadomp_benchmark_log.txt");
    if (!log_file.is_open()) {
        std::cerr << "Warning: Could not open log file for writing. Continuing without logging.\n";
    }

    // Macro to write to both console and log file
    #define LOG_OUTPUT(stream_expr) \
        do { \
            std::cout << stream_expr; \
            if (log_file.is_open()) { log_file << stream_expr; } \
        } while(0)

    // Print benchmark header
    LOG_OUTPUT("BSIM4 LoadOMP Benchmark - Variable Instance Count\n");
    LOG_OUTPUT("==================================================\n");
    LOG_OUTPUT("CPU cores: " << omp_get_num_procs() << "\n");
    LOG_OUTPUT("OpenMP max threads: " << omp_get_max_threads() << "\n");
    LOG_OUTPUT("Output will be logged to: loadomp_benchmark_log.txt\n\n");

    // Instance counts to test
    std::vector<size_t> instance_counts = {100, 200, 300, 400, 500, 600, 700, 800, 900, 1000};

    // Thread counts to test
    std::vector<int> thread_counts = {1, 2, 4, 8, 16, 24, 32, 48, 64, 96, 128};

    // Matrix type to test: true for sparse (with O(1) indexed stamping), false for dense
    bool use_sparse = true;

    // Create BSIM4 model (shared across all tests)
    auto bsim4model = CreateBSIM4Model();

    // Storage for detailed results (one row per test configuration)
    struct BenchmarkResult {
        size_t num_instances;
        int num_iterations;
        int num_threads;
        double single_thread_time_s;
        double parallel_time_s;
        double parallel_calc_time_s;
        double parallel_apply_time_s;
        double apply_stamps_pct;
        double speedup;
        double efficiency_pct;
    };
    std::vector<BenchmarkResult> all_results;

    // Outer loop: iterate through different instance counts
    for (size_t num_instances : instance_counts) {
        LOG_OUTPUT("\n");
        LOG_OUTPUT("========================================================================\n");
        LOG_OUTPUT("Testing with " << num_instances << " BSIM4 instances\n");
        LOG_OUTPUT("========================================================================\n\n");

        // Calculate iterations (scale down for larger instances to keep runtime reasonable)
        int num_iterations = std::max(1000, static_cast<int>(2e6 / num_instances));
        LOG_OUTPUT("Iterations per test: " << num_iterations << "\n\n");

        // Allocate solution vectors for this instance count
        HybridMatrix LHS(num_instances * 4, num_instances * 4, use_sparse);
        // For dense mode, seed with random values to match original behavior
        if (!use_sparse) {
            LHS.get_dense().randu();
        }
        arma::vec RHS = arma::vec(num_instances * 4, arma::fill::randu);
        arma::vec pre_NR_solution = arma::vec(num_instances * 4, arma::fill::randu);
        std::vector<bsim4::BSIM4stamp> stamps(num_instances);

        // Create circuit with instances
        CKTcircuit ckt;
        ckt.CKTelements.bsim4 = CreateBSIM4Instances(bsim4model, num_instances);

        // Pattern discovery for sparse matrices (one-time setup)
        if (use_sparse) {
            for (auto &b4 : ckt.CKTelements.bsim4) {
                bsim4::bsim4RecordPattern(b4.bsim4v82Instance, LHS);
            }
            LHS.lock_pattern();
            LOG_OUTPUT("Sparse matrix pattern locked: " << LHS.pattern_size() << " non-zeros\n");
        }

        // Warm-up run to initialize caches
        LOG_OUTPUT("Running warm-up...\n");
        loadsingle(ckt, pre_NR_solution, LHS, RHS);
        LOG_OUTPUT("Warm-up complete.\n\n");

        // Benchmark single-threaded execution
        LOG_OUTPUT("Benchmarking single-threaded execution...\n");
        auto start = std::chrono::high_resolution_clock::now();
        for (int iter = 0; iter < num_iterations; ++iter) {
            LHS.zeros();
            RHS.zeros();
            loadsingle(ckt, pre_NR_solution, LHS, RHS);
        }
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> single_time = end - start;

        LOG_OUTPUT("Single-threaded total time: " << std::fixed << std::setprecision(6)
                  << single_time.count() << " seconds\n\n");

        // Test different thread counts
        double best_speedup = 0;
        int best_threads = 1;

        LOG_OUTPUT("Testing different thread counts:\n");
        LOG_OUTPUT("--------------------------------\n");

        for (int num_threads : thread_counts) {
            if (num_threads > omp_get_num_procs()) {
                continue; // Skip silently
            }

            omp_set_num_threads(num_threads);

            double total_calc_time = 0.0;
            double total_apply_time = 0.0;

            start = std::chrono::high_resolution_clock::now();
            for (int iter = 0; iter < num_iterations; ++iter) {
                LHS.zeros();
                RHS.zeros();
                LoadOMPTiming timing = loadomp(ckt, pre_NR_solution, LHS, RHS, stamps);
                total_calc_time += timing.parallel_calc_time;
                total_apply_time += timing.apply_stamps_time;
            }
            end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> parallel_time = end - start;

            double speedup = single_time.count() / parallel_time.count();
            double efficiency = (speedup / num_threads) * 100.0;
            double total_loadomp_time = total_calc_time + total_apply_time;
            double apply_percentage = (total_apply_time / total_loadomp_time) * 100.0;

            LOG_OUTPUT("Parallel (" << std::setw(3) << num_threads << " threads): "
                      << std::fixed << std::setprecision(6) << parallel_time.count()
                      << " s  |  Speedup: " << std::setprecision(2) << std::setw(5) << speedup
                      << "x  |  Eff: " << std::setprecision(1) << std::setw(4) << efficiency << "%"
                      << "  |  applyStamps: " << std::setprecision(1) << std::setw(4) << apply_percentage << "%\n");
            std::cout.flush();
            if (log_file.is_open()) { log_file.flush(); }

            // Store result for this configuration
            all_results.push_back({
                num_instances,                    // NumInstances
                num_iterations,                   // NumIterations
                num_threads,                      // NumThreads
                single_time.count(),              // SingleThreadTime_s
                parallel_time.count(),            // ParallelTime_s
                total_calc_time,                  // ParallelCalcTime_s
                total_apply_time,                 // ParallelApplyTime_s
                apply_percentage,                 // ApplyStamps_pct
                speedup,                          // Speedup
                efficiency                        // Efficiency_pct
            });

            if (speedup > best_speedup) {
                best_speedup = speedup;
                best_threads = num_threads;
            }
        }

        LOG_OUTPUT("\nBest result: " << best_threads << " threads with "
                  << std::fixed << std::setprecision(2) << best_speedup << "x speedup\n");
    }

    // Write results to CSV file
    LOG_OUTPUT("\n\n");
    LOG_OUTPUT("================================================================================\n");
    LOG_OUTPUT("Writing results to CSV file...\n");
    LOG_OUTPUT("================================================================================\n");

    std::ofstream csv_file("loadomp_benchmark_results.csv");
    if (!csv_file.is_open()) {
        std::cerr << "Error: Could not open CSV file for writing!\n";
        return 1;
    }

    // Write CSV header
    csv_file << "NumInstances,NumIterations,NumThreads,SingleThreadTime_s,"
             << "ParallelTime_s,ParallelCalcTime_s,ParallelApplyTime_s,"
             << "ApplyStamps_pct,Speedup,Efficiency_pct\n";

    // Write all results
    for (const auto& r : all_results) {
        csv_file << r.num_instances << ","
                 << r.num_iterations << ","
                 << r.num_threads << ","
                 << std::fixed << std::setprecision(6) << r.single_thread_time_s << ","
                 << r.parallel_time_s << ","
                 << r.parallel_calc_time_s << ","
                 << r.parallel_apply_time_s << ","
                 << std::setprecision(2) << r.apply_stamps_pct << ","
                 << r.speedup << ","
                 << r.efficiency_pct << "\n";
    }

    csv_file.close();
    LOG_OUTPUT("Results written to: loadomp_benchmark_results.csv\n");
    LOG_OUTPUT("Total result rows: " << all_results.size() << "\n");

    // Print brief summary
    LOG_OUTPUT("\n");
    LOG_OUTPUT("================================================================================\n");
    LOG_OUTPUT("SUMMARY - Best Speedups by Instance Count\n");
    LOG_OUTPUT("================================================================================\n");
    LOG_OUTPUT("Instances | Best Threads | Best Speedup | Best Efficiency\n");
    LOG_OUTPUT("----------|--------------|--------------|----------------\n");

    // Group results by instance count and find best speedups
    for (size_t num_inst : instance_counts) {
        double best_speedup = 0;
        int best_threads = 0;
        double best_efficiency = 0;

        for (const auto& r : all_results) {
            if (r.num_instances == num_inst && r.speedup > best_speedup) {
                best_speedup = r.speedup;
                best_threads = r.num_threads;
                best_efficiency = r.efficiency_pct;
            }
        }

        LOG_OUTPUT(std::setw(9) << num_inst << " | "
                  << std::setw(12) << best_threads << " | "
                  << std::fixed << std::setprecision(2) << std::setw(12) << best_speedup << "x | "
                  << std::setprecision(1) << std::setw(15) << best_efficiency << "%\n");
    }
    LOG_OUTPUT("================================================================================\n");

    // Close log file
    if (log_file.is_open()) {
        log_file.close();
        std::cout << "\nLog file written to: loadomp_benchmark_log.txt\n";
    }

    #undef LOG_OUTPUT

    return 0;
}

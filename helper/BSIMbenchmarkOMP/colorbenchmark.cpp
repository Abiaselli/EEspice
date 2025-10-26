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
#include <map>
#include <omp.h>

// EEspice includes
#include "CKT.hpp"
#include "global.hpp"
#include "bsim4v82/bsim4v82setup.hpp"
#include "bsim4v82/bsim4v82temp.hpp"

// Benchmarking includes
#include "color.hpp"
#include "single.hpp"


// Removed global num_instances - now variable in main loop
// Removed global stamps vector - now allocated per instance count


// Generate node distribution count sequence: 2 (for ~num_instances/2 colors) and num_instances (for 1 color)
std::vector<int> generateNodeDistributionSequence(int max_node_count) {
    return {2, max_node_count};
}

std::shared_ptr<bsim4::BSIM4model> CreateBSIM4Model(){
    auto model = std::make_shared<bsim4::BSIM4model>();
    // Setup the model using the default temperature
    bsim4::modelSetup(*model, nomTemp);
    // Setup the model using the actual temperature
    bsim4::modelTemp(*model, nomTemp);
    return model;
}

void CreateNodes(std::vector<BSIM4> &bsim4, int node_distribution_count){
    // Create nodes for the circuit, distributing instances across node_distribution_count distinct node values
    for (int i = 0; i < bsim4.size(); ++i) {
        // Distribute instances evenly across node_distribution_count different node values
        const int node = (i % node_distribution_count) + 1;
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

std::vector<BSIM4> CreateBSIM4Instances(std::shared_ptr<bsim4::BSIM4model> &model, size_t num_instances, int node_distribution_count){
    std::vector<BSIM4> bsim4(num_instances);
    CKTcircuit ckt; // Dummy circuit for instance setup
    CreateNodes(bsim4, node_distribution_count); // Create nodes for each instance
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

    // Print benchmark header
    std::cout << "BSIM4 Graph Coloring Benchmark - Color Scan with Variable Instances\n";
    std::cout << "===================================================================\n";
    std::cout << "CPU cores: " << omp_get_num_procs() << "\n";
    std::cout << "OpenMP max threads: " << omp_get_max_threads() << "\n\n";

    // Instance counts to test
    std::vector<size_t> instance_counts = {100, 200, 300, 400, 500, 600, 700, 800, 900, 1000};

    // Thread counts to test
    std::vector<int> thread_counts = {1, 2, 4, 8, 16, 24, 32, 48, 64, 96, 128};

    // Create BSIM4 model (shared across all tests)
    auto bsim4model = CreateBSIM4Model();

    // Storage for detailed results (one row per test configuration)
    struct BenchmarkResult {
        size_t num_instances;
        int node_distribution_count;
        int actual_colors;
        double coloring_time_s;
        double single_thread_time_s;
        int num_threads;
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
        std::cout << "\n";
        std::cout << "========================================================================\n";
        std::cout << "Testing with " << num_instances << " BSIM4 instances\n";
        std::cout << "========================================================================\n\n";

        // Calculate iterations (scale down for larger instances to keep runtime reasonable)
        int num_iterations = std::max(1000, static_cast<int>(2e6 / num_instances));
        std::cout << "Iterations per test: " << num_iterations << "\n\n";

        // Generate node distribution sequence for this instance count
        std::vector<int> node_distribution_sequence = generateNodeDistributionSequence(num_instances);

        // Allocate solution vectors for this instance count
        arma::mat LHS = arma::mat(num_instances * 4, num_instances * 4, arma::fill::randu);
        arma::vec RHS = arma::vec(num_instances * 4, arma::fill::randu);
        arma::vec pre_NR_solution = arma::vec(num_instances * 4, arma::fill::randu);
        std::vector<bsim4::BSIM4stamp> stamps(num_instances);

        // Warm-up run to initialize caches
        std::cout << "Running warm-up...\n";
        CKTcircuit warmup_ckt;
        warmup_ckt.CKTelements.bsim4 = CreateBSIM4Instances(bsim4model, num_instances, node_distribution_sequence[0]);
        loadsingle(warmup_ckt, pre_NR_solution, LHS, RHS);
        std::cout << "Warm-up complete.\n\n";

        // Middle loop: iterate through each node distribution count
        for (int node_distribution_count : node_distribution_sequence) {
            std::cout << "\n";
            std::cout << "----------------------------------------\n";
            std::cout << "Node Distribution Count: " << node_distribution_count << "\n";
            std::cout << "----------------------------------------\n\n";

            // Create circuit with specified node distribution count
            CKTcircuit ckt;
            ckt.CKTelements.bsim4 = CreateBSIM4Instances(bsim4model, num_instances, node_distribution_count);

            // Compute graph coloring
            std::cout << "Computing graph coloring...\n";
            auto start_color = std::chrono::high_resolution_clock::now();
            BSIM4Coloring coloring;
            coloring.computeColoring(ckt.CKTelements.bsim4);
            auto end_color = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> coloring_time = end_color - start_color;

            // Output coloring results
            size_t actual_colors = coloring.getNumColors();
            std::cout << "Actual colors: " << actual_colors << "\n";
            std::cout << "Coloring time: " << std::fixed << std::setprecision(6)
                      << coloring_time.count() << " seconds\n";

            // Benchmark single-threaded execution
            std::cout << "Benchmarking single-threaded execution...\n";
            auto start = std::chrono::high_resolution_clock::now();
            for (int iter = 0; iter < num_iterations; ++iter) {
                loadsingle(ckt, pre_NR_solution, LHS, RHS);
            }
            auto end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> single_time = end - start;

            std::cout << "Single-threaded total time: " << std::fixed << std::setprecision(6)
                      << single_time.count() << " seconds\n\n";

            // Test different thread counts
            double best_speedup = 0;
            int best_threads = 1;

            std::cout << "Testing different thread counts:\n";
            std::cout << "--------------------------------\n";

            for (int num_threads : thread_counts) {
                if (num_threads > omp_get_num_procs()) {
                    continue; // Skip silently
                }

                omp_set_num_threads(num_threads);

                double total_calc_time = 0.0;
                double total_apply_time = 0.0;

                start = std::chrono::high_resolution_clock::now();
                for (int iter = 0; iter < num_iterations; ++iter) {
                    LoadOMPTiming timing = loadompColor(ckt, pre_NR_solution, LHS, RHS, stamps, coloring);
                    total_calc_time += timing.parallel_calc_time;
                    total_apply_time += timing.apply_stamps_time;
                }
                end = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> parallel_time = end - start;

                double speedup = single_time.count() / parallel_time.count();
                double efficiency = (speedup / num_threads) * 100.0;
                double total_loadomp_time = total_calc_time + total_apply_time;
                double apply_percentage = (total_apply_time / total_loadomp_time) * 100.0;

                std::cout << "Parallel (" << std::setw(3) << num_threads << " threads): "
                          << std::fixed << std::setprecision(6) << parallel_time.count()
                          << " s  |  Speedup: " << std::setprecision(2) << std::setw(5) << speedup
                          << "x  |  Eff: " << std::setprecision(1) << std::setw(4) << efficiency << "%"
                          << "  |  applyStamps: " << std::setprecision(1) << std::setw(4) << apply_percentage << "%\n";
                std::cout.flush();

                // Store result for this configuration
                all_results.push_back({
                    num_instances,                    // NumInstances
                    node_distribution_count,          // NodeDistributionCount
                    static_cast<int>(actual_colors),  // ActualColors
                    coloring_time.count(),            // ColoringTime_s
                    single_time.count(),              // SingleThreadTime_s
                    num_threads,                      // NumThreads
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

            std::cout << "\nBest result: " << best_threads << " threads with "
                      << std::fixed << std::setprecision(2) << best_speedup << "x speedup\n";
        }
    }

    // Write results to CSV file
    std::cout << "\n\n";
    std::cout << "================================================================================\n";
    std::cout << "Writing results to CSV file...\n";
    std::cout << "================================================================================\n";

    std::ofstream csv_file("coloring_benchmark_results.csv");
    if (!csv_file.is_open()) {
        std::cerr << "Error: Could not open CSV file for writing!\n";
        return 1;
    }

    // Write CSV header
    csv_file << "NumInstances,NodeDistributionCount,ActualColors,ColoringTime_s,"
             << "SingleThreadTime_s,NumThreads,ParallelTime_s,ParallelCalcTime_s,"
             << "ParallelApplyTime_s,ApplyStamps_pct,Speedup,Efficiency_pct\n";

    // Write all results
    for (const auto& r : all_results) {
        csv_file << r.num_instances << ","
                 << r.node_distribution_count << ","
                 << r.actual_colors << ","
                 << std::fixed << std::setprecision(6) << r.coloring_time_s << ","
                 << r.single_thread_time_s << ","
                 << r.num_threads << ","
                 << r.parallel_time_s << ","
                 << r.parallel_calc_time_s << ","
                 << r.parallel_apply_time_s << ","
                 << std::setprecision(2) << r.apply_stamps_pct << ","
                 << r.speedup << ","
                 << r.efficiency_pct << "\n";
    }

    csv_file.close();
    std::cout << "Results written to: coloring_benchmark_results.csv\n";
    std::cout << "Total result rows: " << all_results.size() << "\n";

    // Print brief summary
    std::cout << "\n";
    std::cout << "================================================================================\n";
    std::cout << "SUMMARY - Best Speedups by Instance Count\n";
    std::cout << "================================================================================\n";
    std::cout << "Instances | Colors (2-node) | Best Speedup | Colors (N-node) | Best Speedup\n";
    std::cout << "----------|-----------------|--------------|-----------------|-------------\n";

    // Group results by instance count and find best speedups
    for (size_t num_inst : instance_counts) {
        double best_speedup_2node = 0;
        int colors_2node = 0;
        double best_speedup_nnode = 0;
        int colors_nnode = 0;

        for (const auto& r : all_results) {
            if (r.num_instances == num_inst) {
                if (r.node_distribution_count == 2 && r.speedup > best_speedup_2node) {
                    best_speedup_2node = r.speedup;
                    colors_2node = r.actual_colors;
                } else if (r.node_distribution_count == static_cast<int>(num_inst) && r.speedup > best_speedup_nnode) {
                    best_speedup_nnode = r.speedup;
                    colors_nnode = r.actual_colors;
                }
            }
        }

        std::cout << std::setw(9) << num_inst << " | "
                  << std::setw(15) << colors_2node << " | "
                  << std::fixed << std::setprecision(2) << std::setw(12) << best_speedup_2node << "x | "
                  << std::setw(15) << colors_nnode << " | "
                  << std::setw(12) << best_speedup_nnode << "x\n";
    }
    std::cout << "================================================================================\n";

    return 0;
}

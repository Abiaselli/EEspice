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


// Generate node distribution sequence to achieve target color counts from 1 to max_colors
// The relationship is approximately: num_colors ≈ ceil(num_instances / node_distribution_count)
// So for target_colors, we use: node_distribution_count ≈ ceil(num_instances / target_colors)
std::vector<int> generateNodeDistributionSequence(int num_instances, int max_target_colors, int step_size) {
    std::vector<int> sequence;
    
    // Add configurations for different target color counts
    for (int target_colors = 1; target_colors <= max_target_colors; target_colors += step_size) {
        // Calculate node distribution count to approximately achieve target_colors
        // Use ceiling division: node_dist_count = ceil(num_instances / target_colors)
        // Enforce minimum of 2 to cap at ~N/2 colors (spec requirement)
        int node_dist_count = std::max(2, (num_instances + target_colors - 1) / target_colors);
        
        // Ensure we don't exceed num_instances
        node_dist_count = std::min(node_dist_count, num_instances);
        
        // Avoid duplicates in the sequence
        if (sequence.empty() || sequence.back() != node_dist_count) {
            sequence.push_back(node_dist_count);
        }
    }
    
    // Always ensure we test the extreme of maximum node distribution (1 color expected)
    if (std::find(sequence.begin(), sequence.end(), num_instances) == sequence.end()) {
        sequence.push_back(num_instances);
    }
    
    // Note: We do NOT add node_dist=1 since that would give N colors, exceeding N/2
    // The minimum node_dist is 2, giving approximately N/2 colors
    
    // Sort in descending order (from fewer colors to more colors)
    std::sort(sequence.begin(), sequence.end(), std::greater<int>());
    
    return sequence;
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
    // Open log file for terminal output
    std::ofstream log_file("coloring_benchmark_log.txt");
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
    LOG_OUTPUT("BSIM4 Graph Coloring Benchmark - Color Scan from 1 to N/2 colors\n");
    LOG_OUTPUT("==================================================================\n");
    LOG_OUTPUT("CPU cores: " << omp_get_num_procs() << "\n");
    LOG_OUTPUT("OpenMP max threads: " << omp_get_max_threads() << "\n");
    LOG_OUTPUT("Output will be logged to: coloring_benchmark_log.txt\n\n");

    // Instance counts to test
    std::vector<size_t> instance_counts = {100, 200, 300, 400, 500, 600, 700, 800, 900, 1000};

    // Thread counts to test
    std::vector<int> thread_counts = {1, 2, 4, 8, 16, 24, 32, 48, 64, 96, 128};

    // Create BSIM4 model (shared across all tests)
    auto bsim4model = CreateBSIM4Model();

    // Storage for detailed results (one row per test configuration)
    struct BenchmarkResult {
        size_t num_instances;
        int num_iterations;
        int node_distribution_count;
        int target_colors;
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
        LOG_OUTPUT("\n");
        LOG_OUTPUT("========================================================================\n");
        LOG_OUTPUT("Testing with " << num_instances << " BSIM4 instances\n");
        LOG_OUTPUT("========================================================================\n\n");

        // Calculate iterations (scale down for larger instances to keep runtime reasonable)
        int num_iterations = std::max(1000, static_cast<int>(2e6 / num_instances));
        LOG_OUTPUT("Iterations per test: " << num_iterations << "\n");

        // Calculate step size for color iteration
        int max_target_colors = num_instances / 2;
        int step_size = std::max(1, max_target_colors / 20); // Aim for ~20 different color configurations
        
        LOG_OUTPUT("Target colors: 1 to " << max_target_colors << " (step size: " << step_size << ")\n\n");

        // Generate node distribution sequence for different target color counts
        std::vector<int> node_distribution_sequence = generateNodeDistributionSequence(num_instances, max_target_colors, step_size);
        
        LOG_OUTPUT("Node distribution counts to test: ");
        for (int ndc : node_distribution_sequence) {
            LOG_OUTPUT(ndc << " ");
        }
        LOG_OUTPUT("\n");
        LOG_OUTPUT("Expected approximate colors: ");
        for (int ndc : node_distribution_sequence) {
            // Use ceiling division for more accurate expectation
            int expected_colors = (num_instances + ndc - 1) / ndc;
            LOG_OUTPUT(expected_colors << " ");
        }
        LOG_OUTPUT("\n\n");

        // Allocate solution vectors for this instance count
        arma::mat LHS = arma::mat(num_instances * 4, num_instances * 4, arma::fill::randu);
        arma::vec RHS = arma::vec(num_instances * 4, arma::fill::randu);
        arma::vec pre_NR_solution = arma::vec(num_instances * 4, arma::fill::randu);
        std::vector<bsim4::BSIM4stamp> stamps(num_instances);

        // Warm-up run to initialize caches
        LOG_OUTPUT("Running warm-up...\n");
        CKTcircuit warmup_ckt;
        warmup_ckt.CKTelements.bsim4 = CreateBSIM4Instances(bsim4model, num_instances, node_distribution_sequence[0]);
        loadsingle(warmup_ckt, pre_NR_solution, LHS, RHS);
        LOG_OUTPUT("Warm-up complete.\n\n");

        // Middle loop: iterate through each node distribution count
        for (int node_distribution_count : node_distribution_sequence) {
            // Use ceiling division for consistent calculation
            int target_colors = (num_instances + node_distribution_count - 1) / node_distribution_count;
            
            LOG_OUTPUT("\n");
            LOG_OUTPUT("----------------------------------------\n");
            LOG_OUTPUT("Node Distribution Count: " << node_distribution_count << "\n");
            LOG_OUTPUT("Target Colors (approx): " << target_colors << "\n");
            LOG_OUTPUT("----------------------------------------\n\n");

            // Create circuit with specified node distribution count
            CKTcircuit ckt;
            ckt.CKTelements.bsim4 = CreateBSIM4Instances(bsim4model, num_instances, node_distribution_count);

            // Compute graph coloring
            LOG_OUTPUT("Computing graph coloring...\n");
            auto start_color = std::chrono::high_resolution_clock::now();
            BSIM4Coloring coloring;
            coloring.computeColoring(ckt.CKTelements.bsim4);
            auto end_color = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> coloring_time = end_color - start_color;

            // Output coloring results
            size_t actual_colors = coloring.getNumColors();
            LOG_OUTPUT("Actual colors: " << actual_colors << " (target was ~" << target_colors << ")\n");
            LOG_OUTPUT("Coloring time: " << std::fixed << std::setprecision(6)
                      << coloring_time.count() << " seconds\n");

            // Benchmark single-threaded execution
            LOG_OUTPUT("Benchmarking single-threaded execution...\n");
            auto start = std::chrono::high_resolution_clock::now();
            for (int iter = 0; iter < num_iterations; ++iter) {
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
                    node_distribution_count,          // NodeDistributionCount
                    target_colors,                    // TargetColors
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

            LOG_OUTPUT("\nBest result: " << best_threads << " threads with "
                      << std::fixed << std::setprecision(2) << best_speedup << "x speedup\n");
        }
    }

    // Write results to CSV file
    LOG_OUTPUT("\n\n");
    LOG_OUTPUT("================================================================================\n");
    LOG_OUTPUT("Writing results to CSV file...\n");
    LOG_OUTPUT("================================================================================\n");

    std::ofstream csv_file("coloring_benchmark_results.csv");
    if (!csv_file.is_open()) {
        std::cerr << "Error: Could not open CSV file for writing!\n";
        return 1;
    }

    // Write CSV header
    csv_file << "NumInstances,NumIterations,NodeDistributionCount,TargetColors,ActualColors,ColoringTime_s,"
             << "SingleThreadTime_s,NumThreads,ParallelTime_s,ParallelCalcTime_s,"
             << "ParallelApplyTime_s,ApplyStamps_pct,Speedup,Efficiency_pct\n";

    // Write all results
    for (const auto& r : all_results) {
        csv_file << r.num_instances << ","
                 << r.num_iterations << ","
                 << r.node_distribution_count << ","
                 << r.target_colors << ","
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
    LOG_OUTPUT("Results written to: coloring_benchmark_results.csv\n");
    LOG_OUTPUT("Total result rows: " << all_results.size() << "\n");

    // Print brief summary grouped by color count ranges
    LOG_OUTPUT("\n");
    LOG_OUTPUT("================================================================================\n");
    LOG_OUTPUT("SUMMARY - Best Speedups by Instance Count and Color Range\n");
    LOG_OUTPUT("================================================================================\n");
    LOG_OUTPUT("Instances | Color Range | Best Speedups (1-10 colors | 11-50 | 51-100 | 100+)\n");
    LOG_OUTPUT("----------|-------------|----------------------------------------------------\n");

    // Group results by instance count and color ranges
    for (size_t num_inst : instance_counts) {
        std::map<std::string, double> best_speedups;
        std::map<std::string, int> color_counts;
        
        // Define color ranges
        std::vector<std::pair<int, int>> ranges = {{1, 10}, {11, 50}, {51, 100}, {101, INT_MAX}};
        
        for (const auto& r : all_results) {
            if (r.num_instances == num_inst) {
                for (const auto& range : ranges) {
                    if (r.actual_colors >= range.first && r.actual_colors <= range.second) {
                        std::string key = std::to_string(range.first) + "-" + 
                                        (range.second == INT_MAX ? "max" : std::to_string(range.second));
                        if (best_speedups.find(key) == best_speedups.end() || r.speedup > best_speedups[key]) {
                            best_speedups[key] = r.speedup;
                            color_counts[key] = r.actual_colors;
                        }
                        break;
                    }
                }
            }
        }

        LOG_OUTPUT(std::setw(9) << num_inst << " | ");
        LOG_OUTPUT(std::setw(11) << "1-" + std::to_string(num_inst/2) << " | ");
        
        for (const auto& range : ranges) {
            std::string key = std::to_string(range.first) + "-" + 
                            (range.second == INT_MAX ? "max" : std::to_string(range.second));
            if (best_speedups.find(key) != best_speedups.end()) {
                LOG_OUTPUT(std::fixed << std::setprecision(2) << std::setw(4) 
                          << best_speedups[key] << "x(" << color_counts[key] << "c) | ");
            } else {
                LOG_OUTPUT("    N/A     | ");
            }
        }
        LOG_OUTPUT("\n");
    }
    LOG_OUTPUT("================================================================================\n");
    LOG_OUTPUT("Note: Numbers in parentheses show actual color count for that best speedup\n");

    // Close log file
    if (log_file.is_open()) {
        log_file.close();
        std::cout << "\nLog file written to: coloring_benchmark_log.txt\n";
    }

    #undef LOG_OUTPUT

    return 0;
}

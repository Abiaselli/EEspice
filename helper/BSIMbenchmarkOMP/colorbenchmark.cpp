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


constexpr size_t num_instances = 60 * 5; // Number of BSIM4 instances
constexpr int num_iterations = 1e7 / 5; // Number of iterations for benchmarking
std::vector<bsim4::BSIM4stamp> stamps(num_instances);


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
    std::cout << "BSIM4 Graph Coloring Benchmark - Color Scan\n";
    std::cout << "============================================\n";
    std::cout << "Number of BSIM4 instances: " << num_instances << "\n";
    std::cout << "CPU cores: " << omp_get_num_procs() << "\n";
    std::cout << "OpenMP max threads: " << omp_get_max_threads() << "\n";
    std::cout << "Iterations per test: " << num_iterations << "\n\n";

    // 1. Generate node distribution sequence (2 test cases: 2 nodes and num_instances nodes)
    std::vector<int> node_distribution_sequence = generateNodeDistributionSequence(num_instances);
    std::cout << "Testing " << node_distribution_sequence.size() << " color configurations\n\n";

    // 2. Create BSIM4 model (shared across all tests)
    auto bsim4model = CreateBSIM4Model();

    // 3. Pre-allocate solution vectors (reused across all tests)
    arma::mat LHS = arma::mat(num_instances * 4, num_instances * 4, arma::fill::randu);
    arma::vec RHS = arma::vec(num_instances * 4, arma::fill::randu);
    arma::vec pre_NR_solution = arma::vec(num_instances * 4, arma::fill::randu);

    // 4. Warm-up run to initialize caches (using first node distribution count)
    std::cout << "Running warm-up (with " << node_distribution_sequence[0] << " distinct nodes)...\n";
    CKTcircuit warmup_ckt;
    warmup_ckt.CKTelements.bsim4 = CreateBSIM4Instances(bsim4model, num_instances, node_distribution_sequence[0]);
    loadsingle(warmup_ckt, pre_NR_solution, LHS, RHS);
    std::cout << "Warm-up complete.\n\n";

    // 5. Thread counts to test
    std::vector<int> thread_counts = {1, 2, 4, 8, 16, 24, 32, 48, 64, 96, 128};

    // 6. Storage for summary results
    struct ColorResult {
        int node_distribution_count;
        int actual_colors;
        double coloring_time;
        double single_time;
        int best_threads;
        double best_speedup;
    };
    std::vector<ColorResult> results;

    // 7. Main loop: iterate through each node distribution count
    for (int node_distribution_count : node_distribution_sequence) {
        std::cout << "\n";
        std::cout << "========================================\n";
        std::cout << "Color Configuration Test\n";
        std::cout << "========================================\n\n";

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

        // Print color group distribution
        const auto& color_groups = coloring.getColorGroups();
        // std::cout << "Color group sizes: ";
        // for (size_t i = 0; i < actual_colors; ++i) {
        //     std::cout << color_groups[i].size();
        //     if (i < actual_colors - 1) std::cout << ", ";
        // }
        // std::cout << "\n\n";

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

            if (speedup > best_speedup) {
                best_speedup = speedup;
                best_threads = num_threads;
            }
        }

        std::cout << "\nBest for " << actual_colors << " colors: " << best_threads << " threads with "
                  << std::fixed << std::setprecision(2) << best_speedup << "x speedup\n";

        // Store results for summary
        results.push_back({node_distribution_count, static_cast<int>(actual_colors), coloring_time.count(),
                          single_time.count(), best_threads, best_speedup});
    }

    // Print summary table
    std::cout << "\n\n";
    std::cout << "================================================================================\n";
    std::cout << "SUMMARY - Performance vs Colors\n";
    std::cout << "================================================================================\n";
    std::cout << "Actual Colors | Coloring Time (s) | Best Speedup | Best Threads\n";
    std::cout << "--------------|-------------------|--------------|-------------\n";

    for (const auto& r : results) {
        std::cout << std::setw(13) << r.actual_colors << " | "
                  << std::fixed << std::setprecision(6) << std::setw(17) << r.coloring_time << " | "
                  << std::setprecision(2) << std::setw(12) << r.best_speedup << "x | "
                  << std::setw(12) << r.best_threads << "\n";
    }
    std::cout << "================================================================================\n";

    return 0;
}

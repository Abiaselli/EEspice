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
#include "bsim4v82/bsim4v82setup.hpp"
#include "bsim4v82/bsim4v82temp.hpp"

// Benchmarking includes
#include "calomp.hpp"
#include "loadomp.hpp"
#include "single.hpp"


constexpr size_t num_instances = 60 * 5; // Number of BSIM4 instances
constexpr int num_iterations = 1e7 / 5; // Number of iterations for benchmarking
std::vector<bsim4::BSIM4stamp> stamps(num_instances);


std::shared_ptr<bsim4::BSIM4model> CreateBSIM4Model(){
    auto model = std::make_shared<bsim4::BSIM4model>();
    // Setup the model using the default temperature
    bsim4::modelSetup(*model, nomTemp);
    // Setup the model using the actual temperature
    bsim4::modelTemp(*model, nomTemp);
    return model;
}

void CreateNodes(std::vector<BSIM4> &bsim4){
    // Create nodes for the circuit
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

    // 1. Create a circuit
    CKTcircuit ckt;

    // 2. bsim4 model
    auto bsim4model = CreateBSIM4Model();

    // 3. Create some instances
    ckt.CKTelements.bsim4 = CreateBSIM4Instances(bsim4model, num_instances);

    // 4. pre solution vector
    arma::mat LHS = arma::mat(num_instances * 4, num_instances * 4, arma::fill::randu); // Assuming 4 nodes per instance
    arma::vec RHS = arma::vec(num_instances * 4, arma::fill::randu);
    arma::vec pre_NR_solution = arma::vec(num_instances * 4, arma::fill::randu); // Assuming 4 nodes per instance

    // Print benchmark header
    std::cout << "BSIM4 Calculation Benchmark\n";
    std::cout << "============================\n";
    std::cout << "Number of BSIM4 instances: " << num_instances << "\n";
    std::cout << "CPU cores: " << omp_get_num_procs() << "\n";
    std::cout << "OpenMP max threads: " << omp_get_max_threads() << "\n";
    // std::cout << "Thread pool threads: " << pool.get_thread_count() << "\n\n";

    // Warm-up run to initialize caches
    std::cout << "Running warm-up...\n";
    calsingle(ckt, pre_NR_solution, LHS, RHS, stamps);
    std::cout << "Warm-up complete.\n\n";

    // Benchmark single-threaded execution (simulating num_iterations of N-R iterations)
    std::cout << "Benchmarking single-threaded execution (" << num_iterations << " iterations)...\n";
    auto start = std::chrono::high_resolution_clock::now();
    for (int iter = 0; iter < num_iterations; ++iter) {
        calsingle(ckt, pre_NR_solution, LHS, RHS, stamps);
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> single_time = end - start;

    std::cout << "Single-threaded total time: " << std::fixed << std::setprecision(6)
              << single_time.count() << " seconds\n\n";

    // Test different thread counts
    std::vector<int> thread_counts = {1, 2, 4, 8, 16, 24, 32, 48, 64, 96, 128};
    double best_speedup = 0;
    int best_threads = 1;

    std::cout << "Testing different thread counts:\n";
    std::cout << "================================\n";

    for (int num_threads : thread_counts) {
        if (num_threads > omp_get_num_procs()) {
            std::cout << "Skipping " << num_threads << " threads (exceeds available cores)\n";
            continue;
        }

        omp_set_num_threads(num_threads);

        start = std::chrono::high_resolution_clock::now();
        for (int iter = 0; iter < num_iterations; ++iter) {
            calomp(ckt, pre_NR_solution, LHS, RHS, stamps);
        }
        end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> parallel_time = end - start;

        double speedup = single_time.count() / parallel_time.count();
        double efficiency = (speedup / num_threads) * 100.0;

        std::cout << "Parallel (" << std::setw(2) << num_threads << " threads): "
                  << std::fixed << std::setprecision(6) << parallel_time.count()
                  << " seconds  |  Speedup: " << std::setprecision(2) << speedup
                  << "x  |  Efficiency: " << std::setprecision(1) << efficiency << "%\n";
        std::cout.flush();

        if (speedup > best_speedup) {
            best_speedup = speedup;
            best_threads = num_threads;
        }
    }

    std::cout << "\nBest configuration: " << best_threads << " threads with "
              << std::fixed << std::setprecision(2) << best_speedup << "x speedup\n";

    return 0;
}
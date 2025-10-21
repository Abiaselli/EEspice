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
#include "color.hpp"


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

    // 5. Create color groups
    auto start_color = std::chrono::high_resolution_clock::now();
    BSIM4Coloring coloring;
    coloring.computeColoring(ckt.CKTelements.bsim4);
    auto end_color = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> coloring_time = end_color - start_color;

    // Output coloring results
    int num_colors = coloring.getNumColors();
    std::cout << "Number of colors: " << num_colors << "\n";
    std::cout << "Coloring time: " << std::fixed << std::setprecision(6)
              << coloring_time.count() << "s\n\n";

    // Test different thread counts
    std::vector<int> thread_counts = {1, 2, 4, 8, 16, 20, 24};
    double baseline_time = 0.0;
    auto start = std::chrono::high_resolution_clock::now();
    auto end = std::chrono::high_resolution_clock::now();

    std::cout << "Threads  Speedup  Efficiency  LoadOMP Time  Apply %\n";

    for (int num_threads : thread_counts) {
        if (num_threads > omp_get_num_procs()) {
            continue;
        }

        omp_set_num_threads(num_threads);

        double total_calc_time = 0.0;
        double total_apply_time = 0.0;

        start = std::chrono::high_resolution_clock::now();
        for (int iter = 0; iter < num_iterations; ++iter) {
            LoadOMPColorTiming timing = loadompColor(ckt, pre_NR_solution, LHS, RHS, stamps, coloring);
            total_calc_time += timing.parallel_calc_time;
            total_apply_time += timing.apply_stamps_time;
        }
        end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> parallel_time = end - start;

        // Use 1-thread as baseline
        if (num_threads == 1) {
            baseline_time = parallel_time.count();
        }

        double speedup = baseline_time / parallel_time.count();
        double efficiency = (speedup / num_threads) * 100.0;
        double total_loadomp_time = total_calc_time + total_apply_time;
        double apply_percentage = (total_apply_time / total_loadomp_time) * 100.0;

        std::cout << std::setw(7) << num_threads << "  "
                  << std::fixed << std::setprecision(2) << std::setw(7) << speedup << "x  "
                  << std::setprecision(1) << std::setw(10) << efficiency << "%  "
                  << std::setprecision(6) << std::setw(12) << total_loadomp_time << "s  "
                  << std::setprecision(1) << std::setw(7) << apply_percentage << "%\n";
        std::cout.flush();
    }

    return 0;
}
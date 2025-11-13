#define ARMA_DONT_USE_WRAPPER
#define ARMA_USE_MKL_ALLOC

#include <omp.h>
#include "model_setup.hpp"
#include "Transient.hpp"
#include "DC.hpp"
#include "AC.hpp"
#include "saveCSV.hpp"
#include "batch.hpp"
#include "simulation_exceptions.hpp"

constexpr int num_iterations = 2000; // Number of iterations for benchmarking

struct LoadOMPTiming {
    double parallel_calc_time;
    double apply_stamps_time;
};

void loadsingle(CKTcircuit &ckt, const arma::vec &pre_NR_solution, arma::mat &LHS, arma::vec &RHS){
    if (!ckt.CKTelements.bsim4.empty()) {
        const size_t num_devices = ckt.CKTelements.bsim4.size();

        // Step 1: Serial computation - each iteration processes one device
        for (size_t i = 0; i < ckt.CKTelements.bsim4.size(); ++i) {
            const bsim4::BSIM4model &b4model = *ckt.CKTelements.bsim4[i].bsim4v82Instance.BSIM4modPtr;
            // Create a local copy of the instance state
            bsim4::BSIM4V82 &instance = ckt.CKTelements.bsim4[i].bsim4v82Instance;

            // Calculate stamps
            const auto stamp = bsim4::BSIM4calculateStamps(
                ckt, b4model, instance, ckt.spiceCompatible,
                pre_NR_solution, ckt.CKTtemp, ckt.CKTgmin);

            // Apply the calculated stamps to the global LHS and RHS matrices
            bsim4::bsim4applyStamps(instance, stamp, LHS, RHS);
        }
    }
}

LoadOMPTiming loadompColor(CKTcircuit &ckt, const arma::vec &pre_NR_solution, 
                      arma::mat &LHS, arma::vec &RHS, 
                      std::vector<bsim4::BSIM4stamp> &stamps,
                      const BSIM4Coloring &coloring)
{
    LoadOMPTiming timing{0.0, 0.0};
    
    if (!ckt.CKTelements.bsim4.empty()) {
        const size_t n = ckt.CKTelements.bsim4.size();
        
        // Phase 1: Parallel computation of stamps
        auto start_calc = std::chrono::high_resolution_clock::now();
        #pragma omp parallel for
        for (size_t i = 0; i < n; ++i) {
            const bsim4::BSIM4model &b4model = *ckt.CKTelements.bsim4[i].bsim4v82Instance.BSIM4modPtr;
            bsim4::BSIM4V82 &instance = ckt.CKTelements.bsim4[i].bsim4v82Instance;
            stamps[i] = bsim4::BSIM4calculateStamps(ckt, b4model, instance, 
                                                     ckt.spiceCompatible, pre_NR_solution, 
                                                     ckt.CKTtemp, ckt.CKTgmin);
        }
        auto end_calc = std::chrono::high_resolution_clock::now();
        timing.parallel_calc_time = std::chrono::duration<double>(end_calc - start_calc).count();
        
        // Graph coloring
        const auto& color_groups = coloring.getColorGroups();

        // Phase 2: Parallel application of stamps using coloring
        auto start_apply = std::chrono::high_resolution_clock::now();
        for (const auto& group : color_groups) {
            // All instances in this group can be processed in parallel
            #pragma omp parallel for
            for (size_t idx = 0; idx < group.size(); ++idx) {
                size_t i = group[idx];
                bsim4::BSIM4V82 &instance = ckt.CKTelements.bsim4[i].bsim4v82Instance;
                bsim4::bsim4applyStamps(instance, stamps[i], LHS, RHS);
            }
        }
        auto end_apply = std::chrono::high_resolution_clock::now();
        timing.apply_stamps_time = std::chrono::duration<double>(end_apply - start_apply).count();
    }
    
    return timing;
}

LoadOMPTiming loadompColor2(CKTcircuit &ckt, const arma::vec &pre_NR_solution, 
                      arma::mat &LHS, arma::vec &RHS, 
                      std::vector<bsim4::BSIM4stamp> &stamps,
                      const BSIM4Coloring &coloring)
{
    LoadOMPTiming timing{0.0, 0.0};
    
    if (!ckt.CKTelements.bsim4.empty()) {
        const size_t n = ckt.CKTelements.bsim4.size();
        
        // Phase 1: Parallel computation of stamps
        auto start_calc = std::chrono::high_resolution_clock::now();
        #pragma omp parallel for
        for (size_t i = 0; i < n; ++i) {
            const bsim4::BSIM4model &b4model = *ckt.CKTelements.bsim4[i].bsim4v82Instance.BSIM4modPtr;
            bsim4::BSIM4V82 &instance = ckt.CKTelements.bsim4[i].bsim4v82Instance;
            stamps[i] = bsim4::BSIM4calculateStamps(ckt, b4model, instance, 
                                                     ckt.spiceCompatible, pre_NR_solution, 
                                                     ckt.CKTtemp, ckt.CKTgmin);
        }
        auto end_calc = std::chrono::high_resolution_clock::now();
        timing.parallel_calc_time = std::chrono::duration<double>(end_calc - start_calc).count();
        
        // Graph coloring
        const auto& color_groups = coloring.getColorGroups();

        // Phase 2: Parallel application of stamps using coloring
        auto start_apply = std::chrono::high_resolution_clock::now();
        #pragma omp parallel // <-- Create threads ONCE
        {
            for (const auto& group : color_groups) {
                // All instances in this group can be processed in parallel
                #pragma omp for // <-- Distribute work, threads are already active
                for (size_t idx = 0; idx < group.size(); ++idx) {
                    size_t i = group[idx];
                    bsim4::BSIM4V82 &instance = ckt.CKTelements.bsim4[i].bsim4v82Instance;
                    bsim4::bsim4applyStamps(instance, stamps[i], LHS, RHS);
                }
            }
        }
        auto end_apply = std::chrono::high_resolution_clock::now();
        timing.apply_stamps_time = std::chrono::duration<double>(end_apply - start_apply).count();
    }
    
    return timing;
}

LoadOMPTiming loadompColor3(CKTcircuit &ckt, const arma::vec &pre_NR_solution, 
                      arma::mat &LHS, arma::vec &RHS, 
                      std::vector<bsim4::BSIM4stamp> &stamps,
                      const BSIM4Coloring &coloring)
{
    LoadOMPTiming timing{0.0, 0.0};
    
    if (!ckt.CKTelements.bsim4.empty()) {
        const size_t n = ckt.CKTelements.bsim4.size();
        
        // Phase 1: Parallel computation of stamps
        auto start_calc = std::chrono::high_resolution_clock::now();
        #pragma omp parallel for schedule(dynamic)
        for (size_t i = 0; i < n; ++i) {
            const bsim4::BSIM4model &b4model = *ckt.CKTelements.bsim4[i].bsim4v82Instance.BSIM4modPtr;
            bsim4::BSIM4V82 &instance = ckt.CKTelements.bsim4[i].bsim4v82Instance;
            stamps[i] = bsim4::BSIM4calculateStamps(ckt, b4model, instance, 
                                                     ckt.spiceCompatible, pre_NR_solution, 
                                                     ckt.CKTtemp, ckt.CKTgmin);
        }
        auto end_calc = std::chrono::high_resolution_clock::now();
        timing.parallel_calc_time = std::chrono::duration<double>(end_calc - start_calc).count();
        
        // Graph coloring
        const auto& color_groups = coloring.getColorGroups();

        // Phase 2: Parallel application of stamps using coloring
        auto start_apply = std::chrono::high_resolution_clock::now();
        #pragma omp parallel // <-- Create threads ONCE
        {
            for (const auto& group : color_groups) {
                // All instances in this group can be processed in parallel
                #pragma omp for schedule(dynamic) // <-- Distribute work, threads are already active
                for (size_t idx = 0; idx < group.size(); ++idx) {
                    size_t i = group[idx];
                    bsim4::BSIM4V82 &instance = ckt.CKTelements.bsim4[i].bsim4v82Instance;
                    bsim4::bsim4applyStamps(instance, stamps[i], LHS, RHS);
                }
            }
        }
        auto end_apply = std::chrono::high_resolution_clock::now();
        timing.apply_stamps_time = std::chrono::duration<double>(end_apply - start_apply).count();
    }
    
    return timing;
}

LoadOMPTiming loadompColor4(CKTcircuit &ckt, const arma::vec &pre_NR_solution, 
                      arma::mat &LHS, arma::vec &RHS, 
                      const BSIM4Coloring &coloring)
{
    LoadOMPTiming timing{0.0, 0.0};
    // Graph coloring
    const auto& color_groups = coloring.getColorGroups();
    for (const auto& group : color_groups){
        #pragma omp parallel for
        for (size_t idx = 0; idx < group.size(); ++idx) {
            size_t i = group[idx];
            const bsim4::BSIM4model &b4model = *ckt.CKTelements.bsim4[i].bsim4v82Instance.BSIM4modPtr;
            bsim4::BSIM4V82 &instance = ckt.CKTelements.bsim4[i].bsim4v82Instance;
            bsim4::BSIM4stamp stamp = bsim4::BSIM4calculateStamps(ckt, b4model, instance, 
                                                     ckt.spiceCompatible, pre_NR_solution, 
                                                     ckt.CKTtemp, ckt.CKTgmin);
            bsim4::bsim4applyStamps(instance, stamp, LHS, RHS);
        }
    }
    return timing;
}

int main(int argc, const char **argv){
    // Check if a netlist file is provided as command line argument
    if (argc < 2) {
        std::cerr << "Usage: ./eespice <netlist_file>" << std::endl;
        return 1;
    }
    // Parse netlist file
    Modelmap modmap;
    Circuitmap cktmap;
    CircuitParser parser(argv[1]);
    parser_netlist(parser, cktmap, modmap);

    // Model setup using the temperature
    modelSetup(modmap, nomTemp);

    // CKT circuit setup
    CKTcircuit ckt;
    ckt.map = cktmap; // Assign the circuit map to the CKTcircuit
    auto denseMatrixPtr = std::make_shared<DenseMatrix>();  // Create the DenseMatrix as a shared pointer.
    CKTsetup(ckt, parser, denseMatrixPtr, modmap); // Pass the parser to the ckt and the initialise LHS and RHS matrices
    CKTload(ckt);
    ckt.cktdematrix->set_initmatrix(); // Set the initial LHS and RHS matrices
    ckt.sim_stats.MNA_Matrix_size = ckt.cktdematrix->LHS.n_rows;    // Store the size of MNA matrix to the simulation statistics

    // Check for multithreading from parser
    // if (ckt.CKTmultithreaded){
    //     std::cout << "Multithreading enabled for BSIM4 transistors with " << ckt.num_threads << " threads." << std::endl;
    // }

    // Output coloring results
    std::cout << "Number of colors: " << ckt.sim_stats.num_colors << std::endl;

    auto LHS = ckt.cktdematrix->LHS;
    auto RHS = ckt.cktdematrix->RHS;
    arma::vec pre_NR_solution = arma::vec(RHS.n_rows, arma::fill::randu);

    // Test different thread counts
    std::vector<int> thread_counts = {1, 2, 4, 8, 16, 24, 32, 48, 64, 96, 128};
    double baseline_time = 0.0;
    auto start = std::chrono::high_resolution_clock::now();
    auto end = std::chrono::high_resolution_clock::now();

    std::cout << "Threads  Speedup  Efficiency  LoadOMP Time  Apply %\n";

    // use single-threaded as baseline - run same number of iterations
    start = std::chrono::high_resolution_clock::now();
    for (int iter = 0; iter < num_iterations; ++iter) {
        loadsingle(ckt, pre_NR_solution, LHS, RHS);
    }
    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> loadomp_time = end - start;
    baseline_time = loadomp_time.count();

    for (int num_threads : thread_counts) {
        if (num_threads > omp_get_num_procs()) {
            continue;
        }

        omp_set_num_threads(num_threads);

        double total_calc_time = 0.0;
        double total_apply_time = 0.0;

        start = std::chrono::high_resolution_clock::now();
        for (int iter = 0; iter < num_iterations; ++iter) {
            LoadOMPTiming timing = loadompColor4(ckt, pre_NR_solution, LHS, RHS, ckt.b4coloring);
            total_calc_time += timing.parallel_calc_time;
            total_apply_time += timing.apply_stamps_time;
        }
        end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> parallel_time = end - start;

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
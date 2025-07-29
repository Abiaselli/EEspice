#pragma once
#include "batch_calcs.hpp"
#include "CKT.hpp"
#include "global.hpp"
#include "Transient.hpp"
#include "DC.hpp"
#include "AC.hpp"
#include "simulation_exceptions.hpp"
#include "BS_thread_pool/BS_thread_pool.hpp"

BS::synced_stream sync_out;

namespace batch {

void dummy_task()
{
    // std::cout << "dummy task" << std::endl;
    std::this_thread::sleep_for(std::chrono::microseconds(100));
}

// A map to represent a single circuit configuration, e.g., {"V1.value": 1.0, "M1.W": 10e-6}
using CircuitConfig = std::map<std::string, double>;

// Holds the results for a single simulation run from the batch
struct BatchRunResult {
    CircuitConfig config;
    CKTcircuit ckt; // Store the circuit state after simulation
    std::variant<std::vector<dc::DCResult>, std::vector<Transient>, std::vector<AC::AC>, arma::vec> results;
    std::string simulation_type;
    
    // Error tracking fields
    bool success = false;
    std::string error_type = "";
    std::string error_message = "";
};

// This is the worker function that will be executed by each thread in the pool.
// It takes a copy of the circuit template, applies a specific configuration, and runs the simulation.
BatchRunResult simulation_worker(
    const Circuitmap &cktmap,
    CircuitParser parser, // Pass by value to avoid issues with thread safety
    const Modelmap& modmap,
    CircuitConfig config) 
{   
    BatchRunResult run_result;
    run_result.config = config;
    run_result.success = false; // Default to failed
    
    try {
        // ======================== FIX STARTS HERE ========================
        // Create a deep copy of the Modelmap for thread-local use.
        Modelmap local_modmap;

        // Deep copy BSIM4 models by creating new model objects.
        for (const auto& [name, model_ptr] : modmap.bsim4Models) {
            local_modmap.bsim4Models[name] = std::make_shared<bsim4::BSIM4model>(*model_ptr);
        }

        // Shallow copy is sufficient for other models as they are not modified in-place.
        local_modmap.nmosModels = modmap.nmosModels;
        local_modmap.pmosModels = modmap.pmosModels;
        // ========================= FIX ENDS HERE =========================

        // 1. Apply the specific parameter configuration for this simulation run
        apply_circuit_config(parser.elements, config);

        // 2. Setup the matrix with the new configuration
        CKTcircuit ckt_template;
        ckt_template.map = cktmap; // Assign the circuit map to the CKTcircuit
        auto denseMatrixPtr = std::make_shared<DenseMatrix>();
        CKTsetup(ckt_template, parser, denseMatrixPtr, local_modmap);
        CKTload(ckt_template);
        ckt_template.cktdematrix->set_initmatrix();


        // 3. Run the appropriate analysis based on the netlist commands
        if (parser.is_op) {
            run_result.simulation_type = "op";
            bool non_linear = false;
            if(!ckt_template.CKTelements.nmos.empty() || !ckt_template.CKTelements.pmos.empty() || !ckt_template.CKTelements.diodes.empty()){
                non_linear = true;
            }
            run_result.results = OperatingPointAnalysis(ckt_template, local_modmap, non_linear);
        }
        else if (parser.is_transient) {
            run_result.simulation_type = "tran";
            TransientSimulator trans_sim = Transsetup(parser, ckt_template);
            run_result.results = Transient_ops(ckt_template, trans_sim, local_modmap);
        } else if (parser.is_dc) {
            run_result.simulation_type = "dc";
            dc::DCSimulator dcSim = dc::DCsetup(parser, ckt_template);
            run_result.results = dc::DC_ops(ckt_template, dcSim, local_modmap);
        } else if (parser.is_ac) {
            run_result.simulation_type = "ac";
            CKTloadAC(ckt_template);
            ckt_template.cktdematrix->set_init_cxmatrix();
            AC::ACsimulator acSim = AC::ACsetup(parser, ckt_template);
            run_result.results = AC::AC_ops(ckt_template, acSim, local_modmap);
        } else {
            // Fallback or error for no simulation type specified
            run_result.simulation_type = "none";
        }
        
        run_result.ckt = std::move(ckt_template); // Store the circuit state in the result
        run_result.success = true; // Mark as successful
        
    } catch (const SimulationException& e) {
        run_result.error_message = e.what();
        run_result.error_type = e.get_error_type();
        // Fail fast - return immediately to free resources for other tasks
        return run_result;
    } catch (const std::exception& e) {
        run_result.error_message = std::string("Unexpected error: ") + e.what();
        run_result.error_type = "UNKNOWN_ERROR";
        return run_result;
    }
    
    return run_result;
}


// Main function to orchestrate the entire batch simulation process.
std::vector<BatchRunResult> run_batch_simulation(
    const Circuitmap &cktmap,
    const CircuitParser& parser,
    const Modelmap& modmap)
{
    // 1. Find all parameters specified with batch syntax
    auto batch_params = find_batch_params(parser.elements);
    if (batch_params.empty()) {
        std::cout << "No batch parameters found. Nothing to do." << std::endl;
        return {};
    }

    // 2. Generate all unique circuit configurations from the batch parameters
    auto configs = generate_all_configs(batch_params);
    std::cout << "Generated " << configs.size() << " circuit configurations for batch simulation." << std::endl;

    // 3. Set up the thread pool and submit each configuration as a separate task
    BS::thread_pool pool;
    BS::multi_future<BatchRunResult> futures;
    futures.reserve(configs.size()); // Reserve space for the results
    for (const auto& config : configs) {
        // Submit one task for each configuration
        futures.push_back(pool.submit_task(
            [&cktmap, &parser, &modmap, config]() { // Pass const references where possible
                // Inside the lambda, call the worker with its arguments.
                return simulation_worker(cktmap, parser, modmap, config);
            }
        ));
    }
    futures.wait(); // Wait for all tasks to complete

    // Single-threaded execution example (uncomment to use):
    // std::vector<BatchRunResult> single_threaded_results;
    // single_threaded_results.reserve(configs.size());
    // for (const auto& config : configs) {
    //     // Call the worker function directly for single-threaded execution
    //     single_threaded_results.push_back(
    //         simulation_worker(cktmap, parser, modmap, config)
    //     );
    // }

    // 4. Wait for all simulations to complete and return the collected results
    std::cout << "All batch simulation tasks submitted. Waiting for completion..." << std::endl;
    return futures.get();
}


} // namespace batch
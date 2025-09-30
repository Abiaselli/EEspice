#pragma once
#include <iostream>
#include <vector>
#include <future>
#include "CKT.hpp"
#include "gal_variables.hpp"
#include "device.hpp"
#include "bsim4v82/bsim4v82load/bsim4v82load.hpp"

/**
 * @struct ThreadBuf
 * @brief Thread-local storage for BSIM4 parallel processing results.
 * Each thread accumulates its results in its own buffer, avoiding false sharing
 * and reducing upfront memory allocation compared to a pre-sized global vector.
 */
struct alignas(std::hardware_destructive_interference_size) ThreadBuf {
    std::vector<size_t> idx;                    // Device indices processed by this thread
    std::vector<bsim4::BSIM4V82> states;        // Updated instance states
    std::vector<bsim4::BSIM4stamp> stamps;      // Computed stamps
};
thread_local ThreadBuf tls_buffer;


std::pair<arma::mat, arma::vec> NonLinearMutithreaded(CKTcircuit &ckt, const arma::vec &pre_NR_solution, 
    const std::pair<arma::mat, arma::vec> &matrixes, const Modelmap &modmap, double h)
{
    ScopedTimer timer(ckt.sim_stats.simTime.matrix_load_time);
    arma::mat LHS = matrixes.first;
    
    arma::vec RHS = matrixes.second;

    for (const auto &nmos : ckt.CKTelements.nmos) {
        const NMOSModel nmosModel = modmap.nmosModels.at(nmos.modelName);
        NMOS_assigner(nmos.id, nmos.node_vd, nmos.node_vg, nmos.node_vs, nmos.node_vb, nmos.W, nmos.L, pre_NR_solution, ckt.T_nodes, LHS, RHS, nmosModel);
    }
    for (const auto &pmos : ckt.CKTelements.pmos) {
        const PMOSModel pmosModel = modmap.pmosModels.at(pmos.modelName);
        PMOS_assigner(pmos.id, pmos.node_vd, pmos.node_vg, pmos.node_vs, pmos.node_vb, pmos.W, pmos.L, pre_NR_solution, ckt.T_nodes, LHS, RHS, pmosModel);
    }

    // Parallel BSIM4: compute stamps in parallel using TLS with chunked task submission
    if (!ckt.CKTelements.bsim4.empty()) {
        const size_t num_devices = ckt.CKTelements.bsim4.size();
        const size_t num_threads = pool_global.get_thread_count();

        // Thread-safe collection of results from all threads
        std::vector<ThreadBuf> thread_results(num_threads);

        // Step 1: Submit chunked tasks to thread pool
        for (size_t thread_id = 0; thread_id < num_threads; ++thread_id) {
            pool_global.detach_task([&, thread_id]() {
                // Clear thread-local buffers at the start of each task
                tls_buffer.idx.clear();
                tls_buffer.states.clear();
                tls_buffer.stamps.clear();

                // Calculate chunk boundaries for this thread
                const size_t chunk_start = thread_id * num_devices / num_threads;
                const size_t chunk_end = (thread_id + 1) * num_devices / num_threads;

                // reserve tls_buffer
                tls_buffer.idx.reserve(chunk_end - chunk_start);
                tls_buffer.states.reserve(chunk_end - chunk_start);
                tls_buffer.stamps.reserve(chunk_end - chunk_start);

                // Process all devices in this thread's chunk
                for (size_t i = chunk_start; i < chunk_end; ++i) {
                    const bsim4::BSIM4model &b4model = *ckt.CKTelements.bsim4[i].bsim4v82Instance.BSIM4modPtr;
                    // Create a local copy of the instance state
                    bsim4::BSIM4V82 local_instance = ckt.CKTelements.bsim4[i].bsim4v82Instance;

                    // Calculate stamps
                    bsim4::BSIM4stamp stamp = bsim4::BSIM4calculateStamps(
                        ckt, b4model, local_instance, ckt.spiceCompatible,
                        pre_NR_solution, ckt.CKTtemp, ckt.CKTgmin);

                    // Accumulate results in thread-local storage
                    tls_buffer.idx.push_back(i);
                    tls_buffer.states.push_back(local_instance);
                    tls_buffer.stamps.push_back(stamp);
                }

                // Copy thread-local results to the shared collection
                thread_results[thread_id] = std::move(tls_buffer);
            });
        }

        // Step 2: Wait for all calculations to complete
        pool_global.wait();

        // Step 3: Merge results from all thread buffers and apply stamps serially
        // Process threads in order to maintain device ordering
        for (size_t thread_id = 0; thread_id < num_threads; ++thread_id) {
            const ThreadBuf &buf = thread_results[thread_id];

            // Apply all results from this thread
            for (size_t j = 0; j < buf.idx.size(); ++j) {
                size_t device_idx = buf.idx[j];
                // Update the original instance with the modified state
                ckt.CKTelements.bsim4[device_idx].bsim4v82Instance = buf.states[j];
                // Apply the stamps to the matrices
                bsim4::bsim4applyStamps(ckt.CKTelements.bsim4[device_idx].bsim4v82Instance,
                                       buf.stamps[j], LHS, RHS);
            }
        }
    }


    for (const auto &diode : ckt.CKTelements.diodes){
        Diode_assigner(diode.nodePos, diode.nodeNeg, diode.Is, diode.VT, LHS, RHS, pre_NR_solution);
    }
    return {LHS, RHS};
}
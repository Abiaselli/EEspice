#pragma once
#include <iostream>
#include <vector>
#include <future>
#include <mutex>
#include "CKT.hpp"
#include "gal_variables.hpp"
#include "device.hpp"
#include "bsim4v82/bsim4v82load/bsim4v82load.hpp"
#include "bsim4v82/bsim4v82load/bsim4v82applyStampsatomic.hpp"


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

    // Parallel processing of BSIM4 transistors with direct matrix updates
    // if (!ckt.CKTelements.bsim4.empty()) {
    //     std::mutex matrix_mutex;

    //     // Use detach_task to avoid the overhead involved in assigning a future to the task, in order to increase performance

    //     for (size_t i = 0; i < ckt.CKTelements.bsim4.size(); ++i) {
    //         pool_global.detach_task([&, i]() -> void {
    //             const bsim4::BSIM4model &b4model = *ckt.CKTelements.bsim4[i].bsim4v82Instance.BSIM4modPtr;
    //             bsim4::BSIM4V82 &b4instance = ckt.CKTelements.bsim4[i].bsim4v82Instance;
    //             // Calculate all the values needed for the matrices.
    //             const auto stamps = bsim4::BSIM4calculateStamps(ckt, b4model, b4instance, ckt.spiceCompatible, pre_NR_solution, ckt.CKTtemp, ckt.CKTgmin);
    //             // Lock and directly modify shared matrices
    //             std::lock_guard<std::mutex> lock(matrix_mutex);
    //             // Apply the calculated values to the LHS and RHS matrices.
    //             bsim4::bsim4applyStampsfetch(b4instance, stamps, LHS, RHS);
    //         });
    //     }
    //     pool_global.wait();
    // }

    // Parallel BSIM4: compute stamps in parallel with detach_task
    // if (!ckt.CKTelements.bsim4.empty()) {
    //     // Step 1: Create a vector to hold the results from each thread.
    //     std::vector<bsim4::BSIM4stamp> results(ckt.CKTelements.bsim4.size());

    //     for (size_t i = 0; i < ckt.CKTelements.bsim4.size(); ++i) {
    //         // Submit tasks and store futures
    //         pool_global.detach_task([&, i]() {
    //             const bsim4::BSIM4model &b4model = *ckt.CKTelements.bsim4[i].bsim4v82Instance.BSIM4modPtr;
    //             bsim4::BSIM4V82 &b4instance = ckt.CKTelements.bsim4[i].bsim4v82Instance;
                
    //             // Calculate and store the result in the corresponding slot. NO LOCKING.
    //             results[i] = bsim4::BSIM4calculateStamps(ckt, b4model, b4instance, ckt.spiceCompatible, pre_NR_solution, ckt.CKTtemp, ckt.CKTgmin);
    //         });
    //     }

    //     // Wait for all calculations to complete.
    //     pool_global.wait();

    //     // Step 2: Now, apply all the results serially. This is very fast.
    //     for (size_t i = 0; i < ckt.CKTelements.bsim4.size(); ++i) {
    //         bsim4::bsim4applyStamps(ckt.CKTelements.bsim4[i].bsim4v82Instance, results[i], LHS, RHS);
    //     }
    // }

    // Parallel BSIM4: compute stamps in parallel with detach_loop, then apply serially
    const std::size_t N = ckt.CKTelements.bsim4.size();
    if (N > 0) {
        // 1) pick a good num_blocks: ~threads^2/2, clamped to [threads, N]
        const std::size_t threads = pool_global.get_thread_count();
        const std::size_t heuristic = (threads * threads) / 2; // avoids over-blocking vs. threads^2
        // const std::size_t num_blocks =
        //     std::max<std::size_t>(threads, std::min<std::size_t>(heuristic, N));
        const std::size_t num_blocks = 2;

        // 2) compute stamps in parallel, one index per task; each writes to its own slot
        std::vector<bsim4::BSIM4stamp> results(N);
        pool_global.detach_loop(
            std::size_t{0}, N,
            [&](std::size_t i) {
                const auto &dev = ckt.CKTelements.bsim4[i];
                const bsim4::BSIM4model &b4model = *dev.bsim4v82Instance.BSIM4modPtr;
                bsim4::BSIM4V82 &b4instance = ckt.CKTelements.bsim4[i].bsim4v82Instance;

                results[i] = bsim4::BSIM4calculateStamps(
                    ckt, b4model, b4instance, ckt.spiceCompatible,
                    pre_NR_solution, ckt.CKTtemp, ckt.CKTgmin);
            },
            num_blocks);

        // 3) wait for all per-device computations to finish
        pool_global.wait();

        // 4) apply stamps to the global matrices (fast, serial; avoids lock contention)
        for (std::size_t i = 0; i < N; ++i) {
            bsim4::bsim4applyStamps(ckt.CKTelements.bsim4[i].bsim4v82Instance, results[i], LHS, RHS);
        }
    }

    // // Parallel BSIM4 using detach_blocks
    // const std::size_t N = ckt.CKTelements.bsim4.size();
    // if (N > 0) {
    //     const std::size_t threads   = pool_global.get_thread_count();
    //     const std::size_t heuristic = (threads * threads) / 2;            // ~threads^2/2
    //     const std::size_t num_blocks =
    //         std::max<std::size_t>(threads, std::min<std::size_t>(heuristic, N));
        
    //     // Per-device results (one slot per device, avoids locking)
    //     std::vector<bsim4::BSIM4stamp> results(N);

    //     // block(start,end) runs once per block; we loop i inside the block
    //     pool_global.detach_blocks(
    //         std::size_t{0}, N,
    //         [&](std::size_t start, std::size_t end) {
    //             for (std::size_t i = start; i < end; ++i) {
    //                 const auto &dev = ckt.CKTelements.bsim4[i];
    //                 const bsim4::BSIM4model &b4model = *dev.bsim4v82Instance.BSIM4modPtr;
    //                 bsim4::BSIM4V82 &b4instance = ckt.CKTelements.bsim4[i].bsim4v82Instance;

    //                 results[i] = bsim4::BSIM4calculateStamps(
    //                     ckt, b4model, b4instance, ckt.spiceCompatible,
    //                     pre_NR_solution, ckt.CKTtemp, ckt.CKTgmin);
    //             }
    //         },
    //         num_blocks
    //     );

    //     pool_global.wait(); // wait for all blocks

    //     // Apply stamps to global matrices (serial: no contention)
    //     for (std::size_t i = 0; i < N; ++i) {
    //         bsim4::bsim4applyStamps(ckt.CKTelements.bsim4[i].bsim4v82Instance, results[i], LHS, RHS);
    //     }
    // }

    // Parallel BSIM4 using detach_sequence (one task per index)
    // const std::size_t N = ckt.CKTelements.bsim4.size();
    // if (N > 0) {
    //     std::vector<bsim4::BSIM4stamp> results(N);

    //     pool_global.detach_sequence(
    //         std::size_t{0}, N,
    //         [&](std::size_t i) {
    //             const auto &dev = ckt.CKTelements.bsim4[i];
    //             const bsim4::BSIM4model &b4model = *dev.bsim4v82Instance.BSIM4modPtr;
    //             bsim4::BSIM4V82 &b4instance = ckt.CKTelements.bsim4[i].bsim4v82Instance;

    //             results[i] = bsim4::BSIM4calculateStamps(
    //                 ckt, b4model, b4instance, ckt.spiceCompatible,
    //                 pre_NR_solution, ckt.CKTtemp, ckt.CKTgmin);
    //         }
    //     );

    //     // Must wait for all detached tasks to complete before applying
    //     pool_global.wait();  // required when using detach_*()

    //     for (std::size_t i = 0; i < N; ++i) {
    //         bsim4::bsim4applyStamps(ckt.CKTelements.bsim4[i].bsim4v82Instance, results[i], LHS, RHS);
    //     }
    // }

    // Parallel BSIM4 using detach_task with chunking (one task per thread)
    // if (!ckt.CKTelements.bsim4.empty()) {
    //     const size_t num_transistors = ckt.CKTelements.bsim4.size();
    //     const size_t num_threads = pool_global.get_thread_count();
    //     std::vector<bsim4::BSIM4stamp> results(num_transistors);

    //     // Create one task for each thread.
    //     for (size_t i = 0; i < num_threads; ++i) {
    //         pool_global.detach_task([&, i]() {
    //             // Calculate the start and end index for this thread's chunk
    //             const size_t start = i * num_transistors / num_threads;
    //             const size_t end = (i + 1) * num_transistors / num_threads;

    //             // Each thread processes its own large chunk of work
    //             for (size_t j = start; j < end; ++j) {
    //                 const auto& dev = ckt.CKTelements.bsim4[j];
    //                 bsim4::BSIM4V82& b4instance = ckt.CKTelements.bsim4[j].bsim4v82Instance;
    //                 const bsim4::BSIM4model& b4model = *dev.bsim4v82Instance.BSIM4modPtr;
                    
    //                 // No locking needed as each thread writes to a unique part of the vector
    //                 results[j] = bsim4::BSIM4calculateStamps(
    //                     ckt, b4model, b4instance, ckt.spiceCompatible,
    //                     pre_NR_solution, ckt.CKTtemp, ckt.CKTgmin);
    //             }
    //         });
    //     }

    //     // Wait for all chunk-processing tasks to complete.
    //     pool_global.wait();

    //     // Apply results serially (this is fast).
    //     for (size_t i = 0; i < num_transistors; ++i) {
    //         bsim4::bsim4applyStamps(ckt.CKTelements.bsim4[i].bsim4v82Instance, results[i], LHS, RHS);
    //     }
    // }


    for (const auto &diode : ckt.CKTelements.diodes){
        Diode_assigner(diode.nodePos, diode.nodeNeg, diode.Is, diode.VT, LHS, RHS, pre_NR_solution);
    }
    return {LHS, RHS};
}
#pragma once
#include "sim_variables.hpp"
#include "SPICEcompatible.hpp"
#include "devsup.hpp"
#include "bsim4v82/bsim4v82.hpp"
#include "bsim4v82/bsim4v82const.hpp"
#include "bsim4v82/bsim4v82NI.hpp"
#include "bsim4v82calculateStamps.hpp"
#include "bsim4v82applyStamps.hpp"
#include "hybrid_matrix.hpp"

#include <cmath>
#include <armadillo>

namespace bsim4{
int
loadomp(CKTcircuit &ckt, const arma::vec &pre_NR_solution,
                      HybridMatrix &LHS, arma::vec &RHS,
                      std::vector<BSIM4stamp> &stamps)
{
    // Parallel BSIM4: compute stamps in parallel using OpenMP (one device per thread)
    if (!ckt.CKTelements.bsim4.empty()) {
        // Parallel computation - each iteration processes one device
        ScopedTimer bsim4_timer(ckt.sim_stats.simTime.bsim4_time);
        #pragma omp parallel for
        for (size_t i = 0; i < ckt.CKTelements.bsim4.size(); ++i) {
            const bsim4::BSIM4model &b4model = *ckt.CKTelements.bsim4[i].bsim4v82Instance.BSIM4modPtr;
            bsim4::BSIM4V82 &instance = ckt.CKTelements.bsim4[i].bsim4v82Instance;
            // Calculate stamps
            stamps[i] = bsim4::BSIM4calculateStamps(ckt, b4model, instance, ckt.spiceCompatible, pre_NR_solution, ckt.CKTtemp, ckt.CKTgmin);
        }
        for (size_t i = 0; i < ckt.CKTelements.bsim4.size(); ++i) {
            bsim4::BSIM4V82 &instance = ckt.CKTelements.bsim4[i].bsim4v82Instance;
            bsim4::bsim4applyStamps(instance, stamps[i], LHS, RHS);
        }
    }
    return 0; // return success
}
int
loadompColor(CKTcircuit &ckt, const arma::vec &pre_NR_solution,
                      HybridMatrix &LHS, arma::vec &RHS,
                      std::vector<BSIM4stamp> &stamps,
                      const BSIM4Coloring &coloring)
{
    if (!ckt.CKTelements.bsim4.empty()) {
        ScopedTimer bsim4_timer(ckt.sim_stats.simTime.bsim4_time);
        const size_t n = ckt.CKTelements.bsim4.size();
        
        // Phase 1: Parallel computation of stamps
        auto start_calc = std::chrono::high_resolution_clock::now();
        #pragma omp parallel for
        for (size_t i = 0; i < n; ++i) {
            const BSIM4model &b4model = *ckt.CKTelements.bsim4[i].bsim4v82Instance.BSIM4modPtr;
            BSIM4V82 &instance = ckt.CKTelements.bsim4[i].bsim4v82Instance;
            stamps[i] = BSIM4calculateStamps(ckt, b4model, instance, 
                                                     ckt.spiceCompatible, pre_NR_solution, 
                                                     ckt.CKTtemp, ckt.CKTgmin);
        }
        auto end_calc = std::chrono::high_resolution_clock::now();
        // timing.parallel_calc_time = std::chrono::duration<double>(end_calc - start_calc).count();
        
        // Graph coloring
        const auto& color_groups = coloring.getColorGroups();

        // Phase 2: Parallel application of stamps using coloring
        auto start_apply = std::chrono::high_resolution_clock::now();
        for (const auto& group : color_groups) {
            // All instances in this group can be processed in parallel
            #pragma omp parallel for
            for (size_t idx = 0; idx < group.size(); ++idx) {
                size_t i = group[idx];
                BSIM4V82 &instance = ckt.CKTelements.bsim4[i].bsim4v82Instance;
                bsim4applyStamps(instance, stamps[i], LHS, RHS);
            }
        }
        auto end_apply = std::chrono::high_resolution_clock::now();
        // timing.apply_stamps_time = std::chrono::duration<double>(end_apply - start_apply).count();
    }
    return 0; // return success
}

int
loadompColor2(CKTcircuit &ckt, const arma::vec &pre_NR_solution,
                      HybridMatrix &LHS, arma::vec &RHS,
                      std::vector<BSIM4stamp> &stamps,
                      const BSIM4Coloring &coloring)
{
    if (!ckt.CKTelements.bsim4.empty()) {
        ScopedTimer bsim4_timer(ckt.sim_stats.simTime.bsim4_time);
        const size_t n = ckt.CKTelements.bsim4.size();
        
        // Phase 1: Parallel computation of stamps
        auto start_calc = std::chrono::high_resolution_clock::now();
        #pragma omp parallel for
        for (size_t i = 0; i < n; ++i) {
            const BSIM4model &b4model = *ckt.CKTelements.bsim4[i].bsim4v82Instance.BSIM4modPtr;
            BSIM4V82 &instance = ckt.CKTelements.bsim4[i].bsim4v82Instance;
            stamps[i] = BSIM4calculateStamps(ckt, b4model, instance, 
                                                     ckt.spiceCompatible, pre_NR_solution, 
                                                     ckt.CKTtemp, ckt.CKTgmin);
        }
        auto end_calc = std::chrono::high_resolution_clock::now();
        // timing.parallel_calc_time = std::chrono::duration<double>(end_calc - start_calc).count();
        
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
                    BSIM4V82 &instance = ckt.CKTelements.bsim4[i].bsim4v82Instance;
                    bsim4applyStamps(instance, stamps[i], LHS, RHS);
                }
            }
        }
        auto end_apply = std::chrono::high_resolution_clock::now();
        // timing.apply_stamps_time = std::chrono::duration<double>(end_apply - start_apply).count();
    }
    return 0; // return success
}

int
loadompColor3(CKTcircuit &ckt, const arma::vec &pre_NR_solution,
                      HybridMatrix &LHS, arma::vec &RHS,
                      std::vector<BSIM4stamp> &stamps,
                      const BSIM4Coloring &coloring)
{
    if (!ckt.CKTelements.bsim4.empty()) {
        ScopedTimer bsim4_timer(ckt.sim_stats.simTime.bsim4_time);
        const size_t n = ckt.CKTelements.bsim4.size();
        
        // Phase 1: Parallel computation of stamps
        auto start_calc = std::chrono::high_resolution_clock::now();
        #pragma omp parallel for schedule(dynamic)
        for (size_t i = 0; i < n; ++i) {
            const BSIM4model &b4model = *ckt.CKTelements.bsim4[i].bsim4v82Instance.BSIM4modPtr;
            BSIM4V82 &instance = ckt.CKTelements.bsim4[i].bsim4v82Instance;
            stamps[i] = BSIM4calculateStamps(ckt, b4model, instance, 
                                                     ckt.spiceCompatible, pre_NR_solution, 
                                                     ckt.CKTtemp, ckt.CKTgmin);
        }
        auto end_calc = std::chrono::high_resolution_clock::now();
        // timing.parallel_calc_time = std::chrono::duration<double>(end_calc - start_calc).count();
        
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
                    BSIM4V82 &instance = ckt.CKTelements.bsim4[i].bsim4v82Instance;
                    bsim4applyStamps(instance, stamps[i], LHS, RHS);
                }
            }
        }
        auto end_apply = std::chrono::high_resolution_clock::now();
        // timing.apply_stamps_time = std::chrono::duration<double>(end_apply - start_apply).count();
    }
    return 0; // return success
}

int
loadompColor4(CKTcircuit &ckt, const arma::vec &pre_NR_solution,
                      HybridMatrix &LHS, arma::vec &RHS,
                      const BSIM4Coloring &coloring)
{
    if (!ckt.CKTelements.bsim4.empty()) {
        ScopedTimer bsim4_timer(ckt.sim_stats.simTime.bsim4_time);

        // Graph coloring
        const auto& color_groups = coloring.getColorGroups();

        for (const auto& group : color_groups) {
            // All instances in this group can be processed in parallel
            #pragma omp parallel for
            for (size_t idx = 0; idx < group.size(); ++idx) {
                size_t i = group[idx];
                const BSIM4model &b4model = *ckt.CKTelements.bsim4[i].bsim4v82Instance.BSIM4modPtr;
                BSIM4V82 &instance = ckt.CKTelements.bsim4[i].bsim4v82Instance;
                // Phase 1: Parallel computation of stamps
                const BSIM4stamp stamp = BSIM4calculateStamps(ckt, b4model, instance, 
                                                     ckt.spiceCompatible, pre_NR_solution, 
                                                     ckt.CKTtemp, ckt.CKTgmin);
                // Phase 2: Apply stamps
                bsim4applyStamps(instance, stamp, LHS, RHS);
            }
        }
    }
    return 0; // return success
}

} // namespace bsim4
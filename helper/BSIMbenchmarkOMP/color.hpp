#pragma once
#include "CKT.hpp"
#include "hybrid_matrix.hpp"
#include "bsim4v82/bsim4v82load/bsim4v82calculateStamps.hpp"
#include "bsim4v82/bsim4v82load/bsim4v82applyStamps.hpp"
#include <omp.h>
#include <chrono>
#include "time.hpp"

// BSIM4Coloring class is already defined in src/color.hpp (included via CKT.hpp)

LoadOMPTiming loadompColor(CKTcircuit &ckt, const arma::vec &pre_NR_solution,
                      HybridMatrix &LHS, arma::vec &RHS,
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
                      HybridMatrix &LHS, arma::vec &RHS,
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
                      HybridMatrix &LHS, arma::vec &RHS,
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
                      HybridMatrix &LHS, arma::vec &RHS,
                      const BSIM4Coloring &coloring)
{
    LoadOMPTiming timing{0.0, 0.0};

    if (!ckt.CKTelements.bsim4.empty()) {
        // Graph coloring
        const auto& color_groups = coloring.getColorGroups();

        // Fused calculation and application (timed together)
        auto start_fused = std::chrono::high_resolution_clock::now();
        for (const auto& group : color_groups) {
            // All instances in this group can be processed in parallel
            #pragma omp parallel for
            for (size_t idx = 0; idx < group.size(); ++idx) {
                size_t i = group[idx];
                const bsim4::BSIM4model &b4model = *ckt.CKTelements.bsim4[i].bsim4v82Instance.BSIM4modPtr;
                bsim4::BSIM4V82 &instance = ckt.CKTelements.bsim4[i].bsim4v82Instance;
                // Fused: Calculate stamp and apply immediately
                const bsim4::BSIM4stamp stamp = bsim4::BSIM4calculateStamps(ckt, b4model, instance,
                                                     ckt.spiceCompatible, pre_NR_solution,
                                                     ckt.CKTtemp, ckt.CKTgmin);
                bsim4::bsim4applyStamps(instance, stamp, LHS, RHS);
            }
        }
        auto end_fused = std::chrono::high_resolution_clock::now();
        // Report fused time in apply_stamps_time (parallel_calc_time = 0 since no separate phase)
        timing.apply_stamps_time = std::chrono::duration<double>(end_fused - start_fused).count();
    }

    return timing;
}


#pragma once
#include "CKT.hpp"
#include "bsim4v82/bsim4v82load/bsim4v82calculateStamps.hpp"
#include "bsim4v82/bsim4v82load/bsim4v82applyStamps.hpp"
#include <omp.h>

/**
 * @brief Benchmark function for testing OpenMP multi-threaded BSIM4 stamp calculation and MNA application.
 * 
 * This function computes the BSIM4 device stamps in parallel using OpenMP,
 * allowing for efficient multi-threaded execution. Each device's stamp is
 * calculated independently, making it suitable for parallel processing.
 * 
 * @param ckt The circuit containing BSIM4 devices.
 * @param pre_NR_solution The solution vector from the previous Newton-Raphson iteration.
 * @param LHS Left-hand side matrix to be updated with device stamps.
 * @param RHS Right-hand side vector to be updated with device contributions.
 * @param stamps Output vector of calculated BSIM4 stamps (size must match device count).
 */
void loadomp(CKTcircuit &ckt, const arma::vec &pre_NR_solution, arma::mat &LHS, arma::vec &RHS, std::vector<bsim4::BSIM4stamp> &stamps){
    // Parallel BSIM4: compute stamps in parallel using OpenMP (one device per thread)
    if (!ckt.CKTelements.bsim4.empty()) {
        // Parallel computation - each iteration processes one device
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
}
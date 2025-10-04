#pragma once
#include "CKT.hpp"
#include "bsim4v82/bsim4v82load/bsim4v82calculateStamps.hpp"
#include <omp.h>

/** 
 * @brief Benchmark function for testing OpenMP multi-threaded BSIM4 stamp calculation
 * 
 *  Parallelizes stamp computation across BSIM4 devices using OpenMP work-sharing.
 *  Each thread independently calculates nonlinear stamps for one MOSFET instance.
 *  Used for performance testing and validation of parallel stamping strategies.
 * @param ckt Circuit containing BSIM4 device instances and simulation state
 * @param pre_NR_solution Previous Newton-Raphson solution vector for linearization point
 * @param LHS Left-hand side matrix (unused in this benchmark variant)
 * @param RHS Right-hand side vector (unused in this benchmark variant)
 * @param stamps Output vector of calculated BSIM4 stamps (size must match device count)
 */
void calomp(CKTcircuit &ckt, const arma::vec &pre_NR_solution, arma::mat &LHS, arma::vec &RHS, std::vector<bsim4::BSIM4stamp> &stamps){
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
    }
}
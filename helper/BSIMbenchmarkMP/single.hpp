#pragma once
#include "CKT.hpp"
#include "bsim4v82/bsim4v82load/bsim4v82load.hpp"

/**
 * @brief Calculate BSIM4 device stamps in a single-threaded manner.
 * 
 * This function iterates over all BSIM4 devices in the circuit, calculates their
 * stamps based on the current circuit state, and stores the results in the provided
 * stamps vector.
 * 
 * @param ckt Reference to the circuit containing BSIM4 devices
 * @param pre_NR_solution Vector of previous Newton-Raphson solution values
 * @param LHS Left-hand side matrix (unused in this benchmark variant)
 * @param RHS Right-hand side vector (unused in this benchmark variant)
 * @param stamps Output vector of calculated BSIM4 stamps (size must match device count)
 */
void calsingle(CKTcircuit &ckt, const arma::vec &pre_NR_solution, arma::mat &LHS, arma::vec &RHS, std::vector<bsim4::BSIM4stamp> &stamps){
    if (!ckt.CKTelements.bsim4.empty()) {
        const size_t num_devices = ckt.CKTelements.bsim4.size();

        // Step 1: Serial computation - each iteration processes one device
        for (size_t i = 0; i < ckt.CKTelements.bsim4.size(); ++i) {
            const bsim4::BSIM4model &b4model = *ckt.CKTelements.bsim4[i].bsim4v82Instance.BSIM4modPtr;
            // Create a local copy of the instance state
            bsim4::BSIM4V82 &instance = ckt.CKTelements.bsim4[i].bsim4v82Instance;

            // Calculate stamps
            stamps[i] = bsim4::BSIM4calculateStamps(
                ckt, b4model, instance, ckt.spiceCompatible,
                pre_NR_solution, ckt.CKTtemp, ckt.CKTgmin);

        }
    }
}

/**
 * @brief Load BSIM4 device stamps into the global LHS and RHS matrices in a single-threaded manner.
 * 
 * This function iterates over all BSIM4 devices in the circuit, calculates their
 * stamps based on the current circuit state, and directly applies them to the
 * provided LHS and RHS matrices.
 * 
 * @param ckt Reference to the circuit containing BSIM4 devices
 * @param pre_NR_solution Vector of previous Newton-Raphson solution values
 * @param LHS Left-hand side matrix to which device stamps will be applied
 * @param RHS Right-hand side vector to which device contributions will be applied
 */
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
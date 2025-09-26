#pragma once
#include <iostream>
#include <vector>
#include <future>
#include <mutex>
#include "CKT.hpp"
#include "gal_variables.hpp"
#include "device.hpp"
#include "bsim4v82/bsim4v82load/bsim4v82load.hpp"


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
    if (!ckt.CKTelements.bsim4.empty()) {
        std::mutex matrix_mutex;

        // Use detach_task to avoid the overhead involved in assigning a future to the task, in order to increase performance

        for (size_t i = 0; i < ckt.CKTelements.bsim4.size(); ++i) {
            pool_global.detach_task([&, i]() -> void {
                const bsim4::BSIM4model &b4model = *ckt.CKTelements.bsim4[i].bsim4v82Instance.BSIM4modPtr;
                bsim4::BSIM4V82 &b4instance = ckt.CKTelements.bsim4[i].bsim4v82Instance;
                // Calculate all the values needed for the matrices.
                const auto stamps = bsim4::BSIM4calculateStamps(ckt, b4model, b4instance, ckt.spiceCompatible, pre_NR_solution, ckt.CKTtemp, ckt.CKTgmin);
                // Lock and directly modify shared matrices
                std::lock_guard<std::mutex> lock(matrix_mutex);
                // Apply the calculated values to the LHS and RHS matrices.
                bsim4::bsim4applyStamps(b4instance, stamps, LHS, RHS);
            });
        }
        pool_global.wait();
    }

    for (const auto &diode : ckt.CKTelements.diodes){
        Diode_assigner(diode.nodePos, diode.nodeNeg, diode.Is, diode.VT, LHS, RHS, pre_NR_solution);
    }
    return {LHS, RHS};
}
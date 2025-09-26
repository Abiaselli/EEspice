#pragma once
#include <iostream>
#include <vector>
#include <future>
#include <mutex>
#include "CKT.hpp"
#include "gal_variables.hpp"
#include "device.hpp"
#include "bsim4v82/bsim4v82load/bsim4v82load.hpp"

// Helper function to safely merge thread-local matrices
inline void merge_matrices(arma::mat &target_LHS, arma::vec &target_RHS,
                          BS::multi_future<std::pair<arma::mat, arma::vec>> &thread_futures)
{
    for (const auto& [lhs_part, rhs_part] : thread_futures.get()) {
        target_LHS += lhs_part;
        target_RHS += rhs_part;
    }
}

std::pair<arma::mat, arma::vec> NonLinearMutithreaded(CKTcircuit &ckt, const arma::vec &pre_NR_solution, 
    const std::pair<arma::mat, arma::vec> &matrixes, const Modelmap &modmap, double h)
{
    ScopedTimer timer(ckt.sim_stats.simTime.matrix_load_time);
    arma::mat LHS = matrixes.first;
    
    arma::vec RHS = matrixes.second;
    int matrix_size = LHS.n_rows;

    for (const auto &nmos : ckt.CKTelements.nmos) {
        const NMOSModel nmosModel = modmap.nmosModels.at(nmos.modelName);
        NMOS_assigner(nmos.id, nmos.node_vd, nmos.node_vg, nmos.node_vs, nmos.node_vb, nmos.W, nmos.L, pre_NR_solution, ckt.T_nodes, LHS, RHS, nmosModel);
    }
    for (const auto &pmos : ckt.CKTelements.pmos) {
        const PMOSModel pmosModel = modmap.pmosModels.at(pmos.modelName);
        PMOS_assigner(pmos.id, pmos.node_vd, pmos.node_vg, pmos.node_vs, pmos.node_vb, pmos.W, pmos.L, pre_NR_solution, ckt.T_nodes, LHS, RHS, pmosModel);
    }

    // Parallel processing of BSIM4 transistors
    if (!ckt.CKTelements.bsim4.empty()) {
        BS::multi_future<std::pair<arma::mat, arma::vec>> futures;
        futures.reserve(ckt.CKTelements.bsim4.size());

        for (auto& bsim4 : ckt.CKTelements.bsim4) {
            futures.push_back(pool_global.submit_task([&bsim4, &ckt, &pre_NR_solution, matrix_size]() -> std::pair<arma::mat, arma::vec> {
                arma::mat local_LHS = arma::zeros<arma::mat>(matrix_size, matrix_size);
                arma::vec local_RHS = arma::zeros<arma::vec>(matrix_size);

                const bsim4::BSIM4model &b4model = *bsim4.bsim4v82Instance.BSIM4modPtr;
                bsim4::BSIM4V82 &b4instance = bsim4.bsim4v82Instance;
                bsim4::BSIM4load(ckt, b4model, b4instance, ckt.spiceCompatible, pre_NR_solution, ckt.CKTtemp, ckt.CKTgmin, local_LHS, local_RHS);

                return std::make_pair(local_LHS, local_RHS);
            }));
        }

        // Merge results directly
        merge_matrices(LHS, RHS, futures);
    }

    for (const auto &diode : ckt.CKTelements.diodes){
        Diode_assigner(diode.nodePos, diode.nodeNeg, diode.Is, diode.VT, LHS, RHS, pre_NR_solution);
    }
    return {LHS, RHS};
}
#pragma once
#include <iostream>
#include <vector>
#include <future>
#include <mutex>
#include "CKT.hpp"
#include "gal_variables.hpp"
#include "device.hpp"
#include "bsim4v82/bsim4v82stamp.hpp"

// Helper function to safely merge thread-local matrices
inline void merge_matrices(arma::mat &target_LHS, arma::vec &target_RHS,
                          const std::vector<std::pair<arma::mat, arma::vec>> &thread_results)
{
    for (const auto &result : thread_results) {
        target_LHS += result.first;
        target_RHS += result.second;
    }
}

std::pair<arma::mat, arma::vec> NonLinearBatch(CKTcircuit &ckt, const arma::vec &pre_NR_solution, 
    const std::pair<arma::mat, arma::vec> &matrixes, const Modelmap &modmap, double h)
{
    ScopedTimer timer(ckt.sim_stats.simTime.matrix_load_time);
    arma::mat LHS = matrixes.first;
    arma::vec RHS = matrixes.second;
    // LHS.print("in_LHS matrix =");
    // RHS.print("RHS matrix =");

    // Parallel processing of NMOS transistors
    if (!ckt.CKTelements.nmos.empty()) {
        std::vector<std::future<std::pair<arma::mat, arma::vec>>> futures;
        futures.reserve(ckt.CKTelements.nmos.size());

        for (const auto& nmos : ckt.CKTelements.nmos) {
            futures.push_back(pool_global.submit_task([&]() -> std::pair<arma::mat, arma::vec> {
                arma::mat local_LHS = arma::zeros<arma::mat>(LHS.n_rows, LHS.n_cols);
                arma::vec local_RHS = arma::zeros<arma::vec>(RHS.n_elem);

                const NMOSModel nmosModel = modmap.nmosModels.at(nmos.modelName);
                NMOS_assigner(nmos.id, nmos.node_vd, nmos.node_vg, nmos.node_vs, nmos.node_vb, nmos.W, nmos.L, pre_NR_solution, ckt.T_nodes, local_LHS, local_RHS, nmosModel);

                return std::make_pair(local_LHS, local_RHS);
            }));
        }

        // Collect and merge results
        std::vector<std::pair<arma::mat, arma::vec>> nmos_results;
        nmos_results.reserve(futures.size());
        for (auto &future : futures) {
            nmos_results.push_back(future.get());
        }
        merge_matrices(LHS, RHS, nmos_results);
    }
    // Parallel processing of PMOS transistors
    if (!ckt.CKTelements.pmos.empty()) {
        std::vector<std::future<std::pair<arma::mat, arma::vec>>> futures;
        futures.reserve(ckt.CKTelements.pmos.size());

        for (const auto& pmos : ckt.CKTelements.pmos) {
            futures.push_back(pool_global.submit_task([&]() -> std::pair<arma::mat, arma::vec> {
                arma::mat local_LHS = arma::zeros<arma::mat>(LHS.n_rows, LHS.n_cols);
                arma::vec local_RHS = arma::zeros<arma::vec>(RHS.n_elem);

                const PMOSModel pmosModel = modmap.pmosModels.at(pmos.modelName);
                PMOS_assigner(pmos.id, pmos.node_vd, pmos.node_vg, pmos.node_vs, pmos.node_vb, pmos.W, pmos.L, pre_NR_solution, ckt.T_nodes, local_LHS, local_RHS, pmosModel);

                return std::make_pair(local_LHS, local_RHS);
            }));
        }

        // Collect and merge results
        std::vector<std::pair<arma::mat, arma::vec>> pmos_results;
        pmos_results.reserve(futures.size());
        for (auto &future : futures) {
            pmos_results.push_back(future.get());
        }
        merge_matrices(LHS, RHS, pmos_results);
    }

    // Parallel processing of BSIM4 transistors
    if (!ckt.CKTelements.bsim4.empty()) {
        std::vector<std::future<std::pair<arma::mat, arma::vec>>> futures;
        futures.reserve(ckt.CKTelements.bsim4.size());

        for (const auto& bsim4 : ckt.CKTelements.bsim4) {
            futures.push_back(pool_global.submit_task([&]() -> std::pair<arma::mat, arma::vec> {
                arma::mat local_LHS = arma::zeros<arma::mat>(LHS.n_rows, LHS.n_cols);
                arma::vec local_RHS = arma::zeros<arma::vec>(RHS.n_elem);

                const bsim4::BSIM4model &b4model = *bsim4.bsim4v82Instance.BSIM4modPtr;
                bsim4::BSIM4V82 &b4instance = const_cast<bsim4::BSIM4V82&>(bsim4.bsim4v82Instance);
                bsim4::BSIM4load(ckt, b4model, b4instance, ckt.spiceCompatible, pre_NR_solution, ckt.CKTtemp, ckt.CKTgmin, local_LHS, local_RHS);

                return std::make_pair(local_LHS, local_RHS);
            }));
        }

        // Collect and merge results
        std::vector<std::pair<arma::mat, arma::vec>> bsim4_results;
        bsim4_results.reserve(futures.size());
        for (auto &future : futures) {
            bsim4_results.push_back(future.get());
        }
        merge_matrices(LHS, RHS, bsim4_results);
    }

    for (const auto &diode : ckt.CKTelements.diodes){
        Diode_assigner(diode.nodePos, diode.nodeNeg, diode.Is, diode.VT, LHS, RHS, pre_NR_solution);
    }
    return {LHS, RHS};
}
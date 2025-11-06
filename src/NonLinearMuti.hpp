#pragma once
#include <iostream>
#include <vector>
#include <future>
#include "CKT.hpp"
#include "device.hpp"
#include "bsim4v82/bsim4v82load/bsim4v82load.hpp"
#include <omp.h>


std::pair<arma::mat, arma::vec> NonLinearMutithreaded(CKTcircuit &ckt, const arma::vec &pre_NR_solution, 
    const std::pair<arma::mat, arma::vec> &matrixes, const Modelmap &modmap, double h)
{
    // std::vector<bsim4::BSIM4stamp> stamps(ckt.CKTelements.bsim4.size());
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

    // Parallel BSIM4: compute stamps in parallel using OpenMP (one device per thread)
    if (!ckt.CKTelements.bsim4.empty()) {
        const size_t n = ckt.CKTelements.bsim4.size();
        std::vector<bsim4::BSIM4stamp> stamps(n);
        
        // Phase 1: Parallel computation of stamps
        #pragma omp parallel for
        for (size_t i = 0; i < n; ++i) {
            const bsim4::BSIM4model &b4model = *ckt.CKTelements.bsim4[i].bsim4v82Instance.BSIM4modPtr;
            bsim4::BSIM4V82 &instance = ckt.CKTelements.bsim4[i].bsim4v82Instance;
            stamps[i] = bsim4::BSIM4calculateStamps(ckt, b4model, instance, 
                                                     ckt.spiceCompatible, pre_NR_solution, 
                                                     ckt.CKTtemp, ckt.CKTgmin);
        }
        
        // Graph coloring
        const auto& color_groups = ckt.b4coloring.getColorGroups();

        // Phase 2: Parallel application of stamps using coloring
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
    }


    for (const auto &diode : ckt.CKTelements.diodes){
        Diode_assigner(diode.nodePos, diode.nodeNeg, diode.Is, diode.VT, LHS, RHS, pre_NR_solution);
    }
    return {LHS, RHS};
}
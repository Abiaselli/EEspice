#pragma once
#include <iostream>
#include <vector>
#include <future>
#include "CKT.hpp"
#include "device.hpp"
#include "bsim4v82/bsim4v82load/bsim4v82loadomp.hpp"
#include <omp.h>


/**
 * @brief Multithreaded stamping of nonlinear devices into the circuit matrices
 *
 * This function resets to step baseline (state after dynamic elements stamped),
 * then stamps nonlinear devices directly into ckt.cktmatrix->LHS and RHS.
 * Works for both dense and sparse matrices via the two-level baseline API.
 */
void NonLinearMutithreaded(CKTcircuit &ckt, const arma::vec &pre_NR_solution,
    const Modelmap &modmap, double h)
{
    ScopedTimer timer(ckt.sim_stats.simTime.matrix_load_time);

    // Reset to step baseline (state after dynamic elements, before nonlinear)
    ckt.cktmatrix->LHS.reset_to_step_baseline();
    ckt.cktmatrix->reset_to_step_baseline_RHS();

    // Stamp nonlinear devices directly into ckt.cktmatrix->LHS and RHS
    for (const auto &nmos : ckt.CKTelements.nmos) {
        const NMOSModel nmosModel = modmap.nmosModels.at(nmos.modelName);
        NMOS_assigner(nmos.id, nmos.node_vd, nmos.node_vg, nmos.node_vs, nmos.node_vb, nmos.W, nmos.L, pre_NR_solution, ckt.T_nodes, ckt.cktmatrix->LHS, ckt.cktmatrix->RHS, nmosModel);
    }
    for (const auto &pmos : ckt.CKTelements.pmos) {
        const PMOSModel pmosModel = modmap.pmosModels.at(pmos.modelName);
        PMOS_assigner(pmos.id, pmos.node_vd, pmos.node_vg, pmos.node_vs, pmos.node_vb, pmos.W, pmos.L, pre_NR_solution, ckt.T_nodes, ckt.cktmatrix->LHS, ckt.cktmatrix->RHS, pmosModel);
    }

    // Parallel BSIM4: compute stamps in parallel using OpenMP (one device per thread)
    bsim4::loadompColor4(ckt, pre_NR_solution, ckt.cktmatrix->LHS, ckt.cktmatrix->RHS, ckt.b4coloring);

    for (const auto &diode : ckt.CKTelements.diodes){
        Diode_assigner(diode.nodePos, diode.nodeNeg, diode.Is, diode.VT, ckt.cktmatrix->LHS, ckt.cktmatrix->RHS, pre_NR_solution);
    }
}
#pragma once
#include <iostream>
#include <vector>
#include <future>
#include "CKT.hpp"
#include "device.hpp"
#include "bsim4v82/bsim4v82load/bsim4v82loadomp.hpp"
#include <omp.h>


MNA NonLinearMutithreaded(CKTcircuit &ckt, const arma::vec &pre_NR_solution, 
    const MNA &matrixes, const Modelmap &modmap, double h)
{
    // std::vector<bsim4::BSIM4stamp> stamps(ckt.CKTelements.bsim4.size());
    ScopedTimer timer(ckt.sim_stats.simTime.matrix_load_time);
    MNA mna;
    mna.LHS = matrixes.LHS;
    mna.RHS = matrixes.RHS;

    for (const auto &nmos : ckt.CKTelements.nmos) {
        const NMOSModel nmosModel = modmap.nmosModels.at(nmos.modelName);
        NMOS_assigner(nmos.id, nmos.node_vd, nmos.node_vg, nmos.node_vs, nmos.node_vb, nmos.W, nmos.L, pre_NR_solution, ckt.T_nodes, mna.LHS, mna.RHS, nmosModel);
    }
    for (const auto &pmos : ckt.CKTelements.pmos) {
        const PMOSModel pmosModel = modmap.pmosModels.at(pmos.modelName);
        PMOS_assigner(pmos.id, pmos.node_vd, pmos.node_vg, pmos.node_vs, pmos.node_vb, pmos.W, pmos.L, pre_NR_solution, ckt.T_nodes, mna.LHS, mna.RHS, pmosModel);
    }

    // Parallel BSIM4: compute stamps in parallel using OpenMP (one device per thread)
    bsim4::loadompColor4(ckt, pre_NR_solution, mna.LHS, mna.RHS, ckt.b4coloring);

    for (const auto &diode : ckt.CKTelements.diodes){
        Diode_assigner(diode.nodePos, diode.nodeNeg, diode.Is, diode.VT, mna.LHS, mna.RHS, pre_NR_solution);
    }
    return mna;
}
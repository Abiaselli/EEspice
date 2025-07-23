#pragma once
#include <iostream>
#include "DC_calcs.hpp"
#include "circuit_parser.hpp"

namespace dc{

DCSimulator DCsetup(const CircuitParser &parser, const CKTcircuit &ckt){
    DCSimulator dcSim;

    // Check if the DC simulation is non-linear
    if(!ckt.CKTelements.nmos.empty() || !ckt.CKTelements.pmos.empty() || !ckt.CKTelements.diodes.empty()){
        dcSim.non_linear = true;
    }
    else{
        dcSim.non_linear = false;
    }

    // Setup DC sweeps
    // dcSim.dcsweep = std::move(parser.dcSweep_parser);
    dcSim.dcsweep = parser.dcSweep_parser;

    // setup DC voltage points
    double vol = dcSim.dcsweep.vstart;
    while(vol <= dcSim.dcsweep.vend){
        dcSim.dcsweep.sweep_values.push_back(vol);
        vol += dcSim.dcsweep.vstep;
    }

    if(dcSim.dcsweep.vend - dcSim.dcsweep.sweep_values.back() > LargeEpsilon){  //1e-6
        dcSim.dcsweep.sweep_values.push_back(dcSim.dcsweep.vend);
    }

    return dcSim;
}

void DeviceEvaluation(DC &dc, CKTcircuit &ckt, const DCSimulator &dcSim){
    // Modify matrixes for DC sweep
    // If the source's nodes and are not changed, we can use the same LHS
    // We only need to modify the RHS
    dc.LHS = ckt.cktdematrix->get_init_LHS();
    dc.RHS = ckt.cktdematrix->get_init_RHS();
    dc.solution = arma::vec(ckt.cktdematrix->RHS.n_rows, arma::fill::none);

    for (const auto &vol : ckt.CKTelements.voltageSources){
        if(vol.id_str == dcSim.dcsweep.sourceName){
            dc.RHS(ckt.T_nodes + vol.id - 1, 0) = dc.sweepValue; // Index is starting from 0
            return;
        }
    }
    for (const auto &cs : ckt.CKTelements.currentSources){
        if(cs.id_str == dcSim.dcsweep.sourceName){
            Is_assigner_reverse(cs.nodePos, cs.nodeNeg, cs.value, dc.RHS);
            Is_assigner(cs.nodePos, cs.nodeNeg, dc.sweepValue, dc.RHS);
            return;
        }
    }
    std::cerr << "Error: DC sweep is not supported for this device: " << dcSim.dcsweep.sourceName << std::endl;
    exit(1);
}

arma::vec DC_analysis_once(CKTcircuit &ckt, const DCSimulator &dcSim, DC &dc, const Modelmap &modmap)
{
    // Initialize the DC analysis
    ckt.spiceCompatible.setFlagsDC();
    DeviceEvaluation(dc, ckt, dcSim);
    arma::vec solution(ckt.cktdematrix->RHS.n_rows, arma::fill::zeros);
    // Solve the system
    if(dcSim.non_linear)
    {
        solution = NewtonRaphson_system(ckt, dc.LHS, dc.RHS, modmap);
    }
    else
    {
        solution = arma::solve(dc.LHS, dc.RHS);
    }
    return solution;
}

std::vector<DC> DC_ops(CKTcircuit &ckt, DCSimulator &dcSim, const Modelmap &modmap)
{
    // single sweep loop (only one device)
    const auto &sweepVals = dcSim.dcsweep.sweep_values;
    const auto &sweepName = dcSim.dcsweep.sourceName;
    dcSim.vec_dc.reserve(sweepVals.size());
    
    for (double val : sweepVals) {
        DC dc;
        // store that single sweep value and name
        dc.sweepValue = val;
        dc.sweepName = sweepName;

        // Run one DC analysis
        dc.solution = DC_analysis_once(ckt, dcSim, dc, modmap);

        // Push the result
        dcSim.vec_dc.emplace_back(dc);
    }

    return  dcSim.vec_dc;
}
} // namespace dc
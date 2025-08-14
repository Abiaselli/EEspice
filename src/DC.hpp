#pragma once
#include <iostream>
#include "DC_calcs.hpp"
#include "circuit_parser.hpp"
#include "simulation_exceptions.hpp"

namespace dc{

DCSimulator DCsetup(const CircuitParser &parser, const CKTcircuit &ckt){
    DCSimulator dcSim;

    // Check if the DC simulation is non-linear
    dcSim.non_linear = CKTisNonLinear(ckt.CKTelements);

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

void DeviceEvaluation(DCResult &dc, const CKTcircuit &ckt, const DCSimulator &dcSim, DCMat &dcMat){
    // Modify matrixes for DC sweep
    // If the source's nodes and are not changed, we can use the same LHS
    // We only need to modify the RHS
    dcMat.LHS = ckt.cktdematrix->get_init_LHS();
    dcMat.RHS = ckt.cktdematrix->get_init_RHS();
    dc.solution = arma::vec(ckt.cktdematrix->RHS.n_rows, arma::fill::none);

    for (const auto &vol : ckt.CKTelements.voltageSources){
        if(vol.id_str == dcSim.dcsweep.sourceName){
            auto it = ckt.map.map_branch_currents.find(vol.id_str);
            if (it != ckt.map.map_branch_currents.end()) {
                dcMat.RHS(it->second) = dc.sweepValue;
            }
            return;
        }
    }
    for (const auto &cs : ckt.CKTelements.currentSources){
        if(cs.id_str == dcSim.dcsweep.sourceName){
            Is_assigner_reverse(cs.nodePos, cs.nodeNeg, cs.value, dcMat.RHS);
            Is_assigner(cs.nodePos, cs.nodeNeg, dc.sweepValue, dcMat.RHS);
            return;
        }
    }
    throw SetupException("Error: DC sweep is not supported for this device: " + dcSim.dcsweep.sourceName, "UNSUPPORTED_DC_SWEEP_DEVICE");
}

arma::vec DC_analysis_once(CKTcircuit &ckt, const DCSimulator &dcSim, DCResult &dc, DCMat &dcMat, const Modelmap &modmap)
{
    // Initialize the DC analysis
    ckt.spiceCompatible.setFlagsDC();
    DeviceEvaluation(dc, ckt, dcSim, dcMat);
    arma::vec solution(ckt.cktdematrix->RHS.n_rows, arma::fill::zeros);
    // Solve the system
    if(dcSim.non_linear)
    {
        solution = NewtonRaphson_system(ckt, dcMat.LHS, dcMat.RHS, modmap);
    }
    else
    {
        solution = arma::solve(dcMat.LHS, dcMat.RHS);
    }
    return solution;
}

std::vector<DCResult> DC_ops(CKTcircuit &ckt, DCSimulator &dcSim, const Modelmap &modmap)
{
    // single sweep loop (only one device)
    const auto &sweepVals = dcSim.dcsweep.sweep_values;
    const auto &sweepName = dcSim.dcsweep.sourceName;
    dcSim.vec_dc.reserve(sweepVals.size());

    DCMat dcMat;

    for (double val : sweepVals) {
        DCResult dc;
        // store that single sweep value and name
        dc.sweepValue = val;
        dc.sweepName = sweepName;

        // Run one DC analysis
        dc.solution = DC_analysis_once(ckt, dcSim, dc, dcMat, modmap);

        // Push the result
        dcSim.vec_dc.emplace_back(dc);
    }

    return  dcSim.vec_dc;
}
} // namespace dc
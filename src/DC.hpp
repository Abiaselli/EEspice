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

/**
 * @brief Prepares the matrices for a DC sweep by resetting to linear baseline
 *        and modifying the sweep source value
 *
 * Resets to linear baseline, modifies the sweep source value in RHS,
 * then saves step baseline for NR iterations.
 */
void DeviceEvaluation(DCResult &dc, CKTcircuit &ckt, const DCSimulator &dcSim){
    // Reset to linear baseline (state after static linear elements stamped)
    ckt.cktmatrix->LHS.reset_to_linear_baseline();
    ckt.cktmatrix->reset_to_linear_baseline_RHS();

    dc.solution = arma::vec(ckt.cktmatrix->RHS.n_rows, arma::fill::none);

    // Modify the sweep source value in RHS
    for (const auto &vol : ckt.CKTelements.voltageSources){     // TODO: Skip the linear scan and look up the key directly
        if(vol.id_str == dcSim.dcsweep.sourceName){
            auto it = ckt.map.map_branch_currents.find(vol.id_str);
            if (it != ckt.map.map_branch_currents.end()) {
                ckt.cktmatrix->RHS(it->second) = dc.sweepValue;
            }
            // Save step baseline for NR iterations
            ckt.cktmatrix->LHS.save_step_baseline();
            ckt.cktmatrix->save_step_baseline_RHS();
            return;
        }
    }
    for (const auto &cs : ckt.CKTelements.currentSources){
        if(cs.id_str == dcSim.dcsweep.sourceName){
            Is_assigner_reverse(cs.nodePos, cs.nodeNeg, cs.value, ckt.cktmatrix->RHS);
            Is_assigner(cs.nodePos, cs.nodeNeg, dc.sweepValue, ckt.cktmatrix->RHS);
            // Save step baseline for NR iterations
            ckt.cktmatrix->LHS.save_step_baseline();
            ckt.cktmatrix->save_step_baseline_RHS();
            return;
        }
    }
    throw SetupException("Error: DC sweep is not supported for this device: " + dcSim.dcsweep.sourceName, "UNSUPPORTED_DC_SWEEP_DEVICE");
}

arma::vec DC_analysis_once(CKTcircuit &ckt, const DCSimulator &dcSim, DCResult &dc, const Modelmap &modmap)
{
    // Initialize the DC analysis
    ckt.spiceCompatible.setFlagsDC();
    DeviceEvaluation(dc, ckt, dcSim);
    arma::vec solution(ckt.cktmatrix->RHS.n_rows, arma::fill::zeros);
    // Solve the system
    if(dcSim.non_linear)
    {
        // NonLinear() resets to step baseline, stamps nonlinear devices
        solution = NewtonRaphson_system(ckt, modmap);
    }
    else
    {
        // For linear circuits, step baseline is already set in DeviceEvaluation
        // Just solve using solver()
        solution = solver(ckt.cktmatrix->LHS, ckt.cktmatrix->RHS, ckt);
    }
    return solution;
}

std::vector<DCResult> DC_ops(CKTcircuit &ckt, DCSimulator &dcSim, const Modelmap &modmap)
{
    ScopedTimer analysisTimer(ckt.sim_stats.simTime.analysis_time); // Time the analysis
    // single sweep loop (only one device)
    const auto &sweepVals = dcSim.dcsweep.sweep_values;
    const auto &sweepName = dcSim.dcsweep.sourceName;
    dcSim.vec_dc.reserve(sweepVals.size());

    for (double val : sweepVals) {
        DCResult dc;
        // store that single sweep value and name
        dc.sweepValue = val;
        dc.sweepName = sweepName;

        // Run one DC analysis (uses baseline reset mechanism internally)
        dc.solution = DC_analysis_once(ckt, dcSim, dc, modmap);

        // Push the result
        dcSim.vec_dc.emplace_back(dc);
    }

    ckt.sim_stats.num_data_points = static_cast<int>(dcSim.vec_dc.size());

    return dcSim.vec_dc;
}
} // namespace dc
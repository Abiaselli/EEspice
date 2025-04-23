#pragma once
#include <iostream>
#include "DC_calcs.hpp"

namespace dc{

void DeviceEvaluation(DC &dc, CKTcircuit &ckt, const DCSimulator &dcSim){
    // Modify matrixes for DC sweep
    // If the source's nodes and are not changed, we can use the same LHS
    // We only need to modify the RHS
    dc.LHS = ckt.cktdematrix->get_init_LHS();
    dc.RHS = ckt.cktdematrix->get_init_RHS();
    dc.solution = arma::vec(ckt.cktdematrix->RHS.n_rows, arma::fill::none);

    for(size_t i = 0; i < dcSim.sweeps.size(); ++i){
        bool found = false;
        // find the device by name:
        for(const auto &element : ckt.CKTelements){
            std::visit([&](auto &&arg){
                if(arg.id_str == dcSim.sweeps[i].sourceName){
                    found = true;
                    if constexpr (std::is_same_v<std::decay_t<decltype(arg)>, VoltageSource>){
                        dc.RHS(ckt.T_nodes + arg.id - 1, 0) = dc.sweepValues[i]; // Index is starting from 0
                        found = true;
                    }
                    else if constexpr(std::is_same_v<std::decay_t<decltype(arg)>, CurrentSource>){
                       Is_assigner_reverse(arg.nodePos, arg.nodeNeg, arg.value, dc.RHS);
                       Is_assigner(arg.nodePos, arg.nodeNeg, dc.sweepValues[i], dc.RHS);
                       found = true;
                    }
                    else{
                        std::cerr << "Error: DC sweep is not supported for this device: " << arg.id_str << std::endl;
                        exit(1);
                    }
                }
            }, element.element);
            if (found) break;
        }
    }
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
    
    if(dcSim.sweeps.size() == 1){
        // single sweep loop (only one device)
        const auto &sweepVals = dcSim.sweeps[0].sweep_values;
        const auto &sweepName = dcSim.sweeps[0].sourceName;
        for (double val : sweepVals) {
            DC dc;
            // store that single sweep value and name
            dc.sweepValues.push_back(val);
            dc.sweepNames.push_back(sweepName);

            // Run one DC analysis
            dc.solution = DC_analysis_once(ckt, dcSim, dc, modmap);

            // Push the result
            dcSim.vec_dc.push_back(dc);
        }
    }
    else if(dcSim.sweeps.size() == 2) {
        // double sweep loop (two devices)
        const auto &sweepA = dcSim.sweeps[0].sweep_values;
        const auto &sweepB = dcSim.sweeps[1].sweep_values;
        const auto &nameA      = dcSim.sweeps[0].sourceName;
        const auto &nameB      = dcSim.sweeps[1].sourceName;

        // Nested loops over both sweep arrays
        for (double valA : sweepA) {
            for (double valB : sweepB) {
                DC dc;

                dc.sweepValues = {valA, valB};
                dc.sweepNames  = {nameA, nameB};

                // Run one DC analysis
                dc.solution = DC_analysis_once(ckt, dcSim, dc, modmap);

                dcSim.vec_dc.push_back(dc);
            }
        }
    }
    else{
        std::cerr << "Error: DC sweep is not supported for more than 2 devices" << std::endl;
        exit(1);
    }

    return  dcSim.vec_dc;
}
} // namespace dc
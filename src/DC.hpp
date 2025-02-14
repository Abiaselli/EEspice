#pragma once
#include <iostream>
#include "DC_calcs.hpp"

arma::vec DC_analysis(const CKTcircuit &ckt, const DCSimulator &dc_sim, double vol)
{
    arma::mat LHS = ckt.cktdematrix->get_init_LHS();
    arma::vec RHS = ckt.cktdematrix->get_init_RHS();
    arma::vec solution(ckt.cktdematrix->RHS.n_rows, arma::fill::zeros);

    // Modify matrixes for voltage source
    // If the voltage source's nodes are not changed, we can use the same LHS
    for(VoltageSource vs : dc_sim.vec_voltages)
    {
        if(vs.id_str == dc_sim.dc_config.srcnam)
        {
            RHS(ckt.T_nodes + vs.id - 1, 0) = vol; // Index is starting from 0
        }
    }

    // Solve the system
    if(dc_sim.dc_config.non_linear)
    {
        solution = NewtonRaphson_system(ckt, LHS, RHS);
    }
    else
    {
        solution = arma::solve(LHS, RHS);
    }
    return solution;
}


std::vector<DC> DC_ops(const CKTcircuit &ckt, DCSimulator &dc_sim)
{
    
    for (double vol : dc_sim.voltage_points)
    {
        DC dc;
        dc.vol_point = vol;

        dc.solution = DC_analysis(ckt, dc_sim, vol);

        dc_sim.vec_dc.push_back(dc);
    }
    return  dc_sim.vec_dc;

}
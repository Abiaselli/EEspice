#pragma once

#include "Transient_calcs.hpp"
#include "Transient_multi_solvers.hpp"
#include "Transient_single_solvers.hpp"

Transient Fixed_TimeStep(CKTcircuit &ckt, TransientSimulator &trans_sim, const Modelmap &modmap){
    Transient trans;
    trans.mode = 1; // 0 to do OP analysis, 1 to do transient simulation
    trans.h = trans_sim.trans_config.init_h; // Initial time step
    trans.time_trans = trans_sim.vec_trans.back().time_trans;
    trans.C_list = trans_sim.vec_trans.back().C_list;

    std::cout << "Fixed time step" << std::endl;
    arma::vec solution = NewtonRaphson_system(ckt, trans.h, 1, trans.time_trans, trans.C_list, trans_sim.vec_trans.back().solution, modmap);
    // solution.print("The transient analysis of the circuit is: ");
    ARMA_PRINT(solution, "The transient analysis of the circuit is: ");

    if (trans.C_list.size() > 0)
    {

        auto cur_vol = get_currents_voltages(trans.C_list, trans.h, solution, trans_sim.vec_trans.back().solution);
        trans.Capacitance = get_capacitance(trans.C_list);
        trans.C_current = cur_vol.first;
        trans.C_voltage = cur_vol.second;
        trans.C_charge = trans.C_voltage % trans.Capacitance;
        trans.C_list_update();

        /*--------------------------------------*/
        double current = arma::abs(trans.C_current).max();
        double pre_current = arma::abs(trans_sim.vec_trans.back().C_current).max();
        double charge = arma::abs(trans.C_charge).max();
        double pre_charge = arma::abs(trans_sim.vec_trans.back().C_charge).max();
        double LTE_B_C = RELTOL * std::max(current, pre_current) + ABSTOL;
        double LTE_B_Q = RELTOL * std::max({charge, pre_charge, CHGTOL}) / trans.h;
        double LTE_bound = std::max(LTE_B_C, LTE_B_Q);

        double LTE = trans.h / (2 * trans.C_list.at(0).value) * std::abs((current - pre_current));

        // if(LTE > 7 * LTE_bound){
        //     std::cout << "LTE > LTE_bound" << std::endl;
        // }
    }
    trans.time_trans = trans_sim.vec_trans.back().time_trans + trans.h;
    trans.solution = solution;

    trans.trans_count = trans_sim.vec_trans.back().trans_count + 1;
    return trans;
}

Transient Varibale_TimeStep(CKTcircuit &ckt, TransientSimulator &trans_sim, const Modelmap &modmap){
    Transient trans;
    trans.mode = 1; // 0 to do OP analysis, 1 to do transient simulation

    if (trans_sim.vec_trans.size() == 1)
    {
        // The first step of the transient simulation
        trans.C_list = trans_sim.vec_trans.back().C_list;
        trans.h = trans_sim.trans_config.init_h; // Initial time step
        trans.time_trans = trans_sim.trans_config.t_start + trans.h;

        trans.solution = NewtonRaphson_system(ckt, trans.h, 1, trans.time_trans, trans.C_list, trans_sim.vec_trans.back().solution, modmap);
        // solution.print("The transient analysis of the circuit is: ");
        ARMA_PRINT(trans.solution, "The transient analysis of the circuit is: ");

        if (trans.C_list.size() > 0)
        {
            std::pair<arma::vec, arma::vec> cur_vol = get_currents_voltages(trans.C_list, trans.h, trans.solution, trans_sim.vec_trans.back().solution);
            trans.Capacitance = get_capacitance(trans.C_list);

            trans.C_current = cur_vol.first;
            trans.C_voltage = cur_vol.second;
            trans.C_charge = trans.C_voltage % trans.Capacitance;
            trans.C_list_update();
        }

        trans.trans_count += 1;
        trans.next_h = trans.h;
    }

    else
    {
        trans.C_list = trans_sim.vec_trans.back().C_list;
        // Just to prevent single_next_h from directly using time_trans as the time_trans of the previous timestep (not recommended).
        trans.time_trans = trans_sim.vec_trans.back().time_trans; 

        single_timestep single_h = single_next_h(trans, ckt, trans_sim, modmap);
        
        // The current time of the transient simulation. Use to compare with the breakpoints!
        double current_time = single_h.h + trans_sim.vec_trans.back().time_trans;

        /*  If the time difference between the previous simulation time point and the breakpoint is less than 10*hmin,
            it will be counted as 1 time point and no additional breakpoint simulation will be performed.*/
        if (trans_sim.breakpoints.empty() == false && trans_sim.breakpoints.front() - trans_sim.vec_trans.back().time_trans < 10 * trans_sim.trans_config.h_MIN)
        {
            trans_sim.breakpoints.pop_front();
            if (trans_sim.breakpoints.empty())
            {
                trans_sim.trans_end = true;
            }
        }

        if (trans_sim.breakpoints.empty() == false && current_time > trans_sim.breakpoints.front())
        {
            Transient breakpoints_trans;
            breakpoints_trans.mode = 1; // 0 to do OP analysis, 1 to do transient simulation
            breakpoints_trans.C_list = trans_sim.vec_trans.back().C_list;

            std::cout << "The time step is too large, back to breakpoint" << std::endl;
            breakpoints_trans.time_trans = trans_sim.breakpoints.front();
            breakpoints_trans.h = breakpoints_trans.time_trans - trans_sim.vec_trans.back().time_trans;

            breakpoints_trans.solution = NewtonRaphson_system(ckt, breakpoints_trans.h, 1, breakpoints_trans.time_trans, breakpoints_trans.C_list, trans_sim.vec_trans.back().solution, modmap);
            auto cur_vol = get_currents_voltages(breakpoints_trans.C_list, breakpoints_trans.h, breakpoints_trans.solution, trans_sim.vec_trans.back().solution);
            breakpoints_trans.Capacitance = get_capacitance(breakpoints_trans.C_list);

            breakpoints_trans.C_current = cur_vol.first;
            breakpoints_trans.C_voltage = cur_vol.second;
            breakpoints_trans.C_charge = breakpoints_trans.C_voltage % breakpoints_trans.Capacitance;
            breakpoints_trans.C_list_update();

            breakpoints_trans.trans_count = trans_sim.vec_trans.back().trans_count + 1;
            breakpoints_trans.next_h = breakpoints_trans.h;
            
            std::cout << "time step: " << breakpoints_trans.h << std::endl;
            std::cout << "time_trans: " << breakpoints_trans.time_trans << std::endl;

            trans_sim.breakpoints.pop_front();
            trans = std::move(breakpoints_trans);
        }

        else
        {
            if (trans_sim.breakpoints.empty() == false && current_time == trans_sim.breakpoints.front())
            {
                trans_sim.breakpoints.pop_front();
                if (trans_sim.breakpoints.empty())
                {
                    trans_sim.trans_end = true;
                }
            }

            // trans has already been updated in the multi_next_h function
            trans.h = single_h.h;
            trans.time_trans = trans_sim.vec_trans.back().time_trans + trans.h;
            trans.next_h = single_h.next_h;                          // Sometimes temp_h(next_h) sometimes is 2 * single_h.h!
            trans.solution = std::move(single_h.solution);
            trans.C_list = std::move(single_h.C_list);
            trans.Capacitance = std::move(single_h.Capacitance);
            trans.C_current = std::move(single_h.C_current);
            trans.C_voltage = std::move(single_h.C_voltage);
            trans.C_charge = std::move(single_h.C_charge);
            trans.C_list_update();
            trans.trans_count = trans_sim.vec_trans.back().trans_count + 1;
        }
    }

    return trans;
}

std::vector<Transient> Transient_ops(CKTcircuit &ckt, TransientSimulator &trans_sim, const Modelmap modmap)
{
    Transient trans_op;

    trans_op.LHS = ckt.cktdematrix->get_init_LHS();
    trans_op.RHS = ckt.cktdematrix->get_init_RHS();
    trans_op.C_list = ckt.C_list;               // Pass the capacitance list to the transient analysis
    trans_op.h = 0;

    // OPERATING POINT ANALYSIS SYSTEM
    trans_op.mode = 0;                          // 0 to do OP analysis, 1 to do transient simulation

    // OP analysis used as initial condition for next evaluation
    arma::vec solution = arma::zeros(ckt.cktdematrix->RHS.n_rows, ckt.cktdematrix->RHS.n_cols);

    // Benchmarking for OP analysis
    auto tstart_op = std::chrono::high_resolution_clock::now();

    solution = NewtonRaphson_system(ckt, 0, 0, 0, trans_op.C_list, solution, modmap);

    if (trans_op.C_list.size() > 0)
    {
        trans_op.Capacitance = get_capacitance(trans_op.C_list); // Get the capacitance matrix from c_list
        trans_op.C_current = arma::vec(trans_op.C_list.size(), arma::fill::zeros);
        trans_op.C_voltage = arma::vec(trans_op.C_list.size(), arma::fill::zeros);
        trans_op.C_charge = arma::vec(trans_op.C_list.size(), arma::fill::zeros);
    }
    auto tstop_op = std::chrono::high_resolution_clock::now();

    trans_op.solution = solution;
    trans_op.next_h = trans_sim.trans_config.init_h;
    history_trans_update(trans_op, trans_sim);

    arma::vec node_volt_print = solution.submat(0, 0, ckt.external_nodes - 1, 0);
    arma::vec current_print = solution.submat(solution.n_rows - ckt.no_of_V_sources, 0, solution.n_rows - 1, 0);
    arma::vec op_solution_print = arma::join_vert(node_volt_print, current_print);
    if (debugMode == false)
    {
        op_solution_print.print("The OP analysis of the circuit is: ");
    }

    ARMA_PRINT(solution, "The OP analysis (include internal nodes) of the circuit is: ");

    /*-----------------------------------------------------------*/
    //                Transient simulation start                   
    /*-----------------------------------------------------------*/
    auto tstart_trans = std::chrono::high_resolution_clock::now();
    std::cout << "transient simulation start" << std::endl;

    if (ckt.pulse_num != 0)
    {
        trans_sim.breakpoints = get_breakpoints(ckt, trans_sim); // Get the time of breakpoints for the transient simulation
        std::cout << "The breakpoints are: ";
        for (double num : trans_sim.breakpoints)
        {
            std::cout << num << " ";
        }
        std::cout << std::endl;
    }

    do
    {
        // std::cout << "-------------------------------------------------------" << std::endl;
        // std::cout << "Next step" << std::endl;

        // Fixed time step
        if (trans_sim.trans_config.timestep_control == false)
        {
           Transient trans = Fixed_TimeStep(ckt, trans_sim, modmap);
           history_trans_update(trans, trans_sim);
            
        }
        // time step control: varibale time step
        else{
            Transient trans = Varibale_TimeStep(ckt, trans_sim, modmap);
            history_trans_update(trans, trans_sim);
        }

    } while (trans_sim.vec_trans.back().time_trans < trans_sim.trans_config.t_end && trans_sim.trans_end == false);

    return trans_sim.vec_trans;
}

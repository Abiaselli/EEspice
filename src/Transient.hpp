#pragma once

#include "Transient_calcs.hpp"
#include "Transient_multi_solvers.hpp"

std::vector<Transient> Transient_ops(CKTcircuit &ckt, DenseMatrix &dematrix, Transient &trans_op)
{

    std::vector<Transient> vec_trans;

    trans_op.LHS = dematrix.get_init_LHS();
    trans_op.RHS = dematrix.get_init_RHS();
    trans_op.C_list = ckt.C_list; // Pass the capacitance list to the transient analysis
    trans_op.h = 0;

    /*----------------------------------------------fixed--------------------------------------------------------*/
    // Checking the LHS and RHS matrices
    // dematrix.LHS.print("Initial LHS matrix =");
    // dematrix.RHS.print("Initial RHS matrix =");
    ARMA_PRINT(dematrix.LHS, "Initial LHS matrix =");
    ARMA_PRINT(dematrix.RHS, "Initial RHS matrix =");
    /*----------------------------------------------fixed--------------------------------------------------------*/
    // OPERATING POINT ANALYSIS SYSTEM
    trans_op.mode = 0; // 0 to do OP analysis, 1 to do transient simulation

    // OP analysis used as initial condition for next evaluation
    arma::vec solution = arma::zeros(dematrix.RHS.n_rows, dematrix.RHS.n_cols);

    // Benchmarking for OP analysis
    auto tstart_op = std::chrono::high_resolution_clock::now();

    solution = NewtonRaphson_system(ckt, solution, 0, 0, 0, trans_op.C_list, vec_trans.back().solution);

    if (trans_op.C_list.size() > 0)
    {
        auto cur_vol_op = get_currents_voltages(trans_op.C_list, trans_op.h, solution, arma::zeros(dematrix.RHS.n_rows, dematrix.RHS.n_cols));
        trans_op.Capacitance = get_capacitance(trans_op.C_list); // Get the capacitance matrix from c_list
        trans_op.C_current = cur_vol_op.first;
        // trans_op.C_current.print("The current matrix at op is: ");
        trans_op.C_voltage = cur_vol_op.second;
        trans_op.C_charge = trans_op.C_voltage % trans_op.Capacitance;
    }
    auto tstop_op = std::chrono::high_resolution_clock::now();

    trans_op.solution = solution;
    history_trans_update(trans_op, vec_trans);

    arma::vec node_volt_print = solution.submat(0, 0, ckt.external_nodes - 1, 0);
    arma::vec current_print = solution.submat(solution.n_rows - ckt.no_of_V_sources, 0, solution.n_rows - 1, 0);
    arma::vec op_solution_print = arma::join_vert(node_volt_print, current_print);
    if (debugMode == false)
    {
        op_solution_print.print("The OP analysis of the circuit is: ");
    }

    ARMA_PRINT(solution, "The OP analysis (include internal nodes) of the circuit is: ");
    /*----------------------------------------------fixed--------------------------------------------------------*/

    /*--------------------------------------------can be changed-------------------------------------------------*/
    // ADDING TRANSIENT SIMULATION LOOP (includes V_pulse or any time dependent sources)

    auto tstart_trans = std::chrono::high_resolution_clock::now();
    std::cout << "transient simulation start" << std::endl;

    std::deque<double> breakpoints;
    bool trans_end = false;

    if (ckt.pulse_num != 0)
    {
        breakpoints = get_breakpoints(ckt, trans_op); // Get the time of breakpoints for the transient simulation
        std::cout << "The breakpoints are: ";
        for (double num : breakpoints)
        {
            std::cout << num << " ";
        }
        std::cout << std::endl;
    }

    do
    {
        // std::cout << "-------------------------------------------------------" << std::endl;
        // std::cout << "Next step" << std::endl;

        Transient trans;
        trans.t_end = trans_op.t_end;
        trans.h_MAX = trans_op.h_MAX;
        trans.h_MIN = trans_op.h_MIN;
        trans.mode = 1; // 0 to do OP analysis, 1 to do transient simulation
        // trans.C_list = trans_op.C_list;

        // Fixed time step
        if (ckt.pulse_num == 0 && trans_op.C_list.size() == 0)
        {
            trans.h = trans_op.init_h; // Initial time step
            trans.time_trans = vec_trans.back().time_trans;
            trans.C_list = vec_trans.back().C_list;

            std::cout << "Fixed time step" << std::endl;
            solution = NewtonRaphson_system(ckt, vec_trans.back().solution, trans.h, 1, trans.time_trans, trans.C_list, vec_trans.back().solution);
            // solution.print("The transient analysis of the circuit is: ");
            ARMA_PRINT(solution, "The transient analysis of the circuit is: ");

            if (trans.C_list.size() > 0)
            {

                auto cur_vol = get_currents_voltages(trans.C_list, trans.h, solution, vec_trans.back().solution);
                trans.Capacitance = get_capacitance(trans.C_list);
                trans.C_current = cur_vol.first;
                trans.C_voltage = cur_vol.second;
                trans.C_charge = trans.C_voltage % trans.Capacitance;
                trans.C_list_update();

                /*--------------------------------------*/
                double current = arma::abs(trans.C_current).max();
                double pre_current = arma::abs(vec_trans.back().C_current).max();
                double charge = arma::abs(trans.C_charge).max();
                double pre_charge = arma::abs(vec_trans.back().C_charge).max();
                double LTE_B_C = RELTOL * std::max(current, pre_current) + ABSTOL;
                double LTE_B_Q = RELTOL * std::max({charge, pre_charge, CHGTOL}) / trans.h;
                double LTE_bound = std::max(LTE_B_C, LTE_B_Q);

                double LTE = trans.h / (2 * trans.C_list.at(0).value) * std::abs((current - pre_current));

                // if(LTE > 7 * LTE_bound){
                //     std::cout << "LTE > LTE_bound" << std::endl;
                // }
            }

            trans.time_trans = vec_trans.back().time_trans + trans.h;
            trans.solution = solution;
            // trans.LHS = matrixs.first;
            // trans.RHS = matrixs.second;
            trans.trans_count = vec_trans.back().trans_count + 1;
            // trans.C_current.print("The current matrix is: ");
            // std::cout << "time trans: " << trans.time_trans << std::endl;
            // std::cout << "time step: " << trans.h << std::endl;

            history_trans_update(trans, vec_trans);

            continue;
        }

        // time step control: varibale time step

        if (vec_trans.size() == 1)
        {
            // The first step of the transient simulation
            trans.C_list = trans_op.C_list;
            trans.h = trans_op.init_h; // Initial time step
            trans.time_trans = trans.t_start + trans.h;

            solution = NewtonRaphson_system(ckt, trans_op.solution, trans.h, 1, trans.time_trans, trans.C_list, vec_trans.back().solution);
            // solution.print("The transient analysis of the circuit is: ");
            ARMA_PRINT(solution, "The transient analysis of the circuit is: ");

            if (trans.C_list.size() > 0)
            {

                std::pair<arma::vec, arma::vec> cur_vol = get_currents_voltages(trans.C_list, trans.h, solution, trans_op.solution);
                trans.Capacitance = get_capacitance(trans.C_list);

                trans.C_current = cur_vol.first;
                trans.C_voltage = cur_vol.second;
                trans.C_charge = trans.C_voltage % trans.Capacitance;
                trans.C_list_update();
            }

            trans.solution = solution;
            // trans.LHS = matrixs.first;
            // trans.RHS = matrixs.second;
            trans.trans_count += 1;
            // trans.C_current.print("The current matrix is: ");

            history_trans_update(trans, vec_trans);
        }

        else
        {

            trans.C_list = vec_trans.back().C_list;

            trans.time_trans = vec_trans.back().time_trans;

            solution = multi_next_h(trans, ckt, vec_trans);

            trans.time_trans = trans.time_trans + trans.h;

            /* If the time difference between the previous simulation time point and the breakpoint is less than 10*hmin,
               it will be counted as 1 time point and no additional breakpoint simulation will be performed.*/
            if (breakpoints.empty() == false && breakpoints.front() - vec_trans.back().time_trans < 10 * trans.h_MIN)
            {

                breakpoints.pop_front();
                if (breakpoints.empty())
                {
                    trans_end = true;
                }
            }

            if (breakpoints.empty() == false && trans.time_trans > breakpoints.front())
            {

                Transient breakpoints_trans;
                breakpoints_trans.t_end = trans_op.t_end;
                breakpoints_trans.h_MAX = trans_op.h_MAX;
                breakpoints_trans.h_MIN = trans_op.h_MIN;
                breakpoints_trans.mode = 1; // 0 to do OP analysis, 1 to do transient simulation
                breakpoints_trans.C_list = vec_trans.back().C_list;

                std::cout << "The time step is too large, back to breakpoint" << std::endl;
                breakpoints_trans.time_trans = breakpoints.front();
                breakpoints_trans.h = breakpoints_trans.time_trans - vec_trans.back().time_trans;
                // std::cout << "breakpoints_trans.time_trans: " << breakpoints_trans.time_trans << std::endl;
                // std::cout << "vec_trans.back().time_trans: " << vec_trans.back().time_trans << std::endl;

                DEBUG_PRINT("breakpoints_trans.time_trans: " << breakpoints_trans.time_trans);
                DEBUG_PRINT("vec_trans.back().time_trans: " << vec_trans.back().time_trans);

                solution = NewtonRaphson_system(ckt, vec_trans.back().solution, breakpoints_trans.h, 1, breakpoints_trans.time_trans, breakpoints_trans.C_list, vec_trans.back().solution);
                auto cur_vol = get_currents_voltages(breakpoints_trans.C_list, breakpoints_trans.h, solution, vec_trans.back().solution);
                breakpoints_trans.Capacitance = get_capacitance(breakpoints_trans.C_list);

                breakpoints_trans.solution = solution;
                breakpoints_trans.C_current = cur_vol.first;
                breakpoints_trans.C_voltage = cur_vol.second;
                breakpoints_trans.C_charge = breakpoints_trans.C_voltage % breakpoints_trans.Capacitance;
                breakpoints_trans.C_list_update();

                // breakpoints_trans.LHS = matrixs.first;
                // breakpoints_trans.RHS = matrixs.second;
                breakpoints_trans.trans_count = vec_trans.back().trans_count + 1;

                history_trans_update(breakpoints_trans, vec_trans);

                // solution.print("The transient analysis of the circuit is: ");
                ARMA_PRINT(solution, "The transient analysis of the circuit is: ");
                std::cout << "time step: " << breakpoints_trans.h << std::endl;
                std::cout << "time_trans: " << breakpoints_trans.time_trans << std::endl;

                breakpoints.pop_front();
                // continue;
            }

            else
            {

                if (breakpoints.empty() == false && trans.time_trans == breakpoints.front())
                {
                    breakpoints.pop_front();
                    if (breakpoints.empty())
                    {
                        trans_end = true;
                    }
                }
                ARMA_PRINT(solution, "The transient analysis of the circuit is: ");

                // trans has already been updated in the multi_next_h function
                trans.C_list_update();
                trans.trans_count = vec_trans.back().trans_count + 1;
                history_trans_update(trans, vec_trans);
            }
        }

    } while (vec_trans.back().time_trans < vec_trans.back().t_end && trans_end == false);

    return vec_trans;
}

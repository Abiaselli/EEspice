#pragma once

#include <iostream>
#include <armadillo>
#include <fstream>
#include <tuple>
#include <cmath>
#include <chrono>
#include <string>
#include <variant>
#include <vector>
#include <algorithm>
#include <deque>
#include <iomanip>
#include <typeinfo>
#include "sim_variables.hpp"
#include "BS_thread_pool/BS_thread_pool.hpp"
#include "BS_thread_pool/BS_thread_pool_utils.hpp"
#include "circuit_parser.hpp"
#include "XB_timer.hpp"
// #include <mutex>
// #include <bandicoot>

#include "Transient_calcs.hpp"

struct multi_timestep
{
    arma::vec C_current_mid;
    arma::vec C_voltage_mid;
    arma::vec C_charge_mid;
    arma::vec Capacitance_mid;

    arma::vec C_current_up;
    arma::vec C_voltage_up;
    arma::vec C_charge_up;
    arma::vec Capacitance_up;

    arma::vec C_current_down;
    arma::vec C_voltage_down;
    arma::vec C_charge_down;
    arma::vec Capacitance_down;

    double h_mid{};
    double h_up{};
    double h_down{};

    double t_mid{};
    double t_up{};
    double t_down{};

    arma::vec RHS_mid;
    arma::vec RHS_up;
    arma::vec RHS_down;

    arma::mat LHS_mid;
    arma::mat LHS_up;
    arma::mat LHS_down;

    arma::vec solution_mid;
    arma::vec solution_up;
    arma::vec solution_down;

    bool TMAX_reach;
    bool TMIN_reach;

    std::vector<Capacitor> C_list_mid;
    std::vector<Capacitor> C_list_up;
    std::vector<Capacitor> C_list_down;
};

// New timestep options
void timestep_options(double &temp_h, double &next_h_up, double &next_h_down, const Transient &trans, bool &TMAX_reach, bool &TMIN_reach)
{
    if (temp_h >= trans.config->h_MAX)
    {
        temp_h = trans.config->h_MAX;
        TMAX_reach = true;
    }
    else if (temp_h <= trans.config->h_MIN)
    {
        temp_h = trans.config->h_MIN;
        TMIN_reach = true;
    }

    next_h_up = temp_h * 2;
    next_h_down = temp_h / 2;

    if (next_h_up > trans.config->h_MAX)
    {
        next_h_up = trans.config->h_MAX;
        // TMAX_reach = true;
    }
    else if (next_h_up < trans.config->h_MIN)
    {
        next_h_down = trans.config->h_MIN;
        // TMAX_reach = true;
    }

    if (next_h_down > trans.config->h_MAX)
    {
        next_h_up = trans.config->h_MAX;
        // TMAX_reach = true;
    }
    else if (next_h_down < trans.config->h_MIN)
    {
        next_h_down = trans.config->h_MIN;
        // TMIN_reach = true;
    }
}

int h_reach_new(const Truncation_error LTE, int &last_decision, const multi_timestep &multi_h)
{

    /* index_h:
        0 = temp_h / 4
        1 = down
        2 = mid
        3 = up
    */
    int reach{};

    if (multi_h.TMAX_reach)
    {
        bool mid = multi_h.h_mid <= TRTOL * LTE.LTE_BE_mid_h.min();

        if (mid == false)
        {
            last_decision = 1;
            reach = 0; // temp_h / 4
            return reach;
        }
        else
        {
            reach = 2; // mid
            return reach;
        }
    }
    else if (multi_h.TMIN_reach)
    {
        bool mid = multi_h.h_mid <= TRTOL * LTE.LTE_BE_mid_h.min();
        bool up = multi_h.h_up <= TRTOL * LTE.LTE_BE_up_h.min();

        if (mid == false)
        {
            std::cerr << "The time step is too small (Reach minimum!)" << std::endl;
            exit(1);
        }
        else if (up == true)
        {
            reach = 3; // up
            return reach;
        }
        else
        {
            reach = 2; // mid
            return reach;
        }
    }
    else
    {
        bool mid = multi_h.h_mid <= TRTOL * LTE.LTE_BE_mid_h.min();
        bool up = multi_h.h_up <= TRTOL * LTE.LTE_BE_up_h.min();
        bool down = multi_h.h_down <= TRTOL * LTE.LTE_BE_down_h.min();

        if (down == false)
        {
            last_decision = 1;
            reach = 0; // temp_h / 4
            return reach;
        }

        else if (mid == true && up == true && down == true)
        {
            reach = 3; // up
            return reach;
        }

        else if (mid == false && up == false && down == false)
        {
            reach = 0; // temp_h / 4
            return reach;
        }

        else if (mid == true && up == false && down == true)
        {

            reach = 2; // mid
            return reach;
        }

        else if (mid == false && up == false && down == true)
        {
            reach = 1; // down
            return reach;
        }

        std::cerr << "h_reach_new function error" << std::endl;
        exit(1);
    }
}

multi_timestep multi_solution_solver(const double &h, const Transient &trans, const CKTcircuit &ckt, std::vector<Transient> &vec_trans)
{

    multi_timestep multi_h;

    multi_h.TMAX_reach = false;
    multi_h.TMIN_reach = false;

    multi_h.h_mid = h;

    timestep_options(multi_h.h_mid, multi_h.h_up, multi_h.h_down, trans, multi_h.TMAX_reach, multi_h.TMIN_reach); // calculate the next time step, up mid down

    // Generic lambda function for solving
    auto solve = [&trans, &ckt](const double h, double &t, arma::vec &solution, arma::vec &C_current, arma::vec &C_voltage,
                                arma::vec &C_charge, arma::vec &capacitance, std::vector<Capacitor> &C_list, std::vector<Transient> &vec_trans)
    {
        t = trans.time_trans + h;
        std::vector<Capacitor> C_list_copy = trans.C_list;
        arma::vec local_solution = NewtonRaphson_system(ckt, h, 1, t, C_list_copy, vec_trans.back().solution);

        solution = std::move(local_solution);
        std::pair<arma::vec, arma::vec> currents_voltages = get_currents_voltages(C_list_copy, h, solution, vec_trans.back().solution);
        C_current = std::move(currents_voltages.first);
        C_voltage = std::move(currents_voltages.second);
        capacitance = get_capacitance(C_list_copy);
        C_charge = capacitance % C_voltage; // element-wise multiplication of two objects (Schur product)
        C_list = std::move(C_list_copy);
    };

    if (multi_h.TMAX_reach)
    {
        // Only Mid time step
        solve(multi_h.h_mid, multi_h.t_mid, multi_h.solution_mid, multi_h.C_current_mid, multi_h.C_voltage_mid, multi_h.C_charge_mid, multi_h.Capacitance_mid, multi_h.C_list_mid, vec_trans);
    }
    else if (multi_h.TMIN_reach)
    {
        // Mid time step
        pool.detach_task(
            [&multi_h, &solve, &vec_trans]
            {
                // auto t1 = std::chrono::high_resolution_clock::now();
                solve(multi_h.h_mid, multi_h.t_mid, multi_h.solution_mid, multi_h.C_current_mid, multi_h.C_voltage_mid, multi_h.C_charge_mid, multi_h.Capacitance_mid, multi_h.C_list_mid, vec_trans);
                // auto t2 = std::chrono::high_resolution_clock::now();
                // begin_mid.push_back(t1);
                // end_mid.push_back(t2);
            });

        // Up time step
        pool.detach_task(
            [&multi_h, &solve, &vec_trans]
            {
                // auto t1 = std::chrono::high_resolution_clock::now();
                solve(multi_h.h_up, multi_h.t_up, multi_h.solution_up, multi_h.C_current_up, multi_h.C_voltage_up, multi_h.C_charge_up, multi_h.Capacitance_up, multi_h.C_list_up, vec_trans);
                // auto t2 = std::chrono::high_resolution_clock::now();
                // begin_up.push_back(t1);
                // end_up.push_back(t2);
            });
    }
    else
    {
        // Mid time step
        pool.detach_task(
            [&multi_h, &solve, &vec_trans]
            {
                // auto t1 = std::chrono::high_resolution_clock::now();
                solve(multi_h.h_mid, multi_h.t_mid, multi_h.solution_mid, multi_h.C_current_mid, multi_h.C_voltage_mid, multi_h.C_charge_mid, multi_h.Capacitance_mid, multi_h.C_list_mid, vec_trans);
                // auto t2 = std::chrono::high_resolution_clock::now();
                // begin_mid.push_back(t1);
                // end_mid.push_back(t2);
            });

        // Up time step
        pool.detach_task(
            [&multi_h, &solve, &vec_trans]
            {
                // auto t1 = std::chrono::high_resolution_clock::now();
                solve(multi_h.h_up, multi_h.t_up, multi_h.solution_up, multi_h.C_current_up, multi_h.C_voltage_up, multi_h.C_charge_up, multi_h.Capacitance_up, multi_h.C_list_up, vec_trans);
                // auto t2 = std::chrono::high_resolution_clock::now();
                // begin_up.push_back(t1);
                // end_up.push_back(t2);
            });

        // Down time step
        pool.detach_task(
            [&multi_h, &solve, &vec_trans]
            {
                // auto t1 = std::chrono::high_resolution_clock::now();
                solve(multi_h.h_down, multi_h.t_down, multi_h.solution_down, multi_h.C_current_down, multi_h.C_voltage_down, multi_h.C_charge_down, multi_h.Capacitance_down, multi_h.C_list_down, vec_trans);
                // auto t2 = std::chrono::high_resolution_clock::now();
                // begin_down.push_back(t1);
                // end_down.push_back(t2);
            });
    }

    pool.wait();

    return multi_h;
}

arma::vec multi_next_h(Transient &trans, const CKTcircuit &ckt, std::vector<Transient> &vec_trans)
{

    // std::cout << "Enter multi_next_h" << std::endl;
    DEBUG_PRINT("Enter multi_next_h");

    arma::vec solution;

    double temp_h = vec_trans.back().h;
    double lats_step = vec_trans.back().h;
    arma::vec pre_current = vec_trans.back().C_current;
    arma::vec pre_charge = vec_trans.back().C_charge;

    int index_h{};
    int pre_decision{}; // 1 = h / 4, 0 = nan
    multi_timestep multi_h;
    multi_timestep pre_multi_h;
    Truncation_error LTE;

    do
    {
        timer.start();

        multi_h = multi_solution_solver(temp_h, trans, ckt, vec_trans);

        timer.stop();
        timer.total();

        if (multi_h.TMAX_reach)
        {
            total_timepoint += 1;
            arma::vec mid_max_current = arma::max(arma::abs(multi_h.C_current_mid), arma::abs(pre_current));

            arma::vec mid_max_charge = arma::max(arma::abs(multi_h.C_charge_mid), arma::abs(pre_charge));

            // LTE current = ABSTOL + RELTOL * max(|ik+1|,|ik|) in SPICE book
            arma::vec mid_LTE_current = ABSTOL + RELTOL * mid_max_current;

            // LTE charge = RELTOL * max(|qk+1|,|qk|, chgtol) / h(k) in SPICE book
            arma::vec CHGTOL_MAT = arma::vec(mid_max_charge.n_rows, mid_max_charge.n_cols).fill(CHGTOL);
            arma::vec mid_LTE_charge = RELTOL * arma::max(mid_max_charge, CHGTOL_MAT) / multi_h.h_mid;

            // LTE bound = max(LTE_current, LTE_charge)
            LTE.LTE_bound_mid = arma::max(mid_LTE_current, mid_LTE_charge);

            arma::vec BEcur_diff_mid = multi_h.C_current_mid - pre_current;

            // tol_h = 2C / |i(k+1) - i(k)| * LTE_bound
            LTE.LTE_BE_mid_h = ((2 * multi_h.Capacitance_mid) / arma::abs(BEcur_diff_mid)) % LTE.LTE_bound_mid;
        }
        else if (multi_h.TMIN_reach)
        {
            total_timepoint += 2;
            arma::vec mid_max_current = arma::max(arma::abs(multi_h.C_current_mid), arma::abs(pre_current));
            arma::vec up_max_current = arma::max(arma::abs(multi_h.C_current_up), arma::abs(pre_current));

            arma::vec mid_max_charge = arma::max(arma::abs(multi_h.C_charge_mid), arma::abs(pre_charge));
            arma::vec up_max_charge = arma::max(arma::abs(multi_h.C_charge_up), arma::abs(pre_charge));

            // LTE current = ABSTOL + RELTOL * max(|ik+1|,|ik|) in SPICE book
            arma::vec mid_LTE_current = ABSTOL + RELTOL * mid_max_current;
            arma::vec up_LTE_current = ABSTOL + RELTOL * up_max_current;

            // LTE charge = RELTOL * max(|qk+1|,|qk|, chgtol) / h(k) in SPICE book
            arma::vec CHGTOL_MAT = arma::vec(mid_max_charge.n_rows, mid_max_charge.n_cols).fill(CHGTOL);
            arma::vec mid_LTE_charge = RELTOL * arma::max(mid_max_charge, CHGTOL_MAT) / multi_h.h_mid;
            arma::vec up_LTE_charge = RELTOL * arma::max(up_max_charge, CHGTOL_MAT) / multi_h.h_up;

            // LTE bound = max(LTE_current, LTE_charge)
            LTE.LTE_bound_mid = arma::max(mid_LTE_current, mid_LTE_charge);
            LTE.LTE_bound_up = arma::max(up_LTE_current, up_LTE_charge);

            arma::vec BEcur_diff_mid = multi_h.C_current_mid - pre_current;
            arma::vec BEcur_diff_up = multi_h.C_current_up - pre_current;

            // tol_h = 2C / |i(k+1) - i(k)| * LTE_bound
            LTE.LTE_BE_mid_h = ((2 * multi_h.Capacitance_mid) / arma::abs(BEcur_diff_mid)) % LTE.LTE_bound_mid;
            LTE.LTE_BE_up_h = ((2 * multi_h.Capacitance_up) / arma::abs(BEcur_diff_up)) % LTE.LTE_bound_up;
        }
        else
        {
            total_timepoint += 3;
            arma::vec mid_max_current = arma::max(arma::abs(multi_h.C_current_mid), arma::abs(pre_current));
            arma::vec up_max_current = arma::max(arma::abs(multi_h.C_current_up), arma::abs(pre_current));
            arma::vec down_max_current = arma::max(arma::abs(multi_h.C_current_down), arma::abs(pre_current));

            arma::vec mid_max_charge = arma::max(arma::abs(multi_h.C_charge_mid), arma::abs(pre_charge));
            arma::vec up_max_charge = arma::max(arma::abs(multi_h.C_charge_up), arma::abs(pre_charge));
            arma::vec down_max_charge = arma::max(arma::abs(multi_h.C_charge_down), arma::abs(pre_charge));

            // LTE current = ABSTOL + RELTOL * max(|ik+1|,|ik|) in SPICE book
            arma::vec mid_LTE_current = ABSTOL + RELTOL * mid_max_current;
            arma::vec up_LTE_current = ABSTOL + RELTOL * up_max_current;
            arma::vec down_LTE_current = ABSTOL + RELTOL * down_max_current;

            // LTE charge = RELTOL * max(|qk+1|,|qk|, chgtol) / h(k) in SPICE book
            arma::vec CHGTOL_MAT = arma::vec(mid_max_charge.n_rows, mid_max_charge.n_cols).fill(CHGTOL);
            arma::vec mid_LTE_charge = RELTOL * arma::max(mid_max_charge, CHGTOL_MAT) / multi_h.h_mid;
            arma::vec up_LTE_charge = RELTOL * arma::max(up_max_charge, CHGTOL_MAT) / multi_h.h_up;
            arma::vec down_LTE_charge = RELTOL * arma::max(down_max_charge, CHGTOL_MAT) / multi_h.h_down;

            // LTE bound = max(LTE_current, LTE_charge)
            LTE.LTE_bound_mid = arma::max(mid_LTE_current, mid_LTE_charge);
            LTE.LTE_bound_up = arma::max(up_LTE_current, up_LTE_charge);
            LTE.LTE_bound_down = arma::max(down_LTE_current, down_LTE_charge);

            arma::vec BEcur_diff_mid = multi_h.C_current_mid - pre_current;
            arma::vec BEcur_diff_up = multi_h.C_current_up - pre_current;
            arma::vec BEcur_diff_down = multi_h.C_current_down - pre_current;

            // tol_h = 2C / |i(k+1) - i(k)| * LTE_bound
            LTE.LTE_BE_mid_h = ((2 * multi_h.Capacitance_mid) / arma::abs(BEcur_diff_mid)) % LTE.LTE_bound_mid;
            LTE.LTE_BE_up_h = ((2 * multi_h.Capacitance_up) / arma::abs(BEcur_diff_up)) % LTE.LTE_bound_up;
            LTE.LTE_BE_down_h = ((2 * multi_h.Capacitance_down) / arma::abs(BEcur_diff_down)) % LTE.LTE_bound_down;
        }

        index_h = h_reach_new(LTE, pre_decision, multi_h);

        /* index_h:
            0 = temp_h / 4
            1 = down
            2 = mid
            3 = up
        */
        switch (index_h)
        {
        case 0:
            temp_h = temp_h / 4.0;
            pre_multi_h = multi_h;
            break;
        case 1:
            trans.h = multi_h.h_down;
            trans.LHS = multi_h.LHS_down;
            trans.RHS = multi_h.RHS_down;
            trans.solution = multi_h.solution_down;
            trans.C_list = multi_h.C_list_down;
            trans.Capacitance = multi_h.Capacitance_down;
            solution = multi_h.solution_down;

            trans.C_current = multi_h.C_current_down;
            trans.C_voltage = multi_h.C_voltage_down;
            trans.C_charge = multi_h.C_charge_down;
            break;
        case 2:
            trans.h = multi_h.h_mid;
            trans.LHS = multi_h.LHS_mid;
            trans.RHS = multi_h.RHS_mid;
            trans.solution = multi_h.solution_mid;
            trans.C_list = multi_h.C_list_mid;
            trans.Capacitance = multi_h.Capacitance_mid;
            solution = multi_h.solution_mid;

            trans.C_current = multi_h.C_current_mid;
            trans.C_voltage = multi_h.C_voltage_mid;
            trans.C_charge = multi_h.C_charge_mid;
            break;
        case 3:
            trans.h = multi_h.h_up;
            trans.LHS = multi_h.LHS_up;
            trans.RHS = multi_h.RHS_up;
            trans.solution = multi_h.solution_up;
            trans.C_list = multi_h.C_list_up;
            trans.Capacitance = multi_h.Capacitance_up;
            solution = multi_h.solution_up;

            trans.C_current = multi_h.C_current_up;
            trans.C_voltage = multi_h.C_voltage_up;
            trans.C_charge = multi_h.C_charge_up;
            break;

        default:
            std::cerr << "Error in multi_next_h function" << std::endl;
            exit(1);
        }

    } while (index_h == 0);

    return solution;
}

// multi_timestep multi_solution_solver(const double &h, const Transient &trans, const CKTcircuit &ckt){

//     multi_timestep multi_h;

//     multi_h.TMAX_reach = false;
//     multi_h.TMIN_reach = false;
//     multi_h.ITE_reach_mid = false;
//     multi_h.ITE_reach_up = false;
//     multi_h.ITE_reach_down = false;

//     multi_h.h_mid = h;

//     timestep_options(multi_h.h_mid, multi_h.h_up, multi_h.h_down, trans, multi_h.TMAX_reach, multi_h.TMIN_reach);  // calculate the next time step, up mid down

//     multi_h.t_mid = trans.time_trans + multi_h.h_mid;
//     multi_h.t_up = trans.time_trans + multi_h.h_up;
//     multi_h.t_down = trans.time_trans + multi_h.h_down;

//     // Double check the time steps
//     // if(multi_h.t_mid == multi_h.t_up){
//     //     if(multi_h.TMAX_reach == false){
//     //         std::cerr << "The TMAX_reach is false in multi_solution_solver function" << std::endl;
//     //         exit(1);
//     //     }
//     // }

//     // if(multi_h.t_mid == multi_h.t_down){
//     //     if(multi_h.TMIN_reach == false){
//     //         std::cerr << "The TMIN_reach is false in multi_solution_solver function" << std::endl;
//     //         exit(1);
//     //     }
//     // }

//     // std::cout << "The time steps are:" << multi_h.h_mid << " " << multi_h.h_up << " " << multi_h.h_down << std::endl;
//     DEBUG_PRINT("The time steps are: " << multi_h.h_mid << " " << multi_h.h_up << " " << multi_h.h_down);

//     // Mid time step
//     pool.detach_task(
//         [&trans, &ckt, &multi_h]
//         {
//             // begin_mid.push_back(std::chrono::high_resolution_clock::now());
//             std::vector<Capacitor> C_list_mid = trans.C_list;
//             arma::mat mid_solution = NewtonRaphson_system(ckt, vec_trans.back().solution, multi_h.h_mid, 1, multi_h.t_mid, C_list_mid);
//             ARMA_PRINT(mid_solution, "mid_solution =");

//                 multi_h.solution_mid = mid_solution;
//                 std::pair<arma::mat, arma::mat> mid_currents_voltages = get_currents_voltages(C_list_mid, multi_h.h_mid, mid_solution, vec_trans.back().solution);
//                 multi_h.C_current_mid = mid_currents_voltages.first;
//                 multi_h.C_voltage_mid = mid_currents_voltages.second;
//                 arma::mat mid_Capacitance = get_capacitance(C_list_mid);       // Update the capacitance values because the MOSFETs
//                 multi_h.C_charge_mid = mid_Capacitance % multi_h.C_voltage_mid;  // element-wise multiplication of two objects (Schur product)
//                 multi_h.Capacitance_mid = mid_Capacitance;
//                 multi_h.C_list_mid = C_list_mid;

//             // end_mid.push_back(std::chrono::high_resolution_clock::now());

//             // sync_out.print("Time mid_solution_solver = ", tmr.current_ms(),"\n");
//         }
//     );

//     // Up time step
//     // if(multi_h.TMAX_reach == false){
//         pool.detach_task(
//             [&trans, &ckt, &multi_h]
//             {
//                 // begin_up.push_back(std::chrono::high_resolution_clock::now());
//                 std::vector<Capacitor> C_list_up = trans.C_list;
//                 arma::mat up_solution = NewtonRaphson_system(ckt, vec_trans.back().solution, multi_h.h_up, 1, multi_h.t_up, C_list_up);

//                 ARMA_PRINT(up_solution, "up_solution =");

//                     multi_h.solution_up = up_solution;
//                     std::pair<arma::mat, arma::mat> up_currents_voltages = get_currents_voltages(C_list_up, multi_h.h_up, up_solution, vec_trans.back().solution);
//                     multi_h.C_current_up = up_currents_voltages.first;
//                     multi_h.C_voltage_up = up_currents_voltages.second;
//                     arma::mat up_Capacitance = get_capacitance(C_list_up);       // Update the capacitance values because the MOSFETs
//                     multi_h.C_charge_up = up_Capacitance % multi_h.C_voltage_up;  // element-wise multiplication of two objects (Schur product)
//                     multi_h.Capacitance_up = up_Capacitance;
//                     multi_h.C_list_up = C_list_up;

//                 // end_up.push_back(std::chrono::high_resolution_clock::now());
//                 // sync_out.print("Time up_solution_solver = ", tmr.current_ms(),"\n");
//             }
//         );
//     // }

//     // Down time step
//     // if(multi_h.TMIN_reach == false){
//         pool.detach_task(
//             [&trans, &ckt, &multi_h]
//             {
//                 // begin_down.push_back(std::chrono::high_resolution_clock::now());
//                 std::vector<Capacitor> C_list_down = trans.C_list;
//                 arma::mat down_solution = NewtonRaphson_system(ckt, vec_trans.back().solution, multi_h.h_down, 1, multi_h.t_down, C_list_down);

//                 ARMA_PRINT(down_solution, "down_solution =");

//                     multi_h.solution_down = down_solution;
//                     std::pair<arma::mat, arma::mat> down_currents_voltages = get_currents_voltages(C_list_down, multi_h.h_down, down_solution, vec_trans.back().solution);
//                     multi_h.C_current_down = down_currents_voltages.first;
//                     multi_h.C_voltage_down = down_currents_voltages.second;
//                     arma::mat down_Capacitance = get_capacitance(C_list_down);       // Update the capacitance values because the MOSFETs
//                     multi_h.C_charge_down = down_Capacitance % multi_h.C_voltage_down;  // element-wise multiplication of two objects (Schur product)
//                     multi_h.Capacitance_down = down_Capacitance;
//                     multi_h.C_list_down = C_list_down;

//                 // end_down.push_back(std::chrono::high_resolution_clock::now());
//             }
//         );
//     // }

//     pool.wait();

//     if(multi_h.h_mid == trans.h_MIN && multi_h.ITE_reach_mid && multi_h.ITE_reach_up && multi_h.ITE_reach_down){
//         std::cerr << "The time step is too small (Reach minimum!)" << std::endl;
//         exit(1);
//     }

//     // std::vector<Capacitor> C_list_mid = trans.C_list;
//     // arma::mat mid_solution = NewtonRaphson_system(ckt, vec_trans.back().solution, multi_h.h_mid, 1, multi_h.t_mid, C_list_mid);
//     // // mid_solution.print("mid_solution =");
//     // ARMA_PRINT(mid_solution, "mid_solution =");
//     // multi_h.solution_mid = mid_solution;
//     // auto mid_currents_voltages = get_currents_voltages(C_list_mid, multi_h.h_mid, mid_solution, vec_trans.back().solution);
//     // multi_h.C_current_mid = mid_currents_voltages.first;
//     // multi_h.C_voltage_mid = mid_currents_voltages.second;
//     // arma::mat mid_Capacitance = get_capacitance(C_list_mid);       // Update the capacitance values because the MOSFETs
//     // multi_h.C_charge_mid = mid_Capacitance % multi_h.C_voltage_mid;  // element-wise multiplication of two objects (Schur product)
//     // // multi_h.LHS_mid = mid_matrixs.first;
//     // // multi_h.RHS_mid = mid_matrixs.second;
//     // multi_h.Capacitance_mid = mid_Capacitance;
//     // multi_h.C_list_mid = C_list_mid;

//     // std::vector<Capacitor> C_list_up = trans.C_list;
//     // arma::mat up_solution = NewtonRaphson_system(ckt, vec_trans.back().solution, multi_h.h_up, 1, multi_h.t_up, C_list_up);
//     // // up_solution.print("up_solution =");
//     // ARMA_PRINT(up_solution, "up_solution =");
//     // multi_h.solution_up = up_solution;
//     // auto up_currents_voltages = get_currents_voltages(C_list_up, multi_h.h_up, up_solution, vec_trans.back().solution);
//     // multi_h.C_current_up = up_currents_voltages.first;
//     // multi_h.C_voltage_up = up_currents_voltages.second;
//     // arma::mat up_Capacitance = get_capacitance(C_list_up);       // Update the capacitance values because the MOSFETs
//     // multi_h.C_charge_up = up_Capacitance % multi_h.C_voltage_up;  // element-wise multiplication of two objects (Schur product)
//     // // multi_h.LHS_up = up_matrixs.first;
//     // // multi_h.RHS_up = up_matrixs.second;
//     // multi_h.Capacitance_up = up_Capacitance;
//     // multi_h.C_list_up = C_list_up;

//     // std::vector<Capacitor> C_list_down = trans.C_list;
//     // arma::mat down_solution = NewtonRaphson_system(ckt, vec_trans.back().solution, multi_h.h_down, 1, multi_h.t_down, C_list_down);
//     // // down_solution.print("down_solution =");
//     // ARMA_PRINT(down_solution, "down_solution =");
//     // multi_h.solution_down = down_solution;
//     // auto down_currents_voltages = get_currents_voltages(C_list_down, multi_h.h_down, down_solution, vec_trans.back().solution);
//     // multi_h.C_current_down = down_currents_voltages.first;
//     // multi_h.C_voltage_down = down_currents_voltages.second;
//     // arma::mat down_Capacitance = get_capacitance(C_list_down);       // Update the capacitance values because the MOSFETs
//     // multi_h.C_charge_down = down_Capacitance % multi_h.C_voltage_down;  // element-wise multiplication of two objects (Schur product)
//     // // multi_h.LHS_down = down_matrixs.first;
//     // // multi_h.RHS_down = down_matrixs.second;
//     // multi_h.Capacitance_down = down_Capacitance;
//     // multi_h.C_list_down = C_list_down;

//     return multi_h;
// }

// void save_threads_time(const std::chrono::time_point<std::chrono::high_resolution_clock> &t1, const std::chrono::time_point<std::chrono::high_resolution_clock> &t2){

//     std::ofstream file("Threads time mid.csv");
//     std::ofstream file2("Threads time up.csv");
//     std::ofstream file3("Threads time down.csv");
//     std::ofstream file4("simulation end time.csv");

//     std::vector<std::chrono::duration<double, std::milli>> bg_mid(begin_mid.size()), bg_up(begin_up.size()), bg_down(begin_down.size()),
//         ed_mid(begin_mid.size()), ed_up(begin_up.size()), ed_down(begin_down.size());

//     for(int i = 0; i < begin_mid.size(); ++i){
//         bg_mid.at(i) = begin_mid.at(i) - t1;
//         ed_mid.at(i) = end_mid.at(i) - t1;
//     }

//     for(int i = 0; i < begin_down.size(); ++i){
//         bg_down.at(i) = begin_down.at(i) - t1;
//         ed_down.at(i) = end_down.at(i) - t1;
//     }

//     for(int i = 0; i < begin_up.size(); ++i){
//         bg_up.at(i) = begin_up.at(i) - t1;
//         ed_up.at(i) = end_up.at(i) - t1;
//     }

//     // file << "Begin Mid, Begin Up, Begin Down, End Mid, End Up, End Down" << std::endl;
//     file << "Begin Mid, End Mid" << std::endl;
//     file2 << "Begin Up, End Up" << std::endl;
//     file3 << "Begin Down, End Down" << std::endl;

//     // for(int i = 0; i < begin_mid.size(); ++i){
//     //     file << bg_mid.at(i).count() << ", " << bg_up.at(i).count() << ", " << bg_down.at(i).count() << ", "
//     //     << ed_mid.at(i).count() << ", " << ed_up.at(i).count() << ", " << ed_down.at(i).count() << std::endl;
//     // }
//     for(int i = 0; i < begin_mid.size(); ++i){
//         file << bg_mid.at(i).count() << ", " << ed_mid.at(i).count() << std::endl;
//     }
//     for(int i = 0; i < begin_up.size(); ++i){
//         file2 << bg_up.at(i).count() << ", " << ed_up.at(i).count() << std::endl;
//     }
//     for(int i = 0; i < begin_down.size(); ++i){
//         file3 << bg_down.at(i).count() << ", " << ed_down.at(i).count() << std::endl;
//     }

//     file4 << "Simulation end time" << std::endl;
//     std::chrono::duration<double, std::milli> total_time = (t2 - t1);
//     file4 << total_time.count() << std::endl;
// }
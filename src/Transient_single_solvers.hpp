#pragma once

#include <iostream>
#include <armadillo>
#include <cmath>
#include <chrono>
#include <string>
#include <variant>
#include <vector>
#include <algorithm>

#include "sim_variables.hpp"
#include "Transient_calcs.hpp"

struct single_timestep
{
    arma::vec C_current;
    arma::vec C_voltage;
    arma::vec C_charge;
    arma::vec Capacitance;

    double h{};             // time step
    double t{};             // time

    arma::vec RHS;
    arma::mat LHS;
    arma::vec solution;

    bool TMAX_reach;
    bool TMIN_reach;

    std::vector<Capacitor> C_list;
};

struct single_Truncation_error
{
    arma::vec LTE_bound{};
    arma::vec LTE_BE_h{};
};

// Calculate the truncation error (LTE)
single_Truncation_error single_LTE_BE_calculation(const single_timestep &single_h, const std::vector<Transient> &vec_trans){
    single_Truncation_error LTE;

    arma::vec pre_current = vec_trans.back().C_current;
    arma::vec pre_charge = vec_trans.back().C_charge;

    // max(|ik+1|,|ik|) and max(|qk+1|,|qk|)
    arma::vec max_current = arma::max(arma::abs(single_h.C_current), arma::abs(pre_current));
    arma::vec max_charge = arma::max(arma::abs(single_h.C_charge), arma::abs(pre_charge));

    // LTE current = ABSTOL + RELTOL * max(|ik+1|,|ik|) in SPICE book
    arma::vec LTE_current = ABSTOL + RELTOL * max_current;

    // LTE charge = RELTOL * max(|qk+1|,|qk|, chgtol) / h(k) in SPICE book
    arma::vec CHGTOL_MAT = arma::vec(max_charge.n_rows, max_charge.n_cols).fill(CHGTOL);
    arma::vec LTE_charge = RELTOL * arma::max(max_charge, CHGTOL_MAT) / single_h.h;

    // LTE bound = max(LTE_current, LTE_charge)
    LTE.LTE_bound = arma::max(LTE_current, LTE_charge);

    // |i(k+1) - i(k)|
    arma::vec BEcur_diff = single_h.C_current - pre_current;
    
    // tol_h = 2C / |i(k+1) - i(k)| * LTE_bound
    LTE.LTE_BE_h = ((2 * single_h.Capacitance) / arma::abs(BEcur_diff)) % LTE.LTE_bound;

    return LTE;
}

single_timestep single_solution_solver(const double &h, const Transient &trans, const CKTcircuit &ckt, const TransientSimulator &trans_sim){
    
    single_timestep single_h;
    single_h.h = h;
    single_h.t = trans_sim.vec_trans.back().time_trans + h;
    std::vector<Capacitor> C_list_copy = trans_sim.trans_config.C_list;

    if(trans_sim.trans_config.non_linear){
        single_h.solution = NewtonRaphson_system(ckt, h, trans.mode, single_h.t, C_list_copy, trans_sim.vec_trans.back().solution);
        std::pair<arma::vec, arma::vec> currents_voltages = get_currents_voltages(C_list_copy, h, single_h.solution, trans_sim.vec_trans.back().solution);
        single_h.C_current = currents_voltages.first;
        single_h.C_voltage = currents_voltages.second;
        single_h.Capacitance = get_capacitance(C_list_copy);
        single_h.C_charge = single_h.Capacitance % single_h.C_voltage; // element-wise multiplication of two objects (Schur product)
        single_h.C_list = C_list_copy;
    }
    else if(trans_sim.trans_config.non_linear == false){
        std::pair<arma::mat, arma::vec> matrices;
        matrices = Dynamic(ckt, h, trans_sim.vec_trans.back().solution, trans.mode, single_h.t);
        single_h.solution = arma::solve(matrices.first, matrices.second, arma::solve_opts::fast);
        std::pair<arma::vec, arma::vec> currents_voltages = get_currents_voltages(C_list_copy, h, single_h.solution, trans_sim.vec_trans.back().solution);
        single_h.C_current = currents_voltages.first;
        single_h.C_voltage = currents_voltages.second;
        single_h.Capacitance = get_capacitance(C_list_copy);
        single_h.C_charge = single_h.Capacitance % single_h.C_voltage; // element-wise multiplication of two objects (Schur product)
        single_h.C_list = C_list_copy;
    }
    else{
        std::cerr << "Error in single_solution_solver function: trans_sim.trans_config.non_linear" << std::endl;
        exit(1);
    }

    return single_h;
}

bool single_LTE_check(single_Truncation_error &LTE, const single_timestep &single_h,
                        double &temp_h, const TransientSimulator &trans_sim){
    if(NR_ITE < ITL4){

       LTE = single_LTE_BE_calculation(single_h, trans_sim.vec_trans);

       if(LTE.LTE_BE_h.min() * TRTOL < 0.9 * single_h.h){
           // Reject the solution, hn = hn+1 and recompute the new hn
           temp_h = LTE.LTE_BE_h.min() * TRTOL;
           return false;
       }
       else{
           // Accept the solution
           temp_h = std::min({(LTE.LTE_BE_h.min() * TRTOL), 2.0 * single_h.h, trans_sim.trans_config.h_MAX});

           return true;
       }
    }

    else{
        // Reject the solution, hn = hn/8 and recompute the new hn
        temp_h = single_h.h / 8.0;

        if(temp_h > trans_sim.trans_config.h_MIN){
            return false;
        }
        else{
            std::cerr << "The time step is too small (Reach minimum!)" << std::endl;
            exit(1);
        }
    }
}

arma::vec single_next_h(Transient &trans, const CKTcircuit &ckt, const TransientSimulator &trans_sim){

    arma::vec solution;
    double temp_h = trans_sim.vec_trans.back().next_h;                     // The temporary time step for this function

    bool LTE_check= false;

    single_timestep single_h;
    single_Truncation_error LTE;

    do{
        single_h = single_solution_solver(temp_h, trans, ckt, trans_sim);

        total_timepoint += 1;

        LTE_check = single_LTE_check(LTE, single_h, temp_h, trans_sim);

    }while(LTE_check == false);

    trans.h = single_h.h;
    trans.next_h = temp_h;
    trans.solution = single_h.solution;
    trans.C_list = single_h.C_list;
    trans.Capacitance = single_h.Capacitance;
    trans.C_current = single_h.C_current;
    trans.C_voltage = single_h.C_voltage;
    trans.C_charge = single_h.C_charge;

    return solution;
}
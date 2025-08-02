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
#include "circuit_parser.hpp"
#include "XB_timer.hpp"
#include "device.hpp"
#include "matrix.hpp"
#include "CKT.hpp"
#include "CKTsetup.hpp"
#include "global.hpp"
#include "Newton.hpp"

// Macro to enable debug code
#define ARMA_PRINT(var, label) \
    if (isDebugMode())         \
    {                          \
        (var).print(label);    \
    }
struct TransientConfig;
struct Transient;
struct Truncation_error;
struct TransientSimulator;

/*
    Clearly distinguishes between data that never changes and data that evolves in Transient simulation.
    TransientConfig is a constant data structure that never changes.
*/
struct TransientConfig {
    double t_start{};
    double t_end{};
    double init_h{};
    double h_MAX{};
    double h_MIN{};
    bool timestep_control;          // true for variable time step, false for fixed time step
    bool non_linear;                // true for non-linear solver, false for linear solver
};

struct CapacitanceState
{
    std::vector<double> CapCharge;
    std::vector<double> CapCurrent;
};

struct Transient
{
    double h{};
    double next_h{};
    double time_trans{};
    int mode{};
    int trans_count{};
    arma::vec solution;
    arma::mat LHS;
    arma::vec RHS;
    CapacitanceState CapState; // Capacitance state (including cap, bsim4, etc.)

};

struct Truncation_error
{
    arma::vec LTE_bound_mid{};
    arma::vec LTE_bound_up{};
    arma::vec LTE_bound_down{};

    arma::vec LTE_BE_mid_h{};
    arma::vec LTE_BE_up_h{};
    arma::vec LTE_BE_down_h{};
};

struct TransientSimulator
{   
    const TransientConfig trans_config;      // Transient configuration(never changes)
    std::vector<Transient> vec_trans;        // Transient history
    std::deque<double> breakpoints;          // Breakpoints for the transient simulation
    bool trans_end = false;                  // End of transient simulation

    TransientSimulator(const TransientConfig& config) : trans_config(config) {}
};

arma::vec get_capacitance(const std::vector<Capacitor> &C_list)
{
    arma::vec Capacitance = arma::zeros(C_list.size(), 1);
    for (size_t i = 0; i < C_list.size(); i++)
    {
        Capacitance(i, 0) = C_list.at(i).value;
    }
    return Capacitance;
}

TransientSimulator Transsetup(const CircuitParser &parser, const CKTcircuit &ckt)
{
    TransientConfig config;
    config.t_end = parser.double_t_end;

    // Only closed when the user turns off step control or there are no dynamic elements
    if(parser.timestep_control == false || ckt.num_of_states == 0){
        config.h_MAX = 1e-9;
        config.init_h = config.h_MAX;
        config.timestep_control = false; // Turn off the time step control
    }
    // Turn on the time step control
    else if (ckt.pulse_num > 0)
    {
        for (const auto &pulse : ckt.CKTelements.pulseVoltages)
        {
            double step = std::min(pulse.tr, pulse.tf);
            config.h_MAX = parser.double_init_h; // 5 is rmax in spice opus
            // trans.h_MIN = std::max(1e-9 * parser.double_init_h, 1e-14);    // 1e-9 is rmin in spice opus
            config.h_MIN = 1e-14;

            if (pulse.td == 0)
            {
                // config.init_h = parser.double_init_h / 100; // initial time step, fs = 0.25 from spice opus
                config.init_h = std::min(parser.double_init_h, parser.double_t_end / 100) / 100; // delta=MIN(ckt->CKTfinalTime/100,ckt->CKTstep)/10; from ngspice
            }
            else
            {
                config.init_h = std::min(pulse.td * 0.25, parser.double_init_h * 0.25);
            }
        }

        config.timestep_control = true;  // Turn on the time step control
    }
    else
    {
        config.init_h = parser.double_init_h / 100; // initial time step
        config.h_MAX = parser.double_init_h;        // maximum time step
        // config.h_MIN = std::max(1e-9 * parser.double_init_h, 1e-14);     // minimum time step
        config.h_MIN = 1e-14;
        config.timestep_control = true;  // Turn on the time step control
    }


    // Check if the transient simulation is non-linear
    config.non_linear = CKTisNonLinear(ckt.CKTelements);

    // Setup the initial capacitance list
    TransientSimulator trans_sim(config);

    return trans_sim;
}

/*
    1. Start of rise time: t=td. This is when the voltage begins to transition from V1 to V2​.
    2. End of rise time: t=td+tr. This is when the voltage finishes transitioning to V2​.
    3. Start of fall time: t=td+tr+tpw. This is when the voltage begins to transition back to V1
    4​. End of fall time: t=td+tr+tpw+tf. This is when the voltage finishes transitioning back to V1​.
*/
std::deque<double> get_breakpoints(const CKTcircuit &ckt, const TransientSimulator &trans_sim)
{

    std::deque<double> breakpoints;

    if (ckt.CKTelements.pulseVoltages.empty())
    {
        return breakpoints;
    }

    for (const auto &pulse : ckt.CKTelements.pulseVoltages){
        for(double cycle_start = 0; cycle_start < trans_sim.trans_config.t_end; cycle_start += pulse.per){

            double cycle_times[] = {
                pulse.td + cycle_start, 
                pulse.td + pulse.tr + cycle_start, 
                pulse.td + pulse.tr + pulse.pw + cycle_start, 
                pulse.td + pulse.tr + pulse.pw + pulse.tf + cycle_start
            };

            for(double t : cycle_times){
                if(t <= trans_sim.trans_config.t_end){
                    breakpoints.push_back(t);
                }
            }
            
        }
    }

    std::sort(breakpoints.begin(), breakpoints.end());

    if (breakpoints.at(0) == 0)
    {
        breakpoints.pop_front();
    }

    if (breakpoints.back() != trans_sim.trans_config.t_end)
    {
        breakpoints.push_back(trans_sim.trans_config.t_end);
    }

    // Ensure the difference between consecutive breakpoints is at least 1e-13
    std::deque<double> filtered_breakpoints;

    if (!breakpoints.empty())
    {
        filtered_breakpoints.push_back(breakpoints.front());
        for (size_t i = 1; i < breakpoints.size(); ++i)
        {
            if (breakpoints.at(i) - filtered_breakpoints.back() >= 1e-13)
            {
                filtered_breakpoints.push_back(breakpoints[i]);
            }
        }
    }

    return filtered_breakpoints;
}

// Assigning the stamp matrices for dynamic and non-linear components and update the LHS and RHS matrices

// std::pair<arma::vec, arma::vec> get_currents_voltages(const std::vector<Capacitor> &pre_C_list, const double h, const arma::vec &solution, const arma::vec &pre_solution)
// {

//     // h/C * i = u(k+1) - u(k)

//     int G_rows = pre_C_list.size();
//     int G_cols = G_rows;
//     double vol{}, pre_vol{}, delta_vol{}; // Delta voltage across the capacitor u(k+1) - u(k)

//     arma::vec current_matrix = arma::zeros(G_rows, 1); // Currents matrix
//     arma::vec delta_v = arma::zeros(G_rows, 1);        // Delta voltage matrix
//     arma::vec G_vec(pre_C_list.size());                // Initialize arma::vec with the size of C_list
//     arma::vec volt = arma::zeros(G_rows, 1);

//     for (size_t i = 0; i < pre_C_list.size(); ++i)
//     {

//         G_vec(i) = h / pre_C_list.at(i).value;

//         if (pre_C_list.at(i).nodePos == 0)
//         {

//             vol = solution(pre_C_list.at(i).nodeNeg - 1, 0);
//             pre_vol = pre_solution(pre_C_list.at(i).nodeNeg - 1, 0);
//             delta_vol = vol - pre_vol; // u(k+1) - u(k)
//         }
//         else if (pre_C_list.at(i).nodeNeg == 0)
//         {

//             vol = solution(pre_C_list.at(i).nodePos - 1, 0);
//             pre_vol = pre_solution(pre_C_list.at(i).nodePos - 1, 0);
//             delta_vol = vol - pre_vol; // u(k+1) - u(k)
//         }
//         else
//         {

//             vol = solution(pre_C_list.at(i).nodePos - 1, 0) - solution(pre_C_list.at(i).nodeNeg - 1, 0);
//             pre_vol = pre_solution(pre_C_list.at(i).nodePos - 1, 0) - pre_solution(pre_C_list.at(i).nodeNeg - 1, 0);
//             delta_vol = vol - pre_vol; // u(k+1) - u(k)
//         }

//         delta_v(i, 0) = delta_vol; // It may be negative!
//         volt(i, 0) = vol;
//         // std::cout << "c: " << pre_C_list.at(i).value << std::endl;
//         // std::cout << "C/h: " << pre_C_list.at(i).value/h << std::endl;
//         // std::cout << "delta_vol: " << delta_vol << std::endl;
//         // std::cout << "C/h * delta_vol" << pre_C_list.at(i).value/h * delta_vol << std::endl;
//     }

//     arma::mat G_matrix = arma::diagmat(G_vec); // Diagonal matrix

//     // G_matrix.print("G_matrix =");
//     // delta_v.print("delta_v =");

//     current_matrix = arma::solve(G_matrix, delta_v, arma::solve_opts::fast);

//     return {current_matrix, volt};
// }

// Function to control debug mode
double cond(double R)
{
    return 1 / R;
}

void history_trans_update(Transient &trans, TransientSimulator &trans_sim)
{
    // Update transient history
    trans_sim.vec_trans.push_back(std::move(trans));
}


double average_vec(std::vector<double> &v)
{
    double sum = std::accumulate(v.begin(), v.end(), 0.0);
    return sum / v.size();
}
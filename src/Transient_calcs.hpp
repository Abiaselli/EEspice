#pragma once
// #define ARMA_DONT_USE_WRAPPER
// #define ARMA_USE_MKL_ALLOC

// #define ARMA_USE_SUPERLU
// #define COOT_DONT_USE_OPENCL
// #define COOT_USE_CUDA

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

#include "device.hpp"
#include "matrix.hpp"
#include "CKT.hpp"
#include "global.hpp"

// Macro to enable debug code

#define ARMA_PRINT(var, label) \
    if (isDebugMode())         \
    {                          \
        (var).print(label);    \
    }

/*void setDebugMode(bool mode);
bool isDebugMode();
void R_assigner(int node_x, int node_y, double G, arma::mat &LHS);
void Vs_assigner(int node_x, int node_y, double V_value, arma::mat &LHS, arma::vec &RHS);
double convertToValue(const std::string &valueStr);
arma::mat branch_ext(arma::mat M, int node_x, int node_y);
void Is_assigner(double node_x, double node_y, double I, arma::vec &RHS);
void C_assigner_BE(int node_x, int node_y, double C, double h, arma::mat &LHS, arma::vec &RHS, const arma::vec &pre_solution, int mode);
int V_pulse_assigner(int node_x, int node_y, double V_value, arma::mat &LHS, arma::vec &RHS);
double V_pulse_value(double V1, double V2, double t1, double td, double tr, double tf, double tpw, double tper);
void VCCS_assigner(int node_x, int node_y, int node_cx, int node_cy, double R, arma::mat &LHS);
double Diode_assigner(int node_x, int node_y, double Is, double VT, arma::mat &LHS, arma::vec &RHS, const arma::vec &solution, int mode);std::pair<arma::mat, arma::vec> DynamicNonLinear(const CKTcircuit &ckt, double h, const arma::vec &pre_solution, int mode, const double time_trans, std::vector<Capacitor> &C_list, int NR_iteration_counter);
arma::vec NewtonRaphson_system(const CKTcircuit &ckt, const arma::vec &pre_solution, const double &h, const int &mode, const double time_trans, std::vector<Capacitor> &C_list);
void dummy_task();*/
/*  Global variables:

*/
// std::vector<Transient> vec_trans;
BS::thread_pool pool(3);
BS::synced_stream sync_out;
XB_Timer timer;
// std::mutex mtx;

/*--------------------------Timeing----------------------------------------------------------------*/
// std::vector<std::chrono::time_point<std::chrono::high_resolution_clock>> begin_mid;
// std::vector<std::chrono::time_point<std::chrono::high_resolution_clock>> begin_up;
// std::vector<std::chrono::time_point<std::chrono::high_resolution_clock>> begin_down;
// std::vector<std::chrono::time_point<std::chrono::high_resolution_clock>> end_mid;
// std::vector<std::chrono::time_point<std::chrono::high_resolution_clock>> end_up;
// std::vector<std::chrono::time_point<std::chrono::high_resolution_clock>> end_down;

/*-----------------------------------------------------------------------------------------------------------------------------------------------------*/

// Function to control debug mode

double cond(double R)
{
    return 1 / R;
}

void history_trans_update(Transient &trans, std::vector<Transient> &vec_trans)
{
    // Update transient history
    vec_trans.push_back(trans);
}

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

struct Transient
{
    int code = 0;
    double t_start = 0;
    double t_end{};
    double h{};
    double init_h{};
    double h_MAX{};
    double h_MIN{};
    double time_trans{};
    int mode{};
    int trans_count{};
    arma::vec solution;
    arma::mat LHS;
    arma::vec RHS;
    arma::vec C_current;
    arma::vec C_voltage;
    arma::vec C_charge;
    arma::vec Capacitance;

    std::vector<Capacitor> C_list;

    void C_list_update()
    {
        for (size_t i = 0; i < C_list.size(); i++)
        {
            C_list.at(i).current = C_current(i, 0);
            C_list.at(i).voltage = C_voltage(i, 0);
            C_list.at(i).charge = C_list.at(i).value * C_list.at(i).voltage;
            C_list.at(i).value = Capacitance(i, 0);
        }
    }

    void get_capacitance()
    {
        Capacitance = arma::zeros(C_list.size(), 1);
        for (size_t i = 0; i < C_list.size(); i++)
        {
            Capacitance(i, 0) = C_list.at(i).value;
        }
    }
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

arma::vec get_capacitance(const std::vector<Capacitor> &C_list)
{
    arma::vec Capacitance = arma::zeros(C_list.size(), 1);
    for (size_t i = 0; i < C_list.size(); i++)
    {
        Capacitance(i, 0) = C_list.at(i).value;
    }
    return Capacitance;
}

void Transsetup(Transient &trans, const CircuitParser &parser, CKTcircuit &ckt)
{
    trans.t_end = parser.double_t_end;

    if (ckt.pulse_num > 0 && (ckt.C_list.size() > 0))
    {
        for (const auto &element : ckt.CKTelements)
        {
            std::visit([&](auto &&arg)
                       {
                           if constexpr (std::is_same_v<std::decay_t<decltype(arg)>, Pulsevoltage>)
                           {

                               double step = std::min(arg.tr, arg.tf);
                               trans.h_MAX = parser.double_init_h; // 5 is rmax in spice opus
                               // trans.h_MIN = std::max(1e-9 * parser.double_init_h, 1e-14);    // 1e-9 is rmin in spice opus
                               trans.h_MIN = 1e-14;

                               if (arg.td == 0)
                               {
                                   trans.init_h = parser.double_init_h / 100; // initial time step, fs = 0.25 from spice opus
                               }
                               else
                               {
                                   trans.init_h = std::min(arg.td * 0.25, parser.double_init_h * 0.25);
                               }
                           } },
                       element.element);
        }
    }
    else if (ckt.no_of_mosfets > 0)
    {
        trans.init_h = parser.double_init_h / 100; // initial time step
        trans.h_MAX = parser.double_init_h;        // maximum time step
        // trans.h_MIN = std::max(1e-9 * parser.double_init_h, 1e-14);     // minimum time step
        trans.h_MIN = 1e-14;
    }
    else
    {

        // trans.h_MAX = trans.t_end / 5000;
        trans.h_MAX = 1e-9;
        trans.init_h = trans.h_MAX;
    }
}

/*
    1. Start of rise time: t=td. This is when the voltage begins to transition from V1 to V2​.
    2. End of rise time: t=td+tr. This is when the voltage finishes transitioning to V2​.
    3. Start of fall time: t=td+tr+tpw. This is when the voltage begins to transition back to V1
    4​. End of fall time: t=td+tr+tpw+tf. This is when the voltage finishes transitioning back to V1​.
*/
std::deque<double> get_breakpoints(const CKTcircuit &ckt, const Transient &trans)
{

    std::deque<double> breakpoints;

    if (ckt.pulse_num == 0)
    {
        return breakpoints;
    }

    for (const auto &element : ckt.CKTelements)
    {
        std::visit([&](auto &&arg)
                   {
            if constexpr(std::is_same_v<std::decay_t<decltype(arg)>, Pulsevoltage>){
                
                for(double cycle_start = 0; cycle_start < trans.t_end; cycle_start += arg.per){

                    double cycle_times[] = {
                        arg.td+cycle_start, 
                        arg.td+arg.tr+cycle_start, 
                        arg.td+arg.tr+arg.pw+cycle_start, 
                        arg.td+arg.tr+arg.pw+arg.tf+cycle_start
                    };

                    for(double t : cycle_times){
                        if(t <= trans.t_end){
                            breakpoints.push_back(t);
                        }
                    }
                    
                }

            } }, element.element);
    }

    std::sort(breakpoints.begin(), breakpoints.end());

    if (breakpoints.at(0) == 0)
    {
        breakpoints.pop_front();
    }

    if (breakpoints.back() != trans.t_end)
    {
        breakpoints.push_back(trans.t_end);
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
std::pair<arma::mat, arma::vec> DynamicNonLinear(const CKTcircuit &ckt, double h, const arma::vec &pre_solution, int mode, const double time_trans, std::vector<Capacitor> &C_list, int NR_iteration_counter, std::vector<Transient> &vec_trans)
{

    arma::mat LHS = ckt.cktdematrix->get_init_LHS();
    arma::vec RHS = ckt.cktdematrix->get_init_RHS();
    // LHS.print("in_LHS matrix =");
    // RHS.print("RHS matrix =");

    for (const auto &element : ckt.CKTelements)
    {
        std::visit([&](auto &&arg)
                   {
                       if constexpr (std::is_same_v<std::decay_t<decltype(arg)>, Capacitor>)
                       {
                           // Linear Capacitor
                           C_assigner_BE(arg.nodePos, arg.nodeNeg, arg.value, h, LHS, RHS, vec_trans.back().solution, mode);
                       }
                       else if constexpr (std::is_same_v<std::decay_t<decltype(arg)>, Pulsevoltage>)
                       {

                           double val_pulse = V_pulse_value(arg.V1, arg.V2, time_trans, arg.td, arg.tr, arg.tf, arg.pw, arg.per);

                           RHS.row(arg.RHS_locate - 1).col(0) += val_pulse;
                       }
                       else if constexpr (std::is_same_v<std::decay_t<decltype(arg)>, NMOS>)
                       {
                           NMOS_assigner(arg.id, arg.node_vd, arg.node_vg, arg.node_vs, arg.node_vb, arg.W, arg.L, h, pre_solution, ckt.T_nodes, LHS, RHS, mode, C_list, NR_iteration_counter);
                       }
                       else if constexpr (std::is_same_v<std::decay_t<decltype(arg)>, Diode>)
                       {
                           Diode_assigner(arg.nodePos, arg.nodeNeg, arg.Is, arg.VT, LHS, RHS, pre_solution, mode);
                       }
                       else if constexpr (std::is_same_v<std::decay_t<decltype(arg)>, PMOS>)
                       {
                           PMOS_assigner(arg.id, arg.node_vd, arg.node_vg, arg.node_vs, arg.node_vb, arg.W, arg.L, h, pre_solution, ckt.T_nodes, LHS, RHS, mode, C_list, NR_iteration_counter);
                       } },
                   element.element);
    }

    return {LHS, RHS};
}

bool isConverge(const std::vector<arma::vec> &NR_solutions, const CKTcircuit &ckt, const int &NR_iteration_counter)
{

    const arma::vec &pre_solution = NR_solutions.at(NR_iteration_counter - 2);
    const arma::vec &current_solution = NR_solutions.at(NR_iteration_counter - 1);
    const arma::vec &next_solution = NR_solutions.at(NR_iteration_counter);

    if (pre_solution.n_rows != current_solution.n_rows || pre_solution.n_rows != next_solution.n_rows || current_solution.n_rows != next_solution.n_rows)
    {
        std::cerr << "The size of pre_solution, current_solution, next_solution is not the same in isConverge function." << std::endl;
        exit(1);
    }

    arma::vec pre_voltages = pre_solution.submat(0, 0, ckt.T_nodes - 1, 0);
    arma::vec current_voltages = current_solution.submat(0, 0, ckt.T_nodes - 1, 0);
    arma::vec next_voltages = next_solution.submat(0, 0, ckt.T_nodes - 1, 0);

    arma::vec pre_current = pre_solution.submat(ckt.T_nodes, 0, pre_solution.n_rows - 1, 0);
    arma::vec current_current = current_solution.submat(ckt.T_nodes, 0, current_solution.n_rows - 1, 0);
    arma::vec next_current = next_solution.submat(ckt.T_nodes, 0, next_solution.n_rows - 1, 0);

    // |v(k+1)-v(k)| <= RELTOL*max(|v(k+1)|,|v(k)|) + VNTOL

    arma::vec Vmax_next_current = arma::max(arma::abs(next_voltages), arma::abs(current_voltages));
    arma::vec Vdelta_next_current = arma::abs(next_voltages - current_voltages);
    arma::vec tolerance_1 = RELTOL * Vmax_next_current + VNTOL;
    // tolerance_1.print("tolerance_1 =");

    // If |v(k+1)-v(k)| is bigger arma::any will return true and the function will return false.
    if (arma::any(Vdelta_next_current > tolerance_1))
    {
        return false;
    }

    // |i(k+1)-i(k)| <= RELTOL*max(|i(k+1)|,|i(k)|) + ABSTOL

    arma::vec Imax_next_current = arma::max(arma::abs(next_current), arma::abs(current_current));
    arma::vec Idelta_next_current = arma::abs(next_current - current_current);
    arma::vec tolerance_2 = RELTOL * Imax_next_current + ABSTOL;
    // tolerance_2.print("tolerance_2 =");

    if (arma::any(Idelta_next_current > tolerance_2))
    {
        return false;
    }

    // |v(k) - v(k-1)| <= RELTOL*max(|v(k)|,|v(k-1)|) + VNTOL

    arma::vec Vmax_current_pre = arma::max(arma::abs(current_voltages), arma::abs(pre_voltages));
    arma::vec Vdelta_current_pre = arma::abs(current_voltages - pre_voltages);
    arma::vec tolerance_3 = RELTOL * Vmax_current_pre + VNTOL;

    if (arma::any(Vdelta_current_pre > tolerance_3))
    {
        return false;
    }

    // |i(k) - i(k-1)| <= RELTOL*max(|i(k)|,|i(k-1)|) + ABSTOL

    arma::vec Imax_current_pre = arma::max(arma::abs(current_current), arma::abs(pre_current));
    arma::vec Idelta_current_pre = arma::abs(current_current - pre_current);
    arma::vec tolerance_4 = RELTOL * Imax_current_pre + ABSTOL;

    if (arma::any(Idelta_current_pre > tolerance_4))
    {
        return false;
    }

    // // |v(k+1) - v(k-1)| <= √|v(k) - v(k-1)|^2 + |v(k+1) - v(k)|^2

    //     arma::mat Vdelta_next_pre = arma::abs(next_voltages - pre_voltages);
    //     arma::mat V_Perpendicular = arma::sqrt(arma::pow(arma::abs(Vdelta_current_pre),2) + arma::pow(arma::abs(Vdelta_next_current),2));

    //     if (arma::any(arma::vectorise(Vdelta_next_pre > V_Perpendicular)) )
    //     {
    //         return false;
    //     }

    // // |i(k+1) - i(k-1)| <= √|i(k) - i(k-1)|^2 + |i(k+1) - i(k)|^2

    //     arma::mat Idelta_next_pre = arma::abs(next_current - pre_current);
    //     arma::mat I_Perpendicular = arma::sqrt(arma::pow(arma::abs(Idelta_current_pre),2) + arma::pow(arma::abs(Idelta_next_current),2));

    //     if (arma::any(arma::vectorise(Idelta_next_pre > I_Perpendicular)))
    //     {
    //         return false;
    //     }

    return true;
}

// Newton Raphson system solver for non-linear and dynamic elements
arma::vec NewtonRaphson_system(const CKTcircuit &ckt, const arma::vec &pre_solution, const double &h, const int &mode, const double time_trans, std::vector<Capacitor> &C_list, std::vector<Transient> &vec_trans)
{

    // std::cout << "Enter NewtonRaphson system" << std::endl;
    DEBUG_PRINT("Enter NewtonRaphson system");

    int NR_iteration_counter = 0;
    bool isconverge = false;
    arma::vec solution = pre_solution;
    std::pair<arma::mat, arma::vec> matrices;

    std::vector<arma::vec> NR_solutions(100);
    NR_solutions[0] = pre_solution;

    for (int i = 1; i < 3; i++)
    {
        matrices = DynamicNonLinear(ckt, h, solution, mode, time_trans, C_list, NR_iteration_counter, vec_trans);
        // const arma::mat &LHS = matrices.first;
        // const arma::mat &RHS = matrices.second;

        // Solve Ax = b
        // J(v) * x(k+1) = [J(v)]x(k) - f(x(k))
        // solution = arma::solve(matrices.first, matrices.second, arma::solve_opts::fast);
        solution = arma::solve(matrices.first, matrices.second);
        NR_iteration_counter += 1;
        NR_solutions.at(NR_iteration_counter) = solution;
    }

    isconverge = isConverge(NR_solutions, ckt, NR_iteration_counter);

    while (!isconverge)
    {

        if (NR_iteration_counter >= 99)
        {
            std::cerr << "The Newton Raphson method did not converge after 100 iterations." << std::endl;
            // exit(1);
            return solution;
        }

        matrices = DynamicNonLinear(ckt, h, solution, mode, time_trans, C_list, NR_iteration_counter, vec_trans);

        // Solve Ax = b
        // J(v) * x(k+1) = [J(v)]x(k) - f(x(k))
        // solution = arma::solve(matrices.first, matrices.second, arma::solve_opts::fast);
        solution = arma::solve(matrices.first, matrices.second);
        NR_iteration_counter += 1;
        NR_solutions.at(NR_iteration_counter) = solution;

        isconverge = isConverge(NR_solutions, ckt, NR_iteration_counter);
    }

    return solution;
}

std::pair<arma::vec, arma::vec> get_currents_voltages(const std::vector<Capacitor> &pre_C_list, const double h, const arma::vec &solution, const arma::vec &pre_solution)
{

    // h/C * i = u(k+1) - u(k)

    int G_rows = pre_C_list.size();
    int G_cols = G_rows;
    double vol{}, pre_vol{}, delta_vol{}; // Delta voltage across the capacitor u(k+1) - u(k)

    arma::vec current_matrix = arma::zeros(G_rows, 1); // Currents matrix
    arma::vec delta_v = arma::zeros(G_rows, 1);        // Delta voltage matrix
    arma::vec G_vec(pre_C_list.size());                // Initialize arma::vec with the size of C_list
    arma::vec volt = arma::zeros(G_rows, 1);

    for (size_t i = 0; i < pre_C_list.size(); ++i)
    {

        G_vec(i) = h / pre_C_list.at(i).value;

        if (pre_C_list.at(i).nodePos == 0)
        {

            vol = solution(pre_C_list.at(i).nodeNeg - 1, 0);
            pre_vol = pre_solution(pre_C_list.at(i).nodeNeg - 1, 0);
            delta_vol = vol - pre_vol; // u(k+1) - u(k)
        }
        else if (pre_C_list.at(i).nodeNeg == 0)
        {

            vol = solution(pre_C_list.at(i).nodePos - 1, 0);
            pre_vol = pre_solution(pre_C_list.at(i).nodePos - 1, 0);
            delta_vol = vol - pre_vol; // u(k+1) - u(k)
        }
        else
        {

            vol = solution(pre_C_list.at(i).nodePos - 1, 0) - solution(pre_C_list.at(i).nodeNeg - 1, 0);
            pre_vol = pre_solution(pre_C_list.at(i).nodePos - 1, 0) - pre_solution(pre_C_list.at(i).nodeNeg - 1, 0);
            delta_vol = vol - pre_vol; // u(k+1) - u(k)
        }

        delta_v(i, 0) = delta_vol; // It may be negative!
        volt(i, 0) = vol;
        // std::cout << "c: " << pre_C_list.at(i).value << std::endl;
        // std::cout << "C/h: " << pre_C_list.at(i).value/h << std::endl;
        // std::cout << "delta_vol: " << delta_vol << std::endl;
        // std::cout << "C/h * delta_vol" << pre_C_list.at(i).value/h * delta_vol << std::endl;
    }

    arma::mat G_matrix = arma::diagmat(G_vec); // Diagonal matrix

    // G_matrix.print("G_matrix =");
    // delta_v.print("delta_v =");

    current_matrix = arma::solve(G_matrix, delta_v, arma::solve_opts::fast);

    return {current_matrix, volt};
}

// New timestep options
void timestep_options(double &temp_h, double &next_h_up, double &next_h_down, const Transient &trans, bool &TMAX_reach, bool &TMIN_reach)
{
    if (temp_h >= trans.h_MAX)
    {
        temp_h = trans.h_MAX;
        TMAX_reach = true;
    }
    else if (temp_h <= trans.h_MIN)
    {
        temp_h = trans.h_MIN;
        TMIN_reach = true;
    }

    next_h_up = temp_h * 2;
    next_h_down = temp_h / 2;

    if (next_h_up > trans.h_MAX)
    {
        next_h_up = trans.h_MAX;
        // TMAX_reach = true;
    }
    else if (next_h_up < trans.h_MIN)
    {
        next_h_down = trans.h_MIN;
        // TMAX_reach = true;
    }

    if (next_h_down > trans.h_MAX)
    {
        next_h_up = trans.h_MAX;
        // TMAX_reach = true;
    }
    else if (next_h_down < trans.h_MIN)
    {
        next_h_down = trans.h_MIN;
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
        arma::vec local_solution = NewtonRaphson_system(ckt, vec_trans.back().solution, h, 1, t, C_list_copy, vec_trans);

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

void dummy_task()
{
    // std::cout << "dummy task" << std::endl;
    std::this_thread::sleep_for(std::chrono::microseconds(100));
}

double average_vec(std::vector<double> &v)
{
    double sum = std::accumulate(v.begin(), v.end(), 0.0);
    return sum / v.size();
}
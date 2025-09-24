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
#include "simulation_exceptions.hpp"

#include "CKT.hpp"
#include "device.hpp"
#include "bsim4v82/bsim4v82.hpp"
#include "bsim4v82/bsim4v82stamp.hpp"
#include "NonLinearBatch.hpp"

std::pair<arma::mat, arma::vec> Dynamic(const CKTcircuit &ckt, const double h, const arma::vec &pre_global_solution, const int mode, const double time_trans, SimulationTime &simTime)
{
    ScopedTimer loadTimer(simTime.matrix_load_time);
    arma::mat LHS = ckt.cktdematrix->get_init_LHS();
    arma::vec RHS = ckt.cktdematrix->get_init_RHS();

    for (const auto &cap : ckt.CKTelements.capacitors)
    {   
        // Linear Capacitor
        C_assigner_BE(cap.nodePos, cap.nodeNeg, cap.value, h, LHS, RHS, pre_global_solution, mode);
    }
    for (const auto &pulse : ckt.CKTelements.pulseVoltages)
    {
        double val_pulse = V_pulse_value(pulse.V1, pulse.V2, time_trans, pulse.td, pulse.tr, pulse.tf, pulse.pw, pulse.per);
        RHS(pulse.RHS_locate) += (val_pulse - pulse.V1);
    }
    for (const auto &sin : ckt.CKTelements.sinVoltages){
        double val_sin = V_sin_value(sin.vo, sin.va, sin.freq, sin.td, sin.theta, sin.phase_rad, time_trans);
        RHS(sin.RHS_locate) += (val_sin - sin.vo);
    }
    return {LHS, RHS};
}

std::pair<arma::mat, arma::vec> NonLinear(CKTcircuit &ckt, const arma::vec &pre_NR_solution, 
    const std::pair<arma::mat, arma::vec> &matrixes, const Modelmap &modmap, double h)
{
    ScopedTimer loadTimer(ckt.sim_stats.simTime.matrix_load_time);
    arma::mat LHS = matrixes.first;
    arma::vec RHS = matrixes.second;
    // LHS.print("in_LHS matrix =");
    // RHS.print("RHS matrix =");

    for (auto &nmos : ckt.CKTelements.nmos){
        if (nmos.modelType == MosfetModelType::LEVEL1){
            const NMOSModel nmosModel = modmap.nmosModels.at(nmos.modelName);
            NMOS_assigner(nmos.id, nmos.node_vd, nmos.node_vg, nmos.node_vs, nmos.node_vb, nmos.W, nmos.L, pre_NR_solution, ckt.T_nodes, LHS, RHS, nmosModel);
        }
        else if (nmos.modelType == MosfetModelType::BSIM4V82){
            const bsim4::BSIM4model &b4model = *nmos.bsim4v82Instance.BSIM4modPtr;
            bsim4::BSIM4V82 &b4instance = nmos.bsim4v82Instance;
            bsim4::BSIM4load(ckt, b4model, b4instance, ckt.spiceCompatible, pre_NR_solution, ckt.CKTtemp, ckt.CKTgmin, LHS, RHS);
        }
        else{
            throw DeviceException("Error: NMOS model type is not supported!", "UNSUPPORTED_NMOS_MODEL");
        }
    }
    for (auto &pmos : ckt.CKTelements.pmos){
        if (pmos.modelType == MosfetModelType::LEVEL1){
            const PMOSModel pmosModel = modmap.pmosModels.at(pmos.modelName);
            PMOS_assigner(pmos.id, pmos.node_vd, pmos.node_vg, pmos.node_vs, pmos.node_vb, pmos.W, pmos.L, pre_NR_solution, ckt.T_nodes, LHS, RHS, pmosModel);
        }
        else if (pmos.modelType == MosfetModelType::BSIM4V82){
            const bsim4::BSIM4model &b4model = *pmos.bsim4v82Instance.BSIM4modPtr;
            bsim4::BSIM4V82 &b4instance = pmos.bsim4v82Instance;
            bsim4::BSIM4load(ckt, b4model, b4instance, ckt.spiceCompatible, pre_NR_solution, ckt.CKTtemp, ckt.CKTgmin, LHS, RHS);
        }
        else{
            throw DeviceException("Error: PMOS model type is not supported!", "UNSUPPORTED_PMOS_MODEL");
        }
    }
    for (const auto &diode : ckt.CKTelements.diodes){
        Diode_assigner(diode.nodePos, diode.nodeNeg, diode.Is, diode.VT, LHS, RHS, pre_NR_solution);
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
        throw MatrixException("The size of pre_solution, current_solution, next_solution is not the same in isConverge function.", "SOLUTION_SIZE_MISMATCH");
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

// Solver for the Newton Raphson method
arma::vec solveDense(const arma::mat &LHS, const arma::vec &RHS, SimulationTime &simTime)
{
    ScopedTimer solveTimer(simTime.solve_time);
    return arma::solve(LHS, RHS, arma::solve_opts::fast);
}

// Transient Simulation
// Newton Raphson system solver for non-linear and dynamic elements
arma::vec NewtonRaphson_system(CKTcircuit &ckt, const double &h, const int &mode, const double time_trans, 
    const arma::vec &pre_global_solution, const Modelmap &modmap)
{
    int NR_iteration_counter = 0;
    bool isconverge = false;
    arma::vec solution = pre_global_solution;
    std::pair<arma::mat, arma::vec> init_matrices;
    std::pair<arma::mat, arma::vec> matrices;

    std::vector<arma::vec> NR_solutions(ITL4+1);
    NR_solutions[0] = pre_global_solution;

    init_matrices = Dynamic(ckt, h, pre_global_solution, mode, time_trans, ckt.sim_stats.simTime);

    for (int i = 1; i < 3; i++)
    {
        matrices = NonLinear(ckt, solution, init_matrices, modmap, h);
        // const arma::mat &LHS = matrices.first;
        // const arma::mat &RHS = matrices.second;

        // Solve Ax = b
        // J(v) * x(k+1) = [J(v)]x(k) - f(x(k))
        solution = solveDense(matrices.first, matrices.second, ckt.sim_stats.simTime);

        NR_iteration_counter += 1;
        NR_solutions.at(NR_iteration_counter) = solution;
        ckt.spiceCompatible.updateStateMachine(false);
    }

    isconverge = isConverge(NR_solutions, ckt, NR_iteration_counter);
    ckt.spiceCompatible.updateStateMachine(isconverge);

    while (!isconverge)
    {
        if (NR_iteration_counter >= ITL4 && mode == 1)
        {   
            // Transient simulation convergence failed
            converged = false;
            NR_ITE = NR_iteration_counter;
            ckt.sim_stats.total_NR_iteration += NR_iteration_counter;
            return solution;
        }
        else if (NR_iteration_counter >= 100 && mode == 0)
        {
            throw ConvergenceException("Transient Simulation did not converge at the op analysis!", "TRANSIENT_OP_CONVERGENCE");
        }
        else {
            matrices = NonLinear(ckt, solution, init_matrices, modmap, h);

            // Solve Ax = b
            // J(v) * x(k+1) = [J(v)]x(k) - f(x(k))
            solution = solveDense(matrices.first, matrices.second, ckt.sim_stats.simTime);

            NR_iteration_counter += 1;
            NR_solutions.at(NR_iteration_counter) = solution;

            isconverge = isConverge(NR_solutions, ckt, NR_iteration_counter);
            ckt.spiceCompatible.updateStateMachine(isconverge);
        }
    }

    converged = true;
    NR_ITE = NR_iteration_counter;
    ckt.sim_stats.total_NR_iteration += NR_iteration_counter;

    return solution;
}


// DC Analysis
arma::vec NewtonRaphson_system(CKTcircuit &ckt, const arma::mat &init_LHS, const arma::vec &init_RHS, const Modelmap &modmap)
{
    int NR_iteration_counter = 0;
    bool isconverge = false;
    arma::vec solution(init_RHS.n_rows, arma::fill::zeros);
    std::pair<arma::mat, arma::vec> init_matrices = {init_LHS, init_RHS};
    std::pair<arma::mat, arma::vec> matrices;

    std::vector<arma::vec> NR_solutions(ITL4+1);
    NR_solutions[0] = solution;

    // DC Analysis does not have dynamic elements (capacitors, inductors)!

    for (int i = 1; i < 3; i++)
    {
        matrices = NonLinear(ckt, solution, init_matrices, modmap, 0.0);

        // Solve Ax = b
        // J(v) * x(k+1) = [J(v)]x(k) - f(x(k))
        solution = solveDense(matrices.first, matrices.second, ckt.sim_stats.simTime);
        
        NR_iteration_counter += 1;
        NR_solutions.at(NR_iteration_counter) = solution;
        ckt.spiceCompatible.updateStateMachine(false);
    }

    isconverge = isConverge(NR_solutions, ckt, NR_iteration_counter);
    ckt.spiceCompatible.updateStateMachine(isconverge);

    while (!isconverge)
    {
        if (NR_iteration_counter >= 100)
        {
            throw ConvergenceException("DC Analysis did not converge!", "DC_CONVERGENCE");
        }

        matrices = NonLinear(ckt, solution, init_matrices, modmap, 0.0);

        // Solve Ax = b
        // J(v) * x(k+1) = [J(v)]x(k) - f(x(k))
        solution = solveDense(matrices.first, matrices.second, ckt.sim_stats.simTime);

        NR_iteration_counter += 1;
        NR_solutions.at(NR_iteration_counter) = solution;

        isconverge = isConverge(NR_solutions, ckt, NR_iteration_counter);
        ckt.spiceCompatible.updateStateMachine(isconverge);
    }

    NR_ITE = NR_iteration_counter;
    ckt.sim_stats.total_NR_iteration += NR_iteration_counter;

    return solution;
}
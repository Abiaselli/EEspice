#pragma once

#include <iostream>
#include <armadillo>
#include <cmath>
#include <chrono>
#include <string>
#include <vector>
#include <algorithm>
#include <array>

#include "sim_variables.hpp"
#include "Transient_calcs.hpp"
#include "integrate.hpp"
#include "eemath.hpp"

struct single_timestep
{
    CapacitanceState CapState; 

    double h{};             // time step
    double t{};             // time
    double next_h{};        // next time step (temp_h)

    arma::vec RHS;
    arma::mat LHS;
    arma::vec solution;

    bool TMAX_reach;
    bool TMIN_reach;
    bool ITL4_reach;

};

struct single_Truncation_error
{
    std::vector<double> LTE_bound{};
    std::vector<double> h_bound{};
};

// Old function for LTE calculation
// Calculate the truncation error (LTE)
// single_Truncation_error single_LTE_BE_calculation(const single_timestep &single_h, const std::vector<Transient> &vec_trans){
//     single_Truncation_error LTE;

//     arma::vec pre_current = vec_trans.back().C_current;
//     arma::vec pre_charge = vec_trans.back().C_charge;

//     // max(|ik+1|,|ik|) and max(|qk+1|,|qk|)
//     arma::vec max_current = arma::max(arma::abs(single_h.C_current), arma::abs(pre_current));
//     arma::vec max_charge = arma::max(arma::abs(single_h.C_charge), arma::abs(pre_charge));

//     // LTE current = ABSTOL + RELTOL * max(|ik+1|,|ik|) in SPICE book
//     arma::vec LTE_current = ABSTOL + RELTOL * max_current;

//     // LTE charge = RELTOL * max(|qk+1|,|qk|, chgtol) / h(k) in SPICE book
//     arma::vec CHGTOL_MAT = arma::vec(max_charge.n_rows, max_charge.n_cols).fill(CHGTOL);
//     arma::vec LTE_charge = RELTOL * arma::max(max_charge, CHGTOL_MAT) / single_h.h;

//     // LTE bound = max(LTE_current, LTE_charge)
//     LTE.LTE_bound = arma::max(LTE_current, LTE_charge);

//     // |i(k+1) - i(k)|
//     arma::vec BEcur_diff = single_h.C_current - pre_current;
    
//     // tol_h = 2C / |i(k+1) - i(k)| * LTE_bound
//     LTE.h_bound = ((2 * single_h.Capacitance) / arma::abs(BEcur_diff)) % LTE.LTE_bound;

//     return LTE;
// }

// LTE calculation from ngspice
single_Truncation_error single_LTE_ngspice(const single_timestep &single_h, const std::vector<Transient> &vec_trans, const CKTcircuit &ckt)
{
    const auto &prev          = vec_trans.back();
    const std::size_t history_sz = vec_trans.size();
    const unsigned order      = ckt.CKTorder;
    const unsigned DIVorder   = order + 1;

    // need DIVorder previous points in history
    if (history_sz < DIVorder) {
        throw std::invalid_argument(
            "Not enough history: need at least " 
            + std::to_string(DIVorder) + " previous steps");
    }

    const std::size_t N = prev.CapState.CapCurrent.size();
    single_Truncation_error LTE;
    LTE.LTE_bound.reserve(N);
    LTE.h_bound.reserve(N);

    // coefficients from SPICE
    static constexpr std::array<double,6> gearCoeff = {
        0.5, 0.22222222, 0.13636364, 0.096,
        0.07299270, 0.05830904
    };
    static constexpr std::array<double,2> trapCoeff = {
        0.5,       // BE (order=1)
        0.0833333  // Trap (order=2)
    };

    double coeff;
    switch (ckt.CKTintegrateMethod) {
      case BACKWARD_EULER:
        coeff = trapCoeff[0];
        break;
      case TRAPEZOIDAL:
        if (order < 1 || order > trapCoeff.size())
            throw std::out_of_range("Unsupported TRAP order " + std::to_string(order));
        coeff = trapCoeff[order - 1];
        break;
      case GEAR:
        if (order < 1 || order > gearCoeff.size())
            throw std::out_of_range("Unsupported GEAR order " + std::to_string(order));
        coeff = gearCoeff[order - 1];
        break;
      default:
        throw std::invalid_argument("Unknown integration method");
    }

    // index of the oldest history point to use
    const std::size_t start_idx = history_sz - DIVorder;

    for (std::size_t i = 0; i < N; ++i) {
        // 1) voltage tolerance
        double volttol = ABSTOL
                       + RELTOL * std::max(
                            std::abs(single_h.CapState.CapCurrent[i]),
                            std::abs(prev.CapState.CapCurrent[i]));

        // 2) charge tolerance
        double chargetol = std::max(
            std::abs(single_h.CapState.CapCharge[i]),
            std::abs(prev.CapState.CapCharge[i]));
        chargetol = RELTOL *
            std::max(chargetol, CHGTOL) / single_h.h;

        // 3) combined tolerance bound
        double tol = std::max(volttol, chargetol);
        LTE.LTE_bound.push_back(tol);

        // 4) collect charge history and time steps
        std::vector<double> qvals(DIVorder + 1);
        std::vector<double> deltas(DIVorder);
        for (unsigned j = 0; j < DIVorder; ++j) {
            qvals[j]  = vec_trans[start_idx + j].CapState.CapCharge[i];
            deltas[j] = vec_trans[start_idx + j + 1].h;
        }
        qvals[DIVorder]        = single_h.CapState.CapCharge[i];
        deltas[DIVorder - 1]   = single_h.h;

        // 5) in-place divided differences to get (order+1)th diff
        for (int level = DIVorder; level > 0; --level) {
            for (int k = 0; k < level; ++k) {
                qvals[k] = (qvals[k] - qvals[k + 1]) / deltas[k];
            }
            if (level > 1) {
                for (int k = 0; k < level - 1; ++k) {
                    deltas[k] = deltas[k + 1] + deltas[k];
                }
            }
        }
        double dd = qvals[0];

        // 6) compute new h bound
        double ratio = (TRTOL * tol)
                     / std::max(ABSTOL, coeff * std::abs(dd));
        double hbound;
        if (order == 2) {
            hbound = std::sqrt(ratio);
        } else if (order > 2) {
            hbound = std::exp(std::log(ratio) / order);
        } else {
            hbound = ratio;
        }
        LTE.h_bound.push_back(hbound);
    }

    return LTE;
}

// New divided differences method for the LTE calculation
single_Truncation_error single_LTE_divided_diff(const single_timestep &single_h, const std::vector<Transient> &vec_trans, const CKTcircuit &ckt){
    const auto& prev        = vec_trans.back();
    const auto history_sz   = vec_trans.size();
    const unsigned CKTorder = ckt.CKTorder;
    const unsigned DIVorder = CKTorder + 1;      // input order to DividedDiff
    const auto fact         = Math::factorial(DIVorder);
    const std::size_t N     = prev.CapState.CapCurrent.size(); // number of capacitors and bsim4 capacitors
    single_Truncation_error LTE;
    
    // --- sanity check: need at least 'DIVorder' previous steps in vec_trans
    if (history_sz < DIVorder) {
        throw std::invalid_argument(
            "Not enough history to compute "
            + std::to_string(DIVorder) + "th DIVorder divided difference");
    }

    // --- choose the error‐coefficient C_{CKTorder+1}
    static constexpr std::array<double, 6> gearCoeff = {
        0.5,        
        0.22222222, 
        0.13636364, 
        0.096,      
        0.07299270, 
        0.05830904  
    };
    static constexpr std::array<double, 2> trapCoeff = {
        0.5,       // Trapezoid with CKTorder=1 (degenerates to BE)
        0.0833333  // Trapezoid with CKTorder=2
    };

    double coeff;
    switch (ckt.CKTintegrateMethod) {
      case BACKWARD_EULER:
        coeff = trapCoeff[0]; // first order of Trapezoidal method is BE
        break;
      case TRAPEZOIDAL:
        if (CKTorder == 0 || CKTorder > trapCoeff.size())
            throw std::invalid_argument("Unsupported Trap DIVorder " + std::to_string(CKTorder));
        coeff = trapCoeff[CKTorder - 1];
        break;
      case GEAR:
        if (CKTorder == 0 || CKTorder > gearCoeff.size())
            throw std::invalid_argument("Unsupported Gear DIVorder " + std::to_string(CKTorder));
        coeff = gearCoeff[CKTorder - 1];
        break;
      default:
        throw std::invalid_argument("Unknown integration method");
    }
    const double denomFactor = fact * coeff;

    LTE.LTE_bound.reserve(N);
    LTE.h_bound .reserve(N);

    // --- work buffers (reuse per element)
    std::vector<double> qValues(DIVorder + 1);
    std::vector<double> hValues(DIVorder);

    // starting index in vec_trans for the oldest time‐point
    const std::size_t start_idx = history_sz - DIVorder;

    for (std::size_t i = 0; i < N; ++i) {
        // 1) compute the per-element LTE bound
        // LTE current = ABSTOL + RELTOL * max(|ik+1|,|ik|) in SPICE book
        double volttol  = ABSTOL + RELTOL * std::max(std::abs(single_h.CapState.CapCurrent[i]), std::abs(prev.CapState.CapCurrent[i]));
        // LTE charge = RELTOL * max(|qk+1|,|qk|, chgtol) / h(k) in SPICE book
        double chargetol = RELTOL * std::max({std::abs(single_h.CapState.CapCharge[i]), std::abs(prev.CapState.CapCharge[i]), CHGTOL}) / single_h.h;
        // LTE bound = max(LTE_current, LTE_charge)
        double bound = std::max(volttol, chargetol);
        LTE.LTE_bound.push_back(bound);

        // 2) build the divided‐difference inputs:
        //    time points t_{k-CKTorder}, …, t_k (in vec_trans), then t_{k+1} (in single_h)
        for (unsigned j = 0; j < DIVorder; ++j) {
            // qValues[0..DIVorder-1] from history
            qValues[j] = vec_trans[start_idx + j].CapState.CapCharge[i];
            // hValues[0..DIVorder-2] from history steps
            hValues[j] = vec_trans[start_idx + j + 1].h;
        }
        // last entries
        qValues[DIVorder]   = single_h.CapState.CapCharge[i];
        hValues[DIVorder-1] = single_h.h;

        // 3) compute the (CKTorder+1)th divided difference DD_{CKTorder+1}[q]
        double dd = Math::DividedDiff(qValues, hValues, DIVorder);

        // 4) compute the h-bound via
        //       h_bound = [ (TRTOL * bound) / max(ABSTOL, (CKTorder+1)!·C_{CKTorder+1}·|dd| ) ]^(1/(CKTorder+1))
        double LocalTruncationError = std::max(ABSTOL, denomFactor * std::abs(dd));
        if(CKTorder == 1){
            LTE.h_bound.push_back(std::sqrt((TRTOL * bound) / LocalTruncationError));
        }
        else{
            double h_bound = std::pow((TRTOL * bound) / LocalTruncationError, 1.0 / (CKTorder + 1));
            LTE.h_bound.push_back(h_bound);
        }
    }

    return LTE;
}

single_timestep single_solution_solver(const double &h, const Transient &trans, CKTcircuit &ckt, const TransientSimulator &trans_sim,
                                        const Modelmap &modmap){

    single_timestep single_h;
    single_h.h = h;
    single_h.t = trans_sim.vec_trans.back().time_trans + h;


    if(trans_sim.trans_config.non_linear){
        single_h.solution = NewtonRaphson_system(ckt, h, trans.mode, single_h.t, trans_sim.vec_trans.back().solution, modmap);
        if(NR_ITE < ITL4){
            single_h.CapState = get_cap_state(ckt, single_h.solution, single_h.h, trans_sim.vec_trans);
        }
    }
    else{
        // Dynamic() resets to linear baseline, stamps capacitors, saves step baseline
        // It works directly on ckt.cktmatrix
        Dynamic(ckt, h, trans_sim.vec_trans.back().solution, trans.mode, single_h.t, ckt.sim_stats.simTime);
        single_h.solution = solver(ckt.cktmatrix->LHS, ckt.cktmatrix->RHS, ckt);
        single_h.CapState = get_cap_state(ckt, single_h.solution, single_h.h, trans_sim.vec_trans);
    }

    return single_h;
}

bool single_LTE_check(single_Truncation_error &LTE, const single_timestep &single_h,
                        double &temp_h, const TransientSimulator &trans_sim, const CKTcircuit &ckt){
    if(converged){

       LTE = single_LTE_ngspice(single_h, trans_sim.vec_trans, ckt);

       if(LTE.h_bound.empty()){
           throw SimulationException("Error in single_LTE_check function: h_bound is empty!", "SINGLE_LTE_CHECK");
       }
       else[[likely]]{
            auto iter_h_min = std::ranges::min_element(LTE.h_bound);
            double h_min = *iter_h_min;
            if(h_min < 0.9 * single_h.h){
                // Reject the solution, hn = hn+1 and recompute the new hn
                temp_h = h_min;
                return false;
            }
            else{
                // Accept the solution
                temp_h = std::min({h_min, 2.0 * single_h.h, trans_sim.trans_config.h_MAX});
                return true;
            }
       }

    }
    else{
        // Reject the solution, hn = hn/8 and recompute the new hn
        temp_h = single_h.h / 8.0;

        if(temp_h > trans_sim.trans_config.h_MIN){
            return false;
        }
        else{
            throw SimulationException("The time step is too small (Reach minimum!)", "SINGLE_LTE_CHECK");
        }
    }
}

single_timestep single_next_h(const Transient &trans, CKTcircuit &ckt, const TransientSimulator &trans_sim, const Modelmap &modmap){

    double temp_h = trans_sim.vec_trans.back().next_h;                     // The temporary time step for this function

    bool LTE_check= false;

    single_timestep single_h;
    single_Truncation_error LTE;

    do{
        NIcomCof(ckt, temp_h); // Calculate the timestep-dependent terms used in the numerical integration
        single_h = single_solution_solver(temp_h, trans, ckt, trans_sim, modmap);

        total_timepoint += 1;

        LTE_check = single_LTE_check(LTE, single_h, temp_h, trans_sim, ckt);

    }while(LTE_check == false);

    // Sometimes temp_h(next_h) sometimes is 2 * single_h.h!
    single_h.next_h = temp_h;

    return  single_h;
}
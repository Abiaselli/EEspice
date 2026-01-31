#pragma once
#include <cmath>
#include <armadillo>
#include <iostream>
#include "CKT.hpp"
#include "Transient_calcs.hpp"


// Get the voltage across the capacitor
double get_cap_vol(const arma::vec &solution, int node_x, int node_y){
    double vol = 0.0;
    if (node_x == 0)
    {
        vol = solution(node_y - 1);
    }
    else if (node_y == 0)
    {
        vol = solution(node_x - 1);
    }
    else
    {
        vol = solution(node_x - 1) - solution(node_y - 1);
    }
    return vol;
}

// Get the Charge of the capacitor and the current of the capacitor
CapacitanceState get_cap_state(const CKTcircuit &ckt, const arma::vec &solution, const double h, const std::vector<Transient> &vec_trans){
    CapacitanceState CapState;
    CapState.CapCharge.reserve(ckt.num_of_states);
    CapState.CapCurrent.reserve(ckt.num_of_states);
    auto &CapCharge = CapState.CapCharge;
    auto &CapCurrent = CapState.CapCurrent;

    double pre_charge{}, charge{}, current{};

    if(vec_trans.empty()){
        // OP analysis
        CapState.CapCharge.resize(ckt.num_of_states, 0.0);
        CapState.CapCurrent.resize(ckt.num_of_states, 0.0);
        return CapState;
    }

    // Update the charge and current of the capacitors
    for (int i = 0; i < ckt.CKTelements.capacitors.size(); ++i){
        const auto &cap = ckt.CKTelements.capacitors[i];
        double vol = get_cap_vol(solution, cap.nodePos, cap.nodeNeg);
        pre_charge = vec_trans.back().CapState.CapCharge[i];
        charge = cap.value * vol;
        switch(ckt.CKTintegrateMethod) {
            // i = c/h * (u(k+1) - u(k))
            case BACKWARD_EULER: 
                current = ckt.CKTag[0] * charge + ckt.CKTag[1] * pre_charge;
                break;
            case TRAPEZOIDAL:
                switch (ckt.CKTorder) {
                    case 1:
                        // order 1 is the backward euler
                        current = ckt.CKTag[0] * charge + ckt.CKTag[1] * pre_charge;
                        break;
                    
                    default:
                        throw SimulationException("Illegal integration order", "get_cap_state");
                }
                break;
            case GEAR:
                throw SimulationException("GEAR method not implemented yet.", "get_cap_state");
                break;
            default:
                throw SimulationException("Unknown integration method", "get_cap_state");
                break;
        }
        CapCharge.push_back(std::abs(charge));
        CapCurrent.push_back(std::abs(current));
    }
    // Update the charge and current of the bsim4 capacitors
    // They are calculated in the bsim4load function by NIintegrate
    for (const auto &bsim4 : ckt.CKTelements.bsim4)
    {
        CapCharge.push_back(std::abs(bsim4.bsim4v82Instance.BSIM4states0[bsim4::BSIM4qb])); // charge update
        CapCurrent.push_back(std::abs(bsim4.bsim4v82Instance.BSIM4states0[bsim4::BSIM4qb+1])); // ccap = (qcap+1)
        CapCharge.push_back(std::abs(bsim4.bsim4v82Instance.BSIM4states0[bsim4::BSIM4qg]));
        CapCurrent.push_back(std::abs(bsim4.bsim4v82Instance.BSIM4states0[bsim4::BSIM4qg+1]));
        CapCharge.push_back(std::abs(bsim4.bsim4v82Instance.BSIM4states0[bsim4::BSIM4qd]));
        CapCurrent.push_back(std::abs(bsim4.bsim4v82Instance.BSIM4states0[bsim4::BSIM4qd+1]));
        if (bsim4.bsim4v82Instance.BSIM4trnqsMod){
            CapCharge.push_back(std::abs(bsim4.bsim4v82Instance.BSIM4states0[bsim4::BSIM4qcdump]));
            CapCurrent.push_back(std::abs(bsim4.bsim4v82Instance.BSIM4states0[bsim4::BSIM4qcdump+1]));
        }
        if(bsim4.bsim4v82Instance.BSIM4rbodyMod){
            CapCharge.push_back(std::abs(bsim4.bsim4v82Instance.BSIM4states0[bsim4::BSIM4qbs]));
            CapCurrent.push_back(std::abs(bsim4.bsim4v82Instance.BSIM4states0[bsim4::BSIM4qbs+1]));
        }
        if(bsim4.bsim4v82Instance.BSIM4rgateMod == 3){
            CapCharge.push_back(std::abs(bsim4.bsim4v82Instance.BSIM4states0[bsim4::BSIM4qgmid]));
            CapCurrent.push_back(std::abs(bsim4.bsim4v82Instance.BSIM4states0[bsim4::BSIM4qgmid+1]));
        }
    }
    
    return CapState;
}

// nicomcof.c
/**********
Copyright 1990 Regents of the University of California.  All rights reserved.
Author: 1985 Thomas L. Quarles
**********/

/* xmu=0:    Backward Euler
 * xmu=0.5:  trapezoidal (standard)
 */

void
NIcomCof(CKTcircuit &ckt, double h)
{
    // double mat[8][8];   /* matrix to compute the gear coefficients in */
    // int i,j,k;          /* generic loop indicies */
    // double arg;
    // double arg1;
    ScopedTimer t(ckt.sim_stats.simTime.nicomcof_time);

    /*  this routine calculates the timestep-dependent terms used in the
     *  numerical integration.
     */ 

    /*  
     *  compute coefficients for particular integration method 
     */ 
    double xmu;
    switch(ckt.CKTintegrateMethod) {
        case BACKWARD_EULER:
            ckt.CKTag[0] = 1 / h;
            ckt.CKTag[1] = -1 / h;
            break;

        case TRAPEZOIDAL:
            switch(ckt.CKTorder) {

                case 1:
                    ckt.CKTag[0] = 1 / h;
                    ckt.CKTag[1] = -1 / h;
                    break;

                case 2:
                    xmu = 0.5;
                    ckt.CKTag[0] = 1.0 / h / (1.0 - xmu);
                    ckt.CKTag[1] = xmu / (1.0 - xmu);
                    break;

                default:
                    throw SimulationException("Trapezoidal method only supports order 1 or 2.", "NIcomCof");
            }
            break;

        case GEAR:
            throw SimulationException("GEAR method not implemented yet.", "NIcomCof");
            break;

    }
}
#pragma once
#include "CKT.hpp"
#include "XB_timer.hpp"
#include "simulation_exceptions.hpp"

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

// niinteg.c
/**********
Copyright 1990 Regents of the University of California.  All rights reserved.
Author: 1985 Thomas L. Quarles
**********/
/* NIintegrate(ckt,geq,ceq,cap,qcap)
*  integrate the specified capacitor - method and order in the
*  ckt structure, ccap follows qcap.
*/
#pragma once
#include <iostream>
#include <string>
#include "bsim4v82.hpp"
#include "CKT.hpp"

namespace bsim4{
 
#define ccap qcap+1

int
NIintegrate(const CKTcircuit &ckt, BSIM4V82 &inst, double *geq, double *ceq, double cap, int qcap, double h)
{
    const std::string ordmsg = "Illegal integration order";
    const std::string methodmsg = "Unknown integration method";

    switch(ckt.CKTintegrateMethod) {

    case TRAPEZOIDAL:
        switch(ckt.CKTorder) {
        case 1:
            inst.BSIM4states0[ccap] = ckt.CKTag[0] * inst.BSIM4states0[qcap]
                    + ckt.CKTag[1] * inst.BSIM4states1[qcap];
            break;
        case 2:
            inst.BSIM4states0[ccap] = - inst.BSIM4states1[ccap] * ckt.CKTag[1] +
                    ckt.CKTag[0] *
                    ( inst.BSIM4states0[qcap] - inst.BSIM4states1[qcap] );
            break;
        default:
            std::cerr << ordmsg << std::endl;
            exit(1);
        }
        break;
    // case GEAR:
    //     inst.BSIM4states0[ccap]=0;
    //     switch(ckt.CKTorder) {

    //     case 6:
    //         inst.BSIM4states0[ccap] += ckt.CKTag[6]* ckt.CKTstate6[qcap];
    //         /* fall through */
    //     case 5:
    //         inst.BSIM4states0[ccap] += ckt.CKTag[5]* ckt.CKTstate5[qcap];
    //         /* fall through */
    //     case 4:
    //         inst.BSIM4states0[ccap] += ckt.CKTag[4]* ckt.CKTstate4[qcap];
    //         /* fall through */
    //     case 3:
    //         inst.BSIM4states0[ccap] += ckt.CKTag[3]* ckt.CKTstate3[qcap];
    //         /* fall through */
    //     case 2:
    //         inst.BSIM4states0[ccap] += ckt.CKTag[2]* ckt.CKTstate2[qcap];
    //         /* fall through */
    //     case 1:
    //         inst.BSIM4states0[ccap] += ckt.CKTag[1]* inst.BSIM4states1[qcap];
    //         inst.BSIM4states0[ccap] += ckt.CKTag[0]* inst.BSIM4states0[qcap];
    //         break;

    //     default:
    //         return(E_ORDER);

    //     }
    //     break;

    default:
        std::cerr << methodmsg << std::endl;
        exit(1);
    }
    *ceq = inst.BSIM4states0[ccap] - ckt.CKTag[0] * inst.BSIM4states0[qcap];
    *geq = ckt.CKTag[0] * cap;
    return 0;
}

} // namespace bsim4
 
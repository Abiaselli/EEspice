/*
   b4ld.c - BSIM4v4.8.2
*/
/* ******************************************************************************
   *  BSIM4 4.8.2 released by Chetan Kumar Dabhi 01/01/2020                     *
   *  BSIM4 Model Equations                                                     *
   ******************************************************************************

   ******************************************************************************
   *  Copyright (c) 2020 University of California                               *
   *                                                                            *
   *  Project Director: Prof. Chenming Hu.                                      *
   *  Current developers: Chetan Kumar Dabhi   (Ph.D. student, IIT Kanpur)      *
   *                      Prof. Yogesh Chauhan (IIT Kanpur)                     *
   *                      Dr. Pragya Kushwaha  (Postdoc, UC Berkeley)           *
   *                      Dr. Avirup Dasgupta  (Postdoc, UC Berkeley)           *
   *                      Ming-Yen Kao         (Ph.D. student, UC Berkeley)     *
   *  Authors: Gary W. Ng, Weidong Liu, Xuemei Xi, Mohan Dunga, Wenwei Yang     *
   *           Ali Niknejad, Chetan Kumar Dabhi, Yogesh Singh Chauhan,          *
   *           Sayeef Salahuddin, Chenming Hu                                   *
   ******************************************************************************/
#pragma once
#include "sim_variables.hpp"
#include "SPICEcompatible.hpp"
#include "devsup.hpp"
#include "bsim4v82/bsim4v82.hpp"
#include "bsim4v82/bsim4v82const.hpp"
#include "bsim4v82/bsim4v82NI.hpp"
#include "bsim4v82calculateStamps.hpp"
#include "bsim4v82applyStamps.hpp"
#include "hybrid_matrix.hpp"

#include <cmath>
#include <armadillo>

namespace bsim4{

/**
 * @brief Main loading function for the BSIM4 model.
 *
 * This function orchestrates the calculation of device physics and the subsequent loading
 * (stamping) of the results into the circuit's system matrices (LHS and RHS).
 * It calls BSIM4calculateStamps to compute all necessary values and then
 * bsim4applyStamps to perform the matrix modifications.
 * @return An integer error code (0 for success).
 */
int
BSIM4load(const CKTcircuit &ckt, const BSIM4model &model, BSIM4V82 &instance, const SPICECompatible &spice, const arma::vec &presolution,
    const double CKTtemp, const double CKTgmin, HybridMatrix &LHS, arma::vec &RHS, const BSIM4StampIndexCache *cache)
{
    // Calculate all the values needed for the matrices.
    const auto stamps = BSIM4calculateStamps(ckt, model, instance, spice, presolution, CKTtemp, CKTgmin);

    // Apply the calculated values to the LHS and RHS matrices.
    if (cache && cache->built && LHS.is_sparse() && LHS.is_pattern_locked()) {
        bsim4applyStampsCached(instance, stamps, *cache, LHS, RHS);
    } else {
        bsim4applyStamps(instance, stamps, LHS, RHS);
    }

    return 0; // success
}

int
BSIM4load(const CKTcircuit &ckt, const BSIM4model &model, BSIM4V82 &instance, const SPICECompatible &spice, const arma::vec &presolution,
    const double CKTtemp, const double CKTgmin, HybridMatrix &LHS, arma::vec &RHS)
{   
    return BSIM4load(ckt, model, instance, spice, presolution, CKTtemp, CKTgmin, LHS, RHS, nullptr);
}

void updateState1(BSIM4V82 &inst){
    inst.BSIM4states1 = inst.BSIM4states0;
}


} // namespace bsim4

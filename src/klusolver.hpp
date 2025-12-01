#pragma once
#include <iostream>
#include <klu.h>
#include "CKT.hpp"

// KLU Solver interface for EEspice

int KLUpreOrder(KLUmatrix& kluMatrix, const std::int64_t n)
{
    kluMatrix.KLUmatrixSymbolic = klu_l_analyze(n, kluMatrix.KLUmatrixAp, kluMatrix.KLUmatrixAi, &kluMatrix.KLUmatrixCommon);
    if (!kluMatrix.KLUmatrixSymbolic) {
        throw MatrixException("Error (PreOrder): KLUsymbolic object is nullptr. A problem occurred", "KLUpreOrder");
    }
    return 0;
}

int KLUreOrder(KLUmatrix& kluMatrix, double PivRel)
{
    kluMatrix.KLUmatrixCommon.tol = PivRel;

    // We must free the old Numeric object before creating a new one.
    if (kluMatrix.KLUmatrixNumeric) {
        klu_l_free_numeric(&kluMatrix.KLUmatrixNumeric, &kluMatrix.KLUmatrixCommon);
    }

    kluMatrix.KLUmatrixNumeric = klu_l_factor(  kluMatrix.KLUmatrixAp, kluMatrix.KLUmatrixAi, 
                                                kluMatrix.KLUmatrixAx, kluMatrix.KLUmatrixSymbolic, 
                                                &kluMatrix.KLUmatrixCommon);
    if (!kluMatrix.KLUmatrixNumeric) {
        if (kluMatrix.KLUmatrixCommon.status == KLU_SINGULAR){
            throw MatrixException("Warning (Factor): KLU Matrix is SINGULAR", "KLUreOrder");
            throw MatrixException(" Numerical Rank: " + std::to_string(kluMatrix.KLUmatrixCommon.numerical_rank), "KLUreOrder");
            throw MatrixException(" Singular Node: " + std::to_string(kluMatrix.KLUmatrixCommon.singular_col + 1), "KLUreOrder");
        }
        if (!kluMatrix.KLUmatrixSymbolic){
            throw MatrixException("Error (Factor): KLUsymbolic object is nullptr. A problem occurred", "KLUreOrder");
        }
        throw MatrixException("Error (Factor): KLUnumeric object is nullptr. A problem occurred", "KLUreOrder");
    }
    return 0;
}

int KLUluFac(KLUmatrix& kluMatrix, double PivRel){
    int success = klu_l_refactor(kluMatrix.KLUmatrixAp, kluMatrix.KLUmatrixAi, 
                                kluMatrix.KLUmatrixAx, kluMatrix.KLUmatrixSymbolic, 
                                kluMatrix.KLUmatrixNumeric, &kluMatrix.KLUmatrixCommon);
    if (!success) {
        // Refactor failed (pivots unstable due to value changes).
        // Note: Common.status will likely contain KLU_SINGULAR or similar.
        KLUreOrder(kluMatrix, PivRel);
    }
    return 0;
}

int KLUsolver(arma::sp_mat& LHS, arma::vec& RHS, KLUmatrix& kluMatrix, CKTcircuit& ckt) 
{
    // --- ZERO-COPY KLU INTEGRATION ---
    // Armadillo exposes 'const uword*'. KLU 'l' functions expect 'SuiteSparse_long*' or 'std::int64_t*'.
    // We use const_cast to strip const (KLU doesn't modify Ap/Ai, but signature requires it)
    // We use reinterpret_cast to treat unsigned long as signed long.

    // 1. Cast Pointers
    kluMatrix.KLUmatrixAp = reinterpret_cast<std::int64_t*>(const_cast<arma::uword*>(LHS.col_ptrs));
    kluMatrix.KLUmatrixAi = reinterpret_cast<std::int64_t*>(const_cast<arma::uword*>(LHS.row_indices));

    // 2. Analyze (Use klu_l_analyze) 
    // Symbolic analyze happens once per matrix
    if (!ckt.NIDIDPREORDER){
        KLUpreOrder(kluMatrix, LHS.n_rows);
        ckt.NIDIDPREORDER = true;
    }

    // 3. Factor (Use klu_l_factor)
    if (ckt.NISHOULDREORDER){
        KLUreOrder(kluMatrix, ckt.CKTpivotRelTol);
        ckt.NISHOULDREORDER = false;
    }
    else {
        // Try Refactoring (Use klu_l_refactor)
        KLUluFac(kluMatrix, ckt.CKTpivotRelTol);
    }

    // 4. Solve (Use klu_l_solve)
    int solve_status = klu_l_solve(kluMatrix.KLUmatrixSymbolic, kluMatrix.KLUmatrixNumeric, 
                                    LHS.n_rows, 1, RHS.memptr(), &kluMatrix.KLUmatrixCommon);
    if (!solve_status){
        if(kluMatrix.KLUmatrixCommon.status == KLU_SINGULAR){
            throw MatrixException("Warning (Solve): KLU Solve failed due to SINGULAR matrix", "KLUsolver");
        }
        throw MatrixException("Error (Solve): KLU Solve failed. A problem occurred", "KLUsolver");
    }

    return 0;
}
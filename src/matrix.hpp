#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <variant>
#include <armadillo>
#include "hybrid_matrix.hpp"

struct Matrix
{
    int Maxi{}; // Size of matrix without branch current
    int Maxj{};

    HybridMatrix LHS; // LHS matrix (can be dense or sparse)
    arma::vec RHS; // RHS matrix

    arma::cx_dmat LHS_cx; // Complex LHS matrix for AC analysis
    arma::cx_dvec RHS_cx; // Complex RHS matrix for AC analysis

    int n_rows{}; // Number of rows and columns
    int n_cols{};
    bool use_sparse{false}; // Flag to determine if sparse matrices should be used

    // =========================================================================
    // Two-Level RHS Baseline Support for Newton-Raphson Reset Optimization
    // =========================================================================

    /** Save current RHS as linear baseline */
    void save_linear_baseline_RHS() {
        baseline_linear_RHS = RHS;
    }

    /** Reset RHS to linear baseline */
    void reset_to_linear_baseline_RHS() {
        if (baseline_linear_RHS.n_elem > 0) {
            RHS = baseline_linear_RHS;
        }
    }

    /** Save current RHS as step baseline */
    void save_step_baseline_RHS() {
        baseline_step_RHS = RHS;
    }

    /** Reset RHS to step baseline */
    void reset_to_step_baseline_RHS() {
        if (baseline_step_RHS.n_elem > 0) {
            RHS = baseline_step_RHS;
        }
    }

    /** Initialize step baseline from linear baseline */
    void copy_linear_to_step_baseline_RHS() {
        baseline_step_RHS = baseline_linear_RHS;
    }

    // =========================================================================
    // AC Analysis Support (unchanged)
    // =========================================================================

    void set_init_cxmatrix()
    {
        init_cxLHS = LHS_cx;
        init_cxRHS = RHS_cx;
    }
    arma::cx_dmat get_init_cxLHS() const
    {
        return init_cxLHS;
    }
    arma::cx_dvec get_init_cxRHS() const
    {
        return init_cxRHS;
    }

private:
    // RHS baselines for two-level reset
    arma::vec baseline_linear_RHS;
    arma::vec baseline_step_RHS;

    // Complex matrices for AC analysis
    arma::cx_dmat init_cxLHS;
    arma::cx_dvec init_cxRHS;
};

struct MNA {
    HybridMatrix LHS;
    arma::vec RHS;
};

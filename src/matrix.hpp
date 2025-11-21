#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <variant>
#include <armadillo>
#include "hybrid_matrix.hpp"

struct DenseMatrix;

struct DenseMatrix
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

    void set_initmatrix()
    {
        init_LHS = LHS;
        init_RHS = RHS;
        n_rows = LHS.rows();
        n_cols = LHS.cols();
    }
    void set_init_cxmatrix()
    {
        init_cxLHS = LHS_cx;
        init_cxRHS = RHS_cx;
    }
    HybridMatrix get_init_LHS() const
    {
        return init_LHS;
    }
    arma::vec get_init_RHS() const
    {
        return init_RHS;
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
    HybridMatrix init_LHS; // The initial LHS and RHS values to be used in the NR-algorithm
    arma::vec init_RHS;
    arma::cx_dmat init_cxLHS;   // The initial complex LHS and RHS values to be used in the NR-algorithm for AC analysis
    arma::cx_dvec init_cxRHS;
};

struct MNA {
    HybridMatrix LHS;
    arma::vec RHS;
};

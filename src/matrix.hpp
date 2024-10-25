#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <variant>
#include <armadillo>

struct DenseMatrix
{
    int Maxi{}; // Size of matrix without branch current
    int Maxj{};

    arma::mat LHS; // LHS matrix
    arma::vec RHS; // RHS matrix

    int n_rows{}; // Number of rows and columns
    int n_cols{};

    void set_initmatrix()
    {
        init_LHS = LHS;
        init_RHS = RHS;
        n_rows = LHS.n_rows;
        n_cols = LHS.n_cols;
    }
    arma::mat get_init_LHS() const
    {
        return init_LHS;
    }
    arma::vec get_init_RHS() const
    {
        return init_RHS;
    }

private:
    arma::mat init_LHS; // The initial LHS and RHS values to be used in the NR-algorithm
    arma::vec init_RHS;
};

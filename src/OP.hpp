#pragma once
#include <iostream>
#include <armadillo>
#include "Newton.hpp"
#include "CKT.hpp"
#include "map.hpp"
#include "SPICEcompatible.hpp"


// Function to perform operating point analysis
arma::vec OperatingPointAnalysis(CKTcircuit &ckt, const Modelmap &modmap, bool non_linear){
    // Set the flags for SPICE compatibility
    ckt.spiceCompatible.setFlagsOP();

    // Get the initial LHS and RHS matrices
    arma::mat init_LHS = ckt.cktdematrix->get_init_LHS();
    arma::vec init_RHS = ckt.cktdematrix->get_init_RHS();

    // Solve the linear system to get the operating point
    if(non_linear){
        // If the circuit is non-linear, use the Newton-Raphson method
        return NewtonRaphson_system(ckt, init_LHS, init_RHS, modmap);
    } else {
        // If the circuit is linear, use the direct solver
        return arma::solve(init_LHS, init_RHS);
    }

    // Todo: Add gmin stepping and stepping sources...
}

void printOperatingPoint(const arma::vec &op_solution, const CKTcircuit &ckt){
    // Print the operating point solution
    arma::vec node_volt_print = op_solution.submat(0, 0, ckt.external_nodes - 1, 0);
    arma::vec current_print = op_solution.submat(op_solution.n_rows - ckt.no_of_V_sources, 0, op_solution.n_rows - 1, 0);
    arma::vec op_solution_print = arma::join_vert(node_volt_print, current_print);
    op_solution_print.print("Operating Point Solution: ");
}
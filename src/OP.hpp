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
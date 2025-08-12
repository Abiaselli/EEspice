#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <variant>
#include <map>
#include <memory>
#include <array>
#include "SPICEcompatible.hpp"


struct CKTcircuit
{
    int pulse_num{};                                // Number of pulse voltages
    // uword is a typedef for an unsigned integer type; it is used for matrix indices as well as all internal counters and loops

    CircuitElements CKTelements;                    // Vector of circuit elements
    int external_nodes{};                           // Number of external nodes excluding ground and nodes inside the MOSFETs
    int internal_nodes{};                           // Number of internal nodes (nodes inside the MOSFETs)
    int num_of_states{};                            // Number of dynamic states (capacitors and bsim4v82)
    int T_nodes{};                                  // Total number of nodes excluding ground (Used to create the initial matrix in CKTsetup)

    std::shared_ptr<DenseMatrix> cktdematrix;       // Dense matrix struct
    Circuitmap map;                                 // Circuit and Model Map struct
    bool ckt_loaded{};                              // To check if the circuit is loaded or not

    double CKTtemp{};                               // Actual temperature of CKT, initialzed to 300.15 K 
    double CKTnomTemp = 300.15;               // Reference temperature 300.15 K
    double CKTgmin = 1.0e-12;                       // Gmin value
    int CKTintegrateMethod{};                       // Integration method (0 for Backward Euler, 1 for Trapezoidal, 2 for Gear)
    int CKTorder{};                                 // Order of the integration method (1 for first order, 2 for second order)
    std::array<double, 7> CKTag;                    // Coefficients for the integration method, 1/h and -1/h for BE

    SPICECompatible spiceCompatible;                // SPICE-compatible (cktmode)
};

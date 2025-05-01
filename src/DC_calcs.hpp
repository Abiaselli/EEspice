#pragma once
#include <vector>
#include <string>
#include <armadillo>
#include <variant>
#include <functional>
#include "CKT.hpp"
#include "device.hpp"
#include "global.hpp"
#include "Newton.hpp"
#include "circuit_parser.hpp"


struct DCSweepSpec {
    // Using in the parser
    std::string sourceName;
    double vstart{};
    double vend{};
    double vstep{};
    std::vector<double> sweep_values;   // All sweep values from vstart to vend
};

// A structure for multi-sweep results
struct DC{
    arma::mat LHS;
    arma::vec RHS;
    double sweepValue;        //  source values
    std::string sweepName;    //  source names
    arma::vec solution;
};

struct DCSimulator {
    DCSweepSpec dcsweep;                   // sweep device data
    std::vector<DC> vec_dc;                // All DC results
    bool non_linear = false;
};

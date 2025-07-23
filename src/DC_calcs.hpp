#pragma once
#include <vector>
#include <string>
#include <armadillo>

namespace dc{
struct DCSweepSpec {
    // Using in the parser
    std::string sourceName;
    double vstart{};
    double vend{};
    double vstep{};
    std::vector<double> sweep_values;   // All sweep values from vstart to vend
};

struct DCMat{
    // Reuse LHS and RHS
    arma::mat LHS;  // Left-hand side matrix
    arma::vec RHS;  // Right-hand side vector
};

// A structure for multi-sweep results
struct DCResult{
    double sweepValue;        //  source values
    std::string sweepName;    //  source names
    arma::vec solution;
};

struct DCSimulator {
    DCSweepSpec dcsweep;              // sweep device data
    std::vector<DCResult> vec_dc;     // All DC results
    bool non_linear = false;
};

} // namespace dc

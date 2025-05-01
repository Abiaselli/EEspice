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

struct DC;
struct DCSimulator;

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
    std::vector<double> sweepValues;        // All source values
    std::vector<std::string> sweepNames;    // All source names
    arma::vec solution;
};

struct DCSimulator {
    // We can store multiple sweeps for multi-dim
    std::vector<DCSweepSpec> sweeps;       // All sweep devices
    std::vector<DC> vec_dc;                // All DC results
    bool non_linear = false;
};

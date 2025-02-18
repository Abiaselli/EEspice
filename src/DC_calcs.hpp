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

// struct DCSweepSpec {
//     // Using in the parser
//     std::string sourceName;
//     double vstart{};
//     double vend{};
//     double vstep{};
//     std::vector<double> sweep_values;   // All sweep values from vstart to vend
// };

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

namespace dc{
DCSimulator DCsetup(const CircuitParser &parser, const CKTcircuit &ckt){
    DCSimulator dcSim;

    // Check if the DC simulation is non-linear
    for(const auto &element : ckt.CKTelements)
    {
        std::visit([&](auto &&arg)
                   {
                       if constexpr (std::is_same_v<std::decay_t<decltype(arg)>, NMOS> || 
                                     std::is_same_v<std::decay_t<decltype(arg)>, PMOS> || 
                                     std::is_same_v<std::decay_t<decltype(arg)>, Diode>)
                       {
                        dcSim.non_linear = true;
                       } },
                   element.element);

        if(dcSim.non_linear == true) break;
    }

    // Setup DC sweeps
    dcSim.sweeps =parser.dcSweeps;

    // setup DC voltage points
    for(auto &sweep : dcSim.sweeps){
        double vol = sweep.vstart;
        while(vol <= sweep.vend){
            sweep.sweep_values.push_back(vol);
            vol += sweep.vstep;
        }

        if(sweep.vend - sweep.sweep_values.back() > LargeEpsilon){  //1e-6
            sweep.sweep_values.push_back(sweep.vend);
        }
    }

    return dcSim;
}
} // namespace dc
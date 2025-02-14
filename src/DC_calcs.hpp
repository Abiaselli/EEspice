#pragma once
#include <vector>
#include <string>
#include <armadillo>
#include "circuit_parser.hpp"
#include "CKT.hpp"
#include "device.hpp"
#include "global.hpp"
#include "Newton.hpp"

struct DCConfig{
    std::string srcnam;                 // Name of the source
    double vstart{};                    // Start voltage
    double vend{};                      // End voltage
    double vincr{};                     // Voltage increment
    bool non_linear = false;            // true for non-linear solver, false for linear solver
};

struct DC{
    double vol_point;                   // Voltage point
    arma::vec solution;                 // Node voltages
};

struct DCSimulator{
    const DCConfig dc_config;                   // DC configuration(never changes)
    std::vector<DC> vec_dc;                     // DC history
    std::vector<double> voltage_points;         // Voltage points
    std::vector<VoltageSource> vec_voltages;    // Voltage sources in the DC simulation

    DCSimulator(const DCConfig &config) : dc_config(config) {}
};

DCSimulator DCsetup(const CircuitParser &parser, const CKTcircuit &ckt){
    DCConfig config;
    config.srcnam = parser.dc_srcnam;
    config.vstart = parser.double_vstart;
    config.vend   = parser.double_vend;
    config.vincr  = parser.double_vincr;

    // Check if the DC simulation is non-linear
    for(const auto &element : ckt.CKTelements)
    {
        std::visit([&](auto &&arg)
                   {
                       if constexpr (std::is_same_v<std::decay_t<decltype(arg)>, NMOS> || 
                                     std::is_same_v<std::decay_t<decltype(arg)>, PMOS> || 
                                     std::is_same_v<std::decay_t<decltype(arg)>, Diode>)
                       {
                           config.non_linear = true;
                       } },
                   element.element);

        if(config.non_linear == true) break;
    }

    DCSimulator dc_sim(config);

    // setup DC voltage points
    double vol = config.vstart;
    while(vol <= config.vend){
        dc_sim.voltage_points.push_back(vol);
        vol += config.vincr;
    }

    if(config.vend - dc_sim.voltage_points.back() > LargeEpsilon){  //1e-6
        dc_sim.voltage_points.push_back(config.vend);
    }

    // setup DC voltage sources
    for(const auto &element : ckt.CKTelements){
        std::visit([&](auto &&arg)
                   {
                       if constexpr (std::is_same_v<std::decay_t<decltype(arg)>, VoltageSource>){
                           if(arg.id_str == config.srcnam){
                                dc_sim.vec_voltages.push_back(arg);
                           }
                       } },
                   element.element);
    }
    if(dc_sim.vec_voltages.empty()){
        std::cerr << "Error: DC source not found for DC simulation." << std::endl;
        exit(1);
    }

    return dc_sim;
}
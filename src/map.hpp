#pragma once
#include <iostream>
#include <string>
#include <map>
#include <unordered_map>
#include "models.hpp"
#include "bsim4v82/bsim4v82.hpp"

// Stores the names and ids of all nodes and devices.
// Only map_branch_currents is strarting from 0, all other maps start from 1.
struct Circuitmap{
    std::map<std::string, int> map_nodes; // (node name, node id) starting from 1, 0 is ground
    std::map<std::string, int> map_branch_currents; // Branch current (device name, index in the solution vector)

    std::map<std::string, int> map_voltages;
    std::map<std::string, int> map_resistors;
    std::map<std::string, int> map_capacitors;
    std::map<std::string, int> map_currents;
    std::map<std::string, int> map_diodes;
    std::map<std::string, int> map_vccs;
    std::map<std::string, int> map_vcvs;
    std::map<std::string, int> map_mosfets;

};

struct Modelmap{
    // Level 1 mosfet models
    std::unordered_map<std::string, NMOSModel> nmosModels;
    std::unordered_map<std::string, PMOSModel> pmosModels;
    // BSIM mosfet models
    std::unordered_map<std::string, std::shared_ptr<bsim4::BSIM4model>> bsim4Models;
};
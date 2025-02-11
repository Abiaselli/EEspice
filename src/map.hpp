#pragma once
#include <iostream>
#include <string>
#include <map>
#include <unordered_map>
#include "models.hpp"

struct Circuitmap{
    std::map<std::string, int> map_nodes;

    std::map<std::string, int> map_voltages;
    std::map<std::string, int> map_resistors;
    std::map<std::string, int> map_capacitors;
    std::map<std::string, int> map_currents;
    std::map<std::string, int> map_diodes;
    std::map<std::string, int> map_vccs;
    std::map<std::string, int> map_mosfets;

    // Model maps:
    std::unordered_map<std::string, NMOSModel> nmosModels;
    std::unordered_map<std::string, PMOSModel> pmosModels;
};
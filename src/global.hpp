#pragma once
#include <iostream>
#include <string>
#include <vector>
#include <variant>
#include <map>

// for debug mode
#include "sim_variables.hpp"

// Forward declarations
struct VoltageSource;
struct Pulsevoltage;
struct Diode;
struct VCCS;
struct NMOS;
struct PMOS;
struct CurrentSource;
struct Resistor;
struct Capacitor;
struct CircuitElement;

struct VoltageSource
{   
    std::string id_str;
    int id;
    std::string nodePos_str, nodeNeg_str;
    int nodePos, nodeNeg;
    double value;
};

struct Pulsevoltage
{   
    std::string id_str;
    int id{};
    std::string nodePos_str, nodeNeg_str;
    int nodePos{}, nodeNeg{};
    double t1_pulse{};
    double V1{};
    double V2{};
    double td{};
    double tr{};
    double tf{};
    double pw{};
    double per{};

    int RHS_locate{};
};

struct Diode
{   
    std::string id_str;
    int id{};
    std::string nodePos_str, nodeNeg_str;
    int nodePos{}, nodeNeg{};
    double Is{};
    double VT{};
};

struct VCCS
{   
    std::string id_str;
    int id{};
    std::string node_x_str, node_y_str, node_cx_str, node_cy_str;
    int node_x{}, node_y{}, node_cx{}, node_cy{};
    double value{};
};

struct NMOS
{   
    std::string id_str;
    int id{};
    std::string node_vd_str, node_vg_str, node_vs_str, node_vb_str;
    int node_vd{}, node_vg{}, node_vs{}, node_vb{};
    double W{}, L{};
    std::string model;
};

struct PMOS
{   
    std::string id_str;
    int id{};
    std::string node_vd_str, node_vg_str, node_vs_str, node_vb_str;
    int node_vd{}, node_vg{}, node_vs{}, node_vb{};
    double W{}, L{};
    std::string model;
};

struct CurrentSource
{   
    std::string id_str;
    int id;
    std::string nodePos_str, nodeNeg_str;
    int nodePos, nodeNeg;
    double value;
};

struct Resistor
{   
    std::string id_str;
    int id;
    std::string nodePos_str, nodeNeg_str;
    int nodePos, nodeNeg;
    double value;
};

struct Capacitor
{   
    std::string id_str;
    int id{};
    std::string name{}; // It's used in MOSFETs Eg: M1.1, M1.2, M1.3, M1.4
    std::string nodePos_str, nodeNeg_str;
    int nodePos{}, nodeNeg{};
    double value{};

    double current{};
    double voltage{};
    double charge{};
};

//////////////////////////////////////////////////////////////////////////////////

struct CircuitElement
{
    std::variant<VoltageSource, CurrentSource, Resistor, Capacitor, Pulsevoltage, Diode, NMOS, PMOS, VCCS> element;
};

double convertToValue(const std::string &valueStr)
{

    size_t unitPos = valueStr.find_first_not_of("0123456789.-eE"); // Find the position of the first non-numeric character
    double value = 0;

    try
    {
        value = std::stod(valueStr.substr(0, unitPos)); // Convert the numeric part of the string to double
    }
    catch (const std::invalid_argument &ia)
    {
        std::cerr << "Error: Invalid argument: " << ia.what() << std::endl;
        return 0;
    }
    catch (const std::out_of_range &oor)
    {
        std::cerr << "Error: Out of Range error: " << oor.what() << std::endl;
        return 0;
    }

    if (unitPos != std::string::npos)
    {
        char unit = valueStr[unitPos]; // Get the unit character
        switch (unit)
        {

        case 'k':
            return value * 1000.0; // Kilo
        case 'm':
            return value * 1.0e-3; // Milli
        case 'u':
            return value * 1.0e-6; // Micro
        case 'n':
            return value * 1.0e-9; // Nano
        case 'p':
            return value * 1.0e-12; // Pico
        case 'f':
            return value * 1.0e-15; // Femto
        default:
            std::cerr << "Error: Unknown unit: " << unit << std::endl;
            exit(1);
        }
    }

    return value * 1.0; // No unit or unrecognized unit, assume the value is in base units
}

int convertToNode(const std::string &nodeStr, std::map<std::string, int> &map_nodes)
{   
    if(nodeStr == "0" || nodeStr == "GND" || nodeStr == "gnd"){
        return 0;
    }
    auto it = map_nodes.find(nodeStr);
    if (it != map_nodes.end()) {
        return it->second;
    }
    int new_node = map_nodes.size() + 1;
    map_nodes.emplace(nodeStr, new_node);
    return new_node;
}

int convertToDevice(const std::string &deviceStr, std::map<std::string, int> &map_device){
    auto it = map_device.find(deviceStr);
    if (it != map_device.end()) {
       std::cerr << "Error: Device " << deviceStr << " already exists!" << std::endl;
       exit(1);
    }
    int new_device = map_device.size() + 1;
    map_device.emplace(deviceStr, new_device);
    return new_device;
}

//////////////////////////////////////////////////////////////

void setDebugMode(bool mode)
{
    debugMode = mode;
}

bool isDebugMode()
{
    return debugMode;
}

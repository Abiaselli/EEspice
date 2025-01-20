#pragma once
#include <iostream>
#include <string>
#include <vector>
#include <variant>

// for debug mode
#include "sim_variables.hpp"

// Forward declarations
struct XB_Timer;
struct CircuitElement;
struct Transient;
struct CKTcircuit;
struct multi_timestep;
struct Truncation_error;
struct Capacitor;

struct VoltageSource
{
    int id;
    std::string nodePos_str, nodeNeg_str;
    int nodePos, nodeNeg;
    double value;
};

struct Pulsevoltage
{
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
    int id{};
    std::string nodePos_str, nodeNeg_str;
    int nodePos{}, nodeNeg{};
    double Is{};
    double VT{};
};

struct VCCS
{
    int id{};
    std::string node_x_str, node_y_str, node_cx_str, node_cy_str;
    int node_x{}, node_y{}, node_cx{}, node_cy{};
    double value{};
};

struct NMOS
{
    int id{};
    std::string node_vd_str, node_vg_str, node_vs_str, node_vb_str;
    int node_vd{}, node_vg{}, node_vs{}, node_vb{};
    double W{}, L{};
};

struct PMOS
{
    int id{};
    std::string node_vd_str, node_vg_str, node_vs_str, node_vb_str;
    int node_vd{}, node_vg{}, node_vs{}, node_vb{};
    double W{}, L{};
};

struct CurrentSource
{
    int id;
    std::string nodePos_str, nodeNeg_str;
    int nodePos, nodeNeg;
    double value;
};

struct Resistor
{
    int id;
    std::string nodePos_str, nodeNeg_str;
    int nodePos, nodeNeg;
    double value;
};

struct Capacitor
{
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

int convertToNode(const std::string &nodeStr, std::vector<Nodes> &vec_nodes)
{
    int node = 0;
    for(int i = 0; i < vec_nodes.size(); i++){
        if(nodeStr == vec_nodes.at(i).name){
            node = vec_nodes.at(i).id;
            return node;
        }
    }

    node = vec_nodes.size() + 1;
    vec_nodes.push_back(Nodes{nodeStr, node});
    return node;
}

struct Nodes{
    std::string name;
    int id;
};

//////////////////////////////////////////////////////////////

void setDebugMode(bool mode)
{
    debugMode = mode;
}

bool isDebugMode()
{
    return debugMode;
}

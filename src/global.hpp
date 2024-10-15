#pragma once
#include <iostream>
#include <string>
#include <vector>
#include <variant>

// for debug mode
#include "sim_variables.hpp"

// Forward declarations
class XB_Timer;
class CircuitElement;
class Transient;
class CKTcircuit;
class multi_timestep;
class Truncation_error;
class Capacitor;

class VoltageSource
{
public:
    int id;
    int nodePos, nodeNeg;
    double value;
};

class Pulsevoltage
{
public:
    int id{};
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

class Diode
{
public:
    int id{};
    int nodePos{}, nodeNeg{};
    double Is{};
    double VT{};
};

class VCCS
{
public:
    int id{};
    int node_x{}, node_y{}, node_cx{}, node_cy{};
    double value{};
};

class NMOS
{
public:
    int id{};
    int node_vd{}, node_vg{}, node_vs{}, node_vb{};
    double W{}, L{};
};

class PMOS
{
public:
    int id{};
    int node_vd{}, node_vg{}, node_vs{}, node_vb{};
    double W{}, L{};
};

class CurrentSource
{
public:
    int id;
    int nodePos, nodeNeg;
    double value;
};

class Resistor
{
public:
    int id;
    int nodePos, nodeNeg;
    double value;
};

class Capacitor
{
public:
    int id{};
    std::string name{}; // It's used in MOSFETs Eg: M1.1, M1.2, M1.3, M1.4
    int nodePos{}, nodeNeg{};
    double value{};

    double current{};
    double voltage{};
    double charge{};
};

//////////////////////////////////////////////////////////////////////////////////

class CircuitElement
{
public:
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

//////////////////////////////////////////////////////////////

void setDebugMode(bool mode)
{
    debugMode = mode;
}

bool isDebugMode()
{
    return debugMode;
}
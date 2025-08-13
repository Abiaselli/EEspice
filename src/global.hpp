#pragma once
#include <iostream>
#include <string>
#include <vector>
#include <variant>
#include <map>
#include "sim_variables.hpp"
#include "simulation_exceptions.hpp"
#include "bsim4v82/bsim4v82.hpp"

// Forward declarations


struct VoltageSource
{   
    std::string id_str;
    int id{};
    std::string nodePos_str, nodeNeg_str;
    int nodePos{}, nodeNeg{};
    double value{};             // The value of the voltage source, can be DC or transient
    double amplitude{};         // For AC analysis, the amplitude of the voltage source
    double phase{};             // For AC analysis, the phase of the voltage source
    double acReal{};            // AC real component
    double acImag{};            // AC imaginary component

    std::vector<double> batchValues; // For batch simulation, e.g., [start:step:end] or (value1 value2 value3 ...)
};

struct Pulsevoltage
{   
    std::string id_str;
    int id{};
    std::string nodePos_str, nodeNeg_str;
    int nodePos{}, nodeNeg{};
    double V1{};
    double V2{};
    double td{};
    double tr{};
    double tf{};
    double pw{};
    double per{};

    int RHS_locate{};   // RHS index (starts from 0)
};

// sinusoidal voltage source
struct Sinvoltage{
    std::string id_str;
    int id{};
    std::string nodePos_str, nodeNeg_str;
    int nodePos{}, nodeNeg{};

    double vo{};        // Offset
    double va{};        // Amplitude
    double freq{};      // Frequency
    double td{};        // Delay
    double theta{};     // Damping factor
    double phase{};     // Phase

    int RHS_locate{};   // RHS index (starts from 0)
};

struct Diode
{   
    std::string id_str;
    int id{};
    std::string nodePos_str, nodeNeg_str;
    int nodePos{}, nodeNeg{};
    double Is{};
    double VT{};
    std::vector<double> batchIs;
    std::vector<double> batchVT;
};

struct VCCS
{   
    std::string id_str;
    int id{};
    std::string node_x_str, node_y_str, node_cx_str, node_cy_str;
    int node_x{}, node_y{}, node_cx{}, node_cy{};
    double value{};
    std::vector<double> batchValues;
};

struct VCVS{
    std::string id_str;
    int id{};
    std::string node_x_str, node_y_str, node_cx_str, node_cy_str;
    int node_x{}, node_y{}, node_cx{}, node_cy{};
    double value{};
    std::vector<double> batchValues;
};

enum class MosfetModelType {
    LEVEL1,
    BSIM4V82
};

struct NMOS
{   
    std::string id_str;
    int id{};
    std::string node_vd_str, node_vg_str, node_vs_str, node_vb_str;
    int node_vd{}, node_vg{}, node_vs{}, node_vb{};
    double W{}, L{};
    std::string modelName;
    MosfetModelType modelType;
    bsim4::BSIM4V82 bsim4v82Instance;
    // NMOS(const std::string& name) : id_str(name), bsim4v82Instance(name) {}
    std::vector<double> batchW; // For batch simulation
    std::vector<double> batchL;
};

struct PMOS
{   
    std::string id_str;
    int id{};
    std::string node_vd_str, node_vg_str, node_vs_str, node_vb_str;
    int node_vd{}, node_vg{}, node_vs{}, node_vb{};
    double W{}, L{};
    std::string modelName;
    MosfetModelType modelType;
    bsim4::BSIM4V82 bsim4v82Instance;
    // PMOS(const std::string& name) : id_str(name), bsim4v82Instance(name) {}
    std::vector<double> batchW; // For batch simulation
    std::vector<double> batchL;
};

struct CurrentSource
{   
    std::string id_str;
    int id{};
    std::string nodePos_str, nodeNeg_str;
    int nodePos{}, nodeNeg{};
    double value{};
    std::vector<double> batchValues; // For batch simulation, e.g., [start:step:end] or (value1 value2 value3 ...)
};

struct Resistor
{   
    std::string id_str;
    int id{};
    std::string nodePos_str, nodeNeg_str;
    int nodePos{}, nodeNeg{};
    double value{};
    std::vector<double> batchValues;
};

struct Capacitor
{   
    std::string id_str;
    int id{};
    std::string name{}; // It's used in MOSFETs Eg: M1.1, M1.2, M1.3, M1.4
    std::string nodePos_str, nodeNeg_str;
    int nodePos{}, nodeNeg{};
    double value{};
    std::vector<double> batchValues;

    double current{};
    double voltage{};
    double charge{};
};

//////////////////////////////////////////////////////////////////////////////////

// Store circuit elements in separate vectors by type
struct CircuitElements
{   
    // std::variant<VoltageSource, CurrentSource, Resistor, Capacitor, Pulsevoltage, Diode, NMOS, PMOS, VCCS> element;
    std::vector<VoltageSource> voltageSources;
    std::vector<CurrentSource> currentSources;
    std::vector<Resistor> resistors;
    std::vector<Capacitor> capacitors;
    std::vector<Pulsevoltage> pulseVoltages;
    std::vector<Sinvoltage> sinVoltages;
    std::vector<Diode> diodes;
    std::vector<NMOS> nmos;
    std::vector<PMOS> pmos;
    std::vector<VCCS> vccs;
    std::vector<VCVS> vcvs;
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
        throw std::invalid_argument("Invalid numeric part in '" + valueStr + "': " + ia.what());
    }
    catch (const std::out_of_range &oor)
    {
        std::cerr << "Error: Out of Range error: " << oor.what() << std::endl;
        throw std::out_of_range("Out of range error in '" + valueStr + "': " + oor.what());
    }

    if (unitPos != std::string::npos)
    {
        char unit = valueStr[unitPos]; // Get the unit character
        switch (unit)
        {
        case 'T':
            return value * 1.0e12; // Tera
        case 'G':
            return value * 1.0e9; // Giga
        case 'M':
            return value * 1.0e6; // Mega
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
            throw ParsingException("Error: Unknown unit: " + std::string(1, unit), "UNKNOWN_UNIT");
        }
    }

    return value; // No unit or unrecognized unit, assume the value is in base units
}

int convertToNode(const std::string &nodeStr, std::unordered_map<std::string, int> &map_nodes)
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

int convertToDevice(const std::string &deviceStr, std::unordered_map<std::string, int> &map_device){
    auto it = map_device.find(deviceStr);
    if (it != map_device.end()) {
       throw SetupException("Error: Device " + deviceStr + " already exists!", "DUPLICATE_DEVICE");
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

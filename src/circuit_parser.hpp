// Parser class
#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <variant>
#include <sstream>
#include <fstream>
#include <tuple>
#include <cmath>
#include <map>
#include <filesystem>
#include <unordered_set>

#include <algorithm>
#include <deque>
#include <iomanip>
#include <typeinfo>

#include "global.hpp"
#include "map.hpp"
#include "model_parser.hpp"
#include "DC_calcs.hpp"
#include "AC_calcs.hpp"
#include "simulation_exceptions.hpp"

struct CircuitParser
{
    std::string filename;
    CircuitElements elements;
    // Simulation tasks
    bool is_batch = false; // If true, the circuit is a batch simulation
    bool is_transient = false;
    bool is_dc = false;
    bool is_ac = false;
    bool is_op = false;
    bool timestep_control = true;
    bool use_sparse = false; // If true, use sparse matrices for LHS
    // CKT parameters
    // Transient simulation parameters
    double double_t_end; // This double_t_end can be passed to the CKTcircuit class
    double double_init_h;
    // DC simulation parameters
    dc::DCSweepSpec dcSweep_parser;
    // AC simulation parameters
    ac::ACSweepSpec acSweep_parser;
    // Simulation options
    bool acct = false;        // If true, will print the statistics of the simulation
    bool multithreaded = false; // If true, will use multithreading for simulation
    int num_threads = 1;      // Number of threads to use if multithreading is enabled
    bool pulse_phase_mode = false; // If true, 8th PULSE param is PHASE not NP (ngspice xs mode)

    // .param parameter table
    std::unordered_map<std::string, double> params;

    // Parser Timer
    XB_Timer parseTimer;

    CircuitParser(const std::string &filename) : filename(filename) {}

};

int getMaxNode(const CircuitElements &elements)
{
    int maxNode = 0;

    auto TwoMaxNode = [&maxNode](int nodePos, int nodeNeg){
        maxNode = std::max(maxNode, nodePos);
        maxNode = std::max(maxNode, nodeNeg);
    };
    auto FourMaxNode = [&maxNode](int node1, int node2, int node3, int node4){
        maxNode = std::max({maxNode, node1, node2, node3, node4});
    };

    for(const auto &vol : elements.voltageSources)
    {
        TwoMaxNode(vol.nodePos, vol.nodeNeg);
    }
    for(const auto &cul : elements.currentSources)
    {
        TwoMaxNode(cul.nodePos, cul.nodeNeg);
    }
    for(const auto &res : elements.resistors)
    {
        TwoMaxNode(res.nodePos, res.nodeNeg);
    }
    for(const auto &cap : elements.capacitors)
    {
        TwoMaxNode(cap.nodePos, cap.nodeNeg);
    }
    for(const auto &pulse : elements.pulseVoltages)
    {
        TwoMaxNode(pulse.nodePos, pulse.nodeNeg);
    }
    for(const auto &sin : elements.sinVoltages)
    {
        TwoMaxNode(sin.nodePos, sin.nodeNeg);
    }
    for(const auto &diode : elements.diodes)
    {
        TwoMaxNode(diode.nodePos, diode.nodeNeg);
    }
    for(const auto &nmos : elements.nmos)
    {
        FourMaxNode(nmos.node_vd, nmos.node_vg, nmos.node_vs, nmos.node_vb);
    }
    for(const auto &pmos : elements.pmos)
    {
        FourMaxNode(pmos.node_vd, pmos.node_vg, pmos.node_vs, pmos.node_vb);
    }
    for(const auto &vccs : elements.vccs)
    {
        FourMaxNode(vccs.node_x, vccs.node_y, vccs.node_cx, vccs.node_cy);
    }

    return maxNode;
}

int getInternalMosfetNodes(const CircuitElements &elements)
{
    int internal_nodes = 0;

    for (const auto &nmos : elements.nmos)
    {
        internal_nodes += 3; // For LEVEL1 model
    }
    for (const auto &pmos : elements.pmos)
    {
        internal_nodes += 3; // For LEVEL1 model
    }
    // Add more cases for different model types if needed...
    return internal_nodes;
}

// Parse [start:step:end] or (value1 value2 value3 ...) for batch simulation
std::vector<double> batchVector(std::string valueStr, const std::string &line) {
    std::vector<double> vec;

    // Parse bracket [start:step:end] or [start:step:end]u
    if (valueStr.front() == '[')
    {
        size_t closingBracketPos = valueStr.find(']');
        if (closingBracketPos == std::string::npos) {
            std::cerr << "Error: Missing closing bracket ']' in batchVector: " << line << std::endl;
            return vec;
        }

        double unitMultiplier = 1.0;
        // Check for a unit after the closing bracket
        if (closingBracketPos < valueStr.length() - 1) {
            std::string unitPart = valueStr.substr(closingBracketPos + 1);
            // Reuse convertToValue to get the multiplier, e.g., "1u" -> 1.0e-6
            unitMultiplier = convertToValue("1" + unitPart);
        }

        // Extract the content inside the brackets
        std::string innerContent = valueStr.substr(1, closingBracketPos - 1);

        std::vector<std::string> parts;
        {
            std::stringstream ss(innerContent);
            std::string piece;
            while (std::getline(ss, piece, ':')) {
                parts.push_back(piece);
            }
        }

        if (parts.size() == 3) {
            // Convert parts to value, apply the global multiplier
            double bracketStart = convertToValue(parts[0]) * unitMultiplier;
            double bracketStep  = convertToValue(parts[1]) * unitMultiplier;
            double bracketEnd   = convertToValue(parts[2]) * unitMultiplier;

            while(bracketStart <= bracketEnd){
                vec.push_back(bracketStart);
                bracketStart += bracketStep;
            }
            // Add the end point if it wasn't reached due to floating point inaccuracies
            if(vec.empty() || (bracketEnd - vec.back() > LargeEpsilon)){
                vec.push_back(bracketEnd);
            }
        } else {
            std::cerr << "Error parsing bracket contents in batchVector: " << line << std::endl;
        }
    }
    // Parse parenthetical list (value1 value2 value3 ...)
    else if (valueStr.front() == '(' && valueStr.back() == ')')
    {
        valueStr.erase(0, 1);  // remove the first bracket (
        valueStr.pop_back();   // remove the last bracket )
        // Now split on ' '
        std::vector<std::string> parts;
        {
            std::stringstream ss(valueStr);
            std::string piece;
            while (std::getline(ss, piece, ' ')) {
                // Ignore empty strings that can result from multiple spaces
                if (!piece.empty()) {
                   parts.push_back(piece);
                }
            }
        }
        for(auto &part : parts){
            vec.push_back(convertToValue(part));
        }
    }
    else {
       std::cerr << "Error: Invalid format for batch specification: " << line << std::endl;
    }
    return vec;
}

// Helper function to trim whitespace
std::string trim(const std::string& str) {
    size_t first = str.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) return "";
    size_t last = str.find_last_not_of(" \t\r\n");
    return str.substr(first, (last - first + 1));
}

// Helper function to extract filename from .INCLUDE directive
std::string extractIncludeFilename(const std::string& line) {
    std::istringstream iss(line);
    std::string directive;
    iss >> directive; // Skip .INCLUDE
    
    std::string remaining;
    std::getline(iss, remaining);
    remaining = trim(remaining);
    
    // Handle quoted filenames
    if (!remaining.empty() && remaining.front() == '\'') {
        size_t closingQuote = remaining.find('\'', 1);
        if (closingQuote != std::string::npos) {
            return remaining.substr(1, closingQuote - 1);
        }
    }
    
    // Otherwise, take the first token
    std::istringstream tokenStream(remaining);
    std::string filename;
    tokenStream >> filename;
    return filename;
}

// Forward declaration for recursive parsing
void parseNetlistFile(const std::string& filename, CircuitParser& parser, 
                     Circuitmap& cktmap, Modelmap& modmap, 
                     std::unordered_set<std::string>& includeStack,
                     bool isMainFile = false);

// Tokenize a string respecting brace groups: {expr with spaces} stays as one token
std::vector<std::string> tokenizeBraceAware(const std::string& str) {
    std::vector<std::string> tokens;
    size_t i = 0;
    while (i < str.size()) {
        // Skip whitespace
        while (i < str.size() && (str[i] == ' ' || str[i] == '\t'))
            ++i;
        if (i >= str.size()) break;

        std::string tok;
        if (str[i] == '{') {
            // Consume until matching '}'
            while (i < str.size() && str[i] != '}') {
                tok += str[i++];
            }
            if (i < str.size()) tok += str[i++]; // consume '}'
        } else {
            // Normal token — but handle {expr} embedded after = (e.g. w={VDD})
            while (i < str.size() && str[i] != ' ' && str[i] != '\t') {
                if (str[i] == '{') {
                    // Consume brace group
                    while (i < str.size() && str[i] != '}') {
                        tok += str[i++];
                    }
                    if (i < str.size()) tok += str[i++]; // consume '}'
                } else {
                    tok += str[i++];
                }
            }
        }
        if (!tok.empty()) tokens.push_back(tok);
    }
    return tokens;
}

// Parse .param line: ".param VDD=1.8" or ".param a=1 b=2"
void parseParamLine(const std::string& line, std::unordered_map<std::string, double>& params) {
    // Skip the ".param" prefix
    std::string rest = line;
    size_t paramPos = rest.find_first_of(" \t");
    if (paramPos == std::string::npos) return;
    rest = rest.substr(paramPos);

    // Tokenize respecting braces
    auto tokens = tokenizeBraceAware(rest);

    // Process tokens: each should be "name=value" or "name = value"
    size_t i = 0;
    while (i < tokens.size()) {
        std::string tok = tokens[i];
        std::string name, valueStr;

        size_t eqPos = tok.find('=');
        if (eqPos != std::string::npos && eqPos > 0) {
            // "name=value" in one token
            name = tok.substr(0, eqPos);
            valueStr = tok.substr(eqPos + 1);
            if (valueStr.empty() && i + 1 < tokens.size()) {
                // "name=" then value as next token
                valueStr = tokens[++i];
            }
        } else if (i + 2 < tokens.size() && tokens[i + 1] == "=") {
            // "name = value" as three tokens
            name = tok;
            valueStr = tokens[i + 2];
            i += 2;
        } else if (eqPos == 0 && tok.size() > 1) {
            // "=value" — previous token was name (shouldn't happen normally)
            ++i;
            continue;
        } else {
            ++i;
            continue;
        }

        // Normalize name to lowercase
        std::transform(name.begin(), name.end(), name.begin(), ::tolower);

        // Reject reserved names
        if (name == "time" || name == "temper" || name == "hertz") {
            throw ParsingException("Error: Cannot redefine reserved parameter '" + name + "'", "RESERVED_PARAMETER");
        }

        // Strip braces from value if present
        std::string exprStr = valueStr;
        if (!exprStr.empty() && exprStr.front() == '{') {
            size_t closeBrace = exprStr.rfind('}');
            if (closeBrace != std::string::npos) {
                exprStr = exprStr.substr(1, closeBrace - 1);
            }
        }

        // Check for self-reference
        {
            std::string lowerExpr = exprStr;
            std::transform(lowerExpr.begin(), lowerExpr.end(), lowerExpr.begin(), ::tolower);
            // Simple check: see if the name appears as an identifier in the expression
            size_t searchPos = 0;
            while ((searchPos = lowerExpr.find(name, searchPos)) != std::string::npos) {
                bool startOk = (searchPos == 0 || !std::isalnum(static_cast<unsigned char>(lowerExpr[searchPos - 1])));
                bool endOk = (searchPos + name.size() >= lowerExpr.size() || !std::isalnum(static_cast<unsigned char>(lowerExpr[searchPos + name.size()])));
                if (startOk && endOk) {
                    throw ParsingException("Error: Self-referencing parameter '" + name + "' in .param", "SELF_REFERENCE_PARAMETER");
                }
                searchPos += name.size();
            }
        }

        double val = evaluateExpression(exprStr, params);
        params[name] = val;
        ++i;
    }
}

void parseLine(const std::string &line, CircuitParser &parser, Circuitmap &cktmap, Modelmap &modmap)
{

    std::istringstream iss(line);
    std::string type, id_str, valueStr;
    std::string v_nodePos_str, v_nodeNeg_str;
    std::string v_type;

    iss >> type; // Automatically skips leading whitespace before reading type

    if (type.empty())
    {
        return; // Skip empty lines or lines with only whitespaces
    }
    if(type == ".model" || type == ".MODEL")
    {
        parseModel(iss, line, modmap);
    }
    else if(type == ".options" || type == ".OPTIONS")
    {
        std::string option_token;
        while (iss >> option_token)
        {
            // Convert to lowercase for case-insensitive comparison
            std::string option_lower = option_token;
            std::transform(option_lower.begin(), option_lower.end(), option_lower.begin(), ::tolower);

            // Check for key=value format
            size_t eq_pos = option_lower.find('=');
            if (eq_pos != std::string::npos)
            {
                std::string key = option_lower.substr(0, eq_pos);
                std::string value = option_token.substr(eq_pos + 1);

                if (key == "num_threads")
                {
                    try {
                        parser.num_threads = std::stoi(value);
                        parser.multithreaded = true;
                    }
                    catch (const std::exception& e) {
                        std::cerr << "Warning: Invalid value for num_threads option: " << value << std::endl;
                    }
                }
                else
                {
                    std::cerr << "Warning: Unknown .options parameter: " << option_token << std::endl;
                }
            }
            else
            {
                // Standalone flag
                if (option_lower == "acct")
                {
                    parser.acct = true;
                }
                else if (option_lower == "pulsephase")
                {
                    parser.pulse_phase_mode = true;
                }
                else
                {
                    std::cerr << "Warning: Unknown .options parameter: " << option_token << std::endl;
                }
            }
        }
    }
    else if (type[0] == 'V' || type[0] == 'v')
    {
        id_str = type;
        int v_id = convertToDevice(id_str, cktmap.map_voltages);
        iss >> v_nodePos_str >> v_nodeNeg_str;

        // Read the first parameter after the nodes to determine the source type
        std::string first_param;
        if (iss >> first_param)
        {
            if (first_param == "pulse" || first_param == "PULSE")
            { // It's a pulse source.
                Pulsevoltage pv;
                std::string pulseParamsString;
                std::getline(iss, pulseParamsString);
                // Remove the parentheses
                pulseParamsString.erase(std::remove(pulseParamsString.begin(), pulseParamsString.end(), '('), pulseParamsString.end());
                pulseParamsString.erase(std::remove(pulseParamsString.begin(), pulseParamsString.end(), ')'), pulseParamsString.end());

                // Use brace-aware tokenizer to handle {expr} tokens
                auto pulseParams = tokenizeBraceAware(pulseParamsString);

                pv.id_str = id_str;
                pv.id = v_id;
                pv.nodePos_str = v_nodePos_str;
                pv.nodeNeg_str = v_nodeNeg_str;
                pv.nodePos = convertToNode(v_nodePos_str, cktmap.map_nodes);
                pv.nodeNeg = convertToNode(v_nodeNeg_str, cktmap.map_nodes);

                // Parse with defaults of 0 (actual ngspice defaults applied later in V_pulse_value)
                pv.V1 = pulseParams.size() > 0 ? evaluateNumeric(pulseParams[0], parser.params) : 0.0;
                pv.V2 = pulseParams.size() > 1 ? evaluateNumeric(pulseParams[1], parser.params) : 0.0;
                pv.td = pulseParams.size() > 2 ? evaluateNumeric(pulseParams[2], parser.params) : 0.0;
                pv.tr = pulseParams.size() > 3 ? evaluateNumeric(pulseParams[3], parser.params) : 0.0;
                pv.tf = pulseParams.size() > 4 ? evaluateNumeric(pulseParams[4], parser.params) : 0.0;
                pv.pw = pulseParams.size() > 5 ? evaluateNumeric(pulseParams[5], parser.params) : 0.0;
                pv.per = pulseParams.size() > 6 ? evaluateNumeric(pulseParams[6], parser.params) : 0.0;
                pv.param8 = pulseParams.size() > 7 ? evaluateNumeric(pulseParams[7], parser.params) : 0.0;

                parser.elements.pulseVoltages.emplace_back(pv);
            }
            else if (first_param == "sin" || first_param == "SIN"){
                Sinvoltage sv;
                std::string sinParamsString;
                std::getline(iss, sinParamsString);
                // Remove the parentheses
                sinParamsString.erase(std::remove(sinParamsString.begin(), sinParamsString.end(), '('), sinParamsString.end());
                sinParamsString.erase(std::remove(sinParamsString.begin(), sinParamsString.end(), ')'), sinParamsString.end());

                // Use brace-aware tokenizer to handle {expr} tokens
                auto sinParams = tokenizeBraceAware(sinParamsString);

                sv.id_str = id_str;
                sv.id = v_id;
                sv.nodePos_str = v_nodePos_str;
                sv.nodeNeg_str = v_nodeNeg_str;
                sv.nodePos = convertToNode(v_nodePos_str, cktmap.map_nodes);
                sv.nodeNeg = convertToNode(v_nodeNeg_str, cktmap.map_nodes);

                sv.vo = sinParams.size() > 0 ? evaluateNumeric(sinParams[0], parser.params) : 0.0;
                sv.va = sinParams.size() > 1 ? evaluateNumeric(sinParams[1], parser.params) : 0.0;
                sv.freq = sinParams.size() > 2 ? evaluateNumeric(sinParams[2], parser.params) : 0.0;
                sv.td = sinParams.size() > 3 ? evaluateNumeric(sinParams[3], parser.params) : 0.0;
                sv.theta = sinParams.size() > 4 ? evaluateNumeric(sinParams[4], parser.params) : 0.0;
                sv.phase = sinParams.size() > 5 ? evaluateNumeric(sinParams[5], parser.params) : 0.0;

                parser.elements.sinVoltages.emplace_back(sv);
            }
            else
            {
                // It's a DC and/or AC source.
                VoltageSource vs;
                vs.id_str = id_str;
                vs.id = v_id;
                vs.nodePos_str = v_nodePos_str;
                vs.nodeNeg_str = v_nodeNeg_str;
                vs.nodePos = convertToNode(v_nodePos_str, cktmap.map_nodes);
                vs.nodeNeg = convertToNode(v_nodeNeg_str, cktmap.map_nodes);

                // Set defaults
                vs.value = 0.0;
                vs.amplitude = 0.0;
                vs.phase = 0.0;

                std::string dc_val_str, ac_keyword, ac_amp_str, ac_phase_str;

                // Case 1: V... DC <dc_val> AC ... or V... DC <dc_val> 
                if (first_param == "dc" || first_param == "DC")
                {
                    iss >> dc_val_str;
                    vs.value = evaluateNumeric(dc_val_str, parser.params);

                    // Check for an optional AC part
                    if (iss >> ac_keyword && (ac_keyword == "ac" || ac_keyword == "AC"))
                    {
                        iss >> ac_amp_str;
                        vs.amplitude = evaluateNumeric(ac_amp_str, parser.params);
                        if (iss >> ac_phase_str)
                        {
                            vs.phase = evaluateNumeric(ac_phase_str, parser.params);
                        }
                    }
                }
                // Case 2: V... AC <amp> [<phase>]
                else if (first_param == "ac" || first_param == "AC")
                {
                    iss >> ac_amp_str;
                    vs.amplitude = evaluateNumeric(ac_amp_str, parser.params);
                    if (iss >> ac_phase_str)
                    {
                        vs.phase = evaluateNumeric(ac_phase_str, parser.params);
                    }
                }
                else if (first_param.front() == '[' || first_param.front() == '(')
                {
                    // Batch simulation case
                    parser.is_batch = true;
                    vs.batchValues = batchVector(first_param, line);
                }
                // Case 3: V... <dc_val>
                else
                {
                    vs.value = evaluateNumeric(first_param, parser.params);
                }

                // Set the AC real and imaginary components based on amplitude and phase
                double radians = vs.phase * M_PI / 180.0;
                vs.acReal = vs.amplitude * std::cos(radians);
                vs.acImag = vs.amplitude * std::sin(radians);

                parser.elements.voltageSources.emplace_back(vs);
            }
        }
        else
        {
            throw ParsingException("Error: Incomplete voltage source definition for " + id_str, "INCOMPLETE_VOLTAGE_SOURCE");
        }
    }
    else if (type[0] == 'R' || type[0] == 'r')
    {
        Resistor r;
        r.id_str = type;
        r.id = convertToDevice(r.id_str, cktmap.map_resistors);


        iss >> r.nodePos_str >> r.nodeNeg_str >> valueStr;

        r.nodePos = convertToNode(r.nodePos_str, cktmap.map_nodes);
        r.nodeNeg = convertToNode(r.nodeNeg_str, cktmap.map_nodes);
        if (valueStr.front() == '[' || valueStr.front() == '(')
        {
            // Batch simulation case
            parser.is_batch = true;
            r.batchValues = batchVector(valueStr, line);
        }
        else
        {
            // Single value resistor
            r.value = evaluateNumeric(valueStr, parser.params);
        }

        parser.elements.resistors.emplace_back(r);
    }
    else if (type[0] == 'C' || type[0] == 'c')
    {
        Capacitor c;
        c.id_str = type;
        c.id = convertToDevice(c.id_str, cktmap.map_capacitors);

        iss >> c.nodePos_str >> c.nodeNeg_str >> valueStr;

        c.nodePos = convertToNode(c.nodePos_str, cktmap.map_nodes);
        c.nodeNeg = convertToNode(c.nodeNeg_str, cktmap.map_nodes);
        if (valueStr.front() == '[' || valueStr.front() == '(')
        {
            // Batch simulation case
            parser.is_batch = true;
            c.batchValues = batchVector(valueStr, line);
        }
        else
        {
            // Single value capacitor
            c.value = evaluateNumeric(valueStr, parser.params);
        }
        parser.elements.capacitors.emplace_back(c);
    }
    else if (type[0] == 'I' || type[0] == 'i')
    {
        CurrentSource cs;
        cs.id_str = type;
        cs.id = convertToDevice(cs.id_str, cktmap.map_currents);

        iss >> cs.nodePos_str >> cs.nodeNeg_str >> valueStr;

        cs.nodePos = convertToNode(cs.nodePos_str, cktmap.map_nodes);
        cs.nodeNeg = convertToNode(cs.nodeNeg_str, cktmap.map_nodes);
        
        if (valueStr.front() == '[' || valueStr.front() == '(')
        {
            // Batch simulation case
            parser.is_batch = true;
            cs.batchValues = batchVector(valueStr, line);
        }
        else
        {
            // Normal Current source
            cs.value = evaluateNumeric(valueStr, parser.params);
        }
        
        parser.elements.currentSources.emplace_back(cs);
    }
    else if (type[0] == 'D' || type[0] == 'd')
    {

        Diode d;
        d.id_str = type;
        d.id = convertToDevice(d.id_str, cktmap.map_diodes);
        iss >> d.nodePos_str >> d.nodeNeg_str >> valueStr;

        d.nodePos = convertToNode(d.nodePos_str, cktmap.map_nodes);
        d.nodeNeg = convertToNode(d.nodeNeg_str, cktmap.map_nodes);
        if (valueStr.front() == '[' || valueStr.front() == '(')
        {
            // Batch simulation case
            parser.is_batch = true;
            d.batchIs = batchVector(valueStr, line);
        }
        else
        {
            // Single value diode
            d.Is = evaluateNumeric(valueStr, parser.params);
        }

        iss >> valueStr;
        if (valueStr.front() == '[' || valueStr.front() == '(')
        {
            // Batch simulation case
            parser.is_batch = true;
            d.batchVT = batchVector(valueStr, line);
        }
        else
        {
            // Single value diode
            d.VT = evaluateNumeric(valueStr, parser.params);
        }

        parser.elements.diodes.emplace_back(d);
    }
    else if (type[0] == 'G' || type[0] == 'g')
    {
        VCCS g;
        g.id_str = type;
        g.id = convertToDevice(g.id_str, cktmap.map_vccs);

        iss >> g.node_x_str >> g.node_y_str >> g.node_cx_str >> g.node_cy_str >> valueStr;

        g.node_x = convertToNode(g.node_x_str, cktmap.map_nodes);
        g.node_y = convertToNode(g.node_y_str, cktmap.map_nodes);
        g.node_cx = convertToNode(g.node_cx_str, cktmap.map_nodes);
        g.node_cy = convertToNode(g.node_cy_str, cktmap.map_nodes);
        if (valueStr.front() == '[' || valueStr.front() == '(')
        {
            // Batch simulation case
            parser.is_batch = true;
            g.batchValues = batchVector(valueStr, line);
        }
        else
        {
            // Single value VCCS
            g.value = evaluateNumeric(valueStr, parser.params);
        }

        parser.elements.vccs.emplace_back(g);
    }
    else if (type[0] == 'E' || type[0] == 'e'){
        VCVS e;
        e.id_str = type;
        e.id = convertToDevice(e.id_str, cktmap.map_vcvs);

        iss >> e.node_x_str >> e.node_y_str >> e.node_cx_str >> e.node_cy_str >> valueStr;

        e.node_x = convertToNode(e.node_x_str, cktmap.map_nodes);
        e.node_y = convertToNode(e.node_y_str, cktmap.map_nodes);
        e.node_cx = convertToNode(e.node_cx_str, cktmap.map_nodes);
        e.node_cy = convertToNode(e.node_cy_str, cktmap.map_nodes);
        if (valueStr.front() == '[' || valueStr.front() == '(')
        {
            // Batch simulation case
            parser.is_batch = true;
            e.batchValues = batchVector(valueStr, line);
        }
        else
        {
            // Single value VCVS
            e.value = evaluateNumeric(valueStr, parser.params);
        }

        parser.elements.vcvs.emplace_back(e);
    }

    else if (type[0] == 'M' || type[0] == 'm')
    {
        id_str = type;
        int M_id = convertToDevice(id_str, cktmap.map_mosfets);

        std::string M_node_vd_str, M_node_vg_str, M_node_vs_str, M_node_vb_str;
        std::string M_modelName, parameter;

        iss >> M_node_vd_str >> M_node_vg_str >> M_node_vs_str >> M_node_vb_str >> M_modelName;

        const auto &iter_nmos  = modmap.nmosModels.find(M_modelName);
        const auto &iter_pmos  = modmap.pmosModels.find(M_modelName);
        const auto &iter_bsim4 = modmap.bsim4Models.find(M_modelName);

        // Level 1 NMOS
        if (iter_nmos != modmap.nmosModels.end())
        {
            NMOS mn;
            mn.id_str = id_str;
            // mn.id_str = id_str;
            mn.id = M_id;
            // mn.modelType = MosfetModelType::LEVEL1;

            mn.node_vd_str = M_node_vd_str;
            mn.node_vg_str = M_node_vg_str;
            mn.node_vs_str = M_node_vs_str;
            mn.node_vb_str = M_node_vb_str;

            mn.node_vd = convertToNode(M_node_vd_str, cktmap.map_nodes);
            mn.node_vg = convertToNode(M_node_vg_str, cktmap.map_nodes);
            mn.node_vs = convertToNode(M_node_vs_str, cktmap.map_nodes);
            mn.node_vb = convertToNode(M_node_vb_str, cktmap.map_nodes);
            mn.modelName = M_modelName;

            // Read and parse the W and L parameters with their prefixes
            while (iss >> parameter)
            {
                size_t pos = parameter.find('=');
                if (pos != std::string::npos)
                {
                    std::string key = parameter.substr(0, pos);
                    std::string value = parameter.substr(pos + 1);

                    if (key == "W" || key == "w")
                    {
                        valueStr = value;
                        if (valueStr.front() == '[' || valueStr.front() == '(')
                        {
                            // Batch simulation case
                            parser.is_batch = true;
                            mn.batchW = batchVector(valueStr, line);
                        }
                        else
                        {
                            // Single value W
                            mn.W = evaluateNumeric(valueStr, parser.params);
                        }
                    }
                    else if (key == "L" || key == "l")
                    {
                        valueStr = value;
                        if (valueStr.front() == '[' || valueStr.front() == '(')
                        {
                            // Batch simulation case
                            parser.is_batch = true;
                            mn.batchL = batchVector(valueStr, line);
                        }
                        else
                        {
                            // Single value L
                            mn.L = evaluateNumeric(valueStr, parser.params);
                        }
                    }
                }
            }

            parser.elements.nmos.emplace_back(mn);
        }
        // Level 1 PMOS
        else if (iter_pmos != modmap.pmosModels.end())
        {
            PMOS mp;
            mp.id_str = id_str;
            // mp.id_str = id_str;
            mp.id = M_id;
            // mp.modelType = MosfetModelType::LEVEL1;

            mp.node_vd_str = M_node_vd_str;
            mp.node_vg_str = M_node_vg_str;
            mp.node_vs_str = M_node_vs_str;
            mp.node_vb_str = M_node_vb_str;

            mp.node_vd = convertToNode(M_node_vd_str, cktmap.map_nodes);  
            mp.node_vg = convertToNode(M_node_vg_str, cktmap.map_nodes);
            mp.node_vs = convertToNode(M_node_vs_str, cktmap.map_nodes);
            mp.node_vb = convertToNode(M_node_vb_str, cktmap.map_nodes);
            mp.modelName = M_modelName;

            while (iss >> parameter)
            {
                size_t pos = parameter.find('=');
                if (pos != std::string::npos)
                {
                    std::string key = parameter.substr(0, pos);
                    std::string value = parameter.substr(pos + 1);

                    if (key == "W" || key == "w")
                    {
                        valueStr = value;
                        if (valueStr.front() == '[' || valueStr.front() == '(')
                        {
                            // Batch simulation case
                            parser.is_batch = true;
                            mp.batchW = batchVector(valueStr, line);
                        }
                        else
                        {
                            // Single value W
                            mp.W = evaluateNumeric(valueStr, parser.params);
                        }
                    }
                    else if (key == "L" || key == "l")
                    {
                        valueStr = value;
                        if (valueStr.front() == '[' || valueStr.front() == '(')
                        {
                            // Batch simulation case
                            parser.is_batch = true;
                            mp.batchL = batchVector(valueStr, line);
                        }
                        else
                        {
                            // Single value L
                            mp.L = evaluateNumeric(valueStr, parser.params);
                        }
                    }
                }
            }

            parser.elements.pmos.emplace_back(mp);
        }
        // BSIM4v82
        else if (iter_bsim4 != modmap.bsim4Models.end())
        {
            // Parse BSIM4 instance
            BSIM4 b4;
            b4.id_str = id_str;
            b4.id = M_id;
            b4.node_vd_str = M_node_vd_str;
            b4.node_vg_str = M_node_vg_str;
            b4.node_vs_str = M_node_vs_str;
            b4.node_vb_str = M_node_vb_str;

            b4.node_vd = convertToNode(M_node_vd_str, cktmap.map_nodes);
            b4.node_vg = convertToNode(M_node_vg_str, cktmap.map_nodes);
            b4.node_vs = convertToNode(M_node_vs_str, cktmap.map_nodes);
            b4.node_vb = convertToNode(M_node_vb_str, cktmap.map_nodes);
            b4.modelName = M_modelName;
            // Read and parse the W and L parameters with their prefixes
            while (iss >> parameter)
            {
                size_t pos = parameter.find('=');
                if (pos != std::string::npos)
                {
                    std::string key = parameter.substr(0, pos);
                    std::string value = parameter.substr(pos + 1);

                    if (key == "W" || key == "w")
                    {
                        valueStr = value;
                        if (valueStr.front() == '[' || valueStr.front() == '(')
                        {
                            // Batch simulation case
                            parser.is_batch = true;
                            b4.batchW = batchVector(valueStr, line);
                        }
                        else
                        {
                            // Single value W
                            b4.W = evaluateNumeric(valueStr, parser.params);
                        }
                    }
                    else if (key == "L" || key == "l")
                    {
                        valueStr = value;
                        if (valueStr.front() == '[' || valueStr.front() == '(')
                        {
                            // Batch simulation case
                            parser.is_batch = true;
                            b4.batchL = batchVector(valueStr, line);
                        }
                        else
                        {
                            // Single value L
                            b4.L = evaluateNumeric(valueStr, parser.params);
                        }
                    }
                }
            }
            // Move to CKTinstanceSetup!!!
            // b4.bsim4v82Instance = bsim4::paserBSIM4instance(id_str, iter_bsim4->second, b4.node_vd, b4.node_vg, b4.node_vs, b4.node_vb, b4.W, b4.L);
            parser.elements.bsim4.emplace_back(b4);
        }
        else
        {
            throw ParsingException("Error: Unknown MOSFET model: " + M_modelName, "UNKNOWN_MOSFET_MODEL");
        }
    }
    else if (type == ".op" || type == ".OP")
    {
        parser.is_op = true;
    }
    else if (type == ".tran" || type == ".TRAN")
    {
        // Use remaining line with brace-aware tokenizer for {expr} support
        std::string tranRest;
        std::getline(iss, tranRest);
        auto tranParams = tokenizeBraceAware(tranRest);
        if (tranParams.size() < 2) {
            throw ParsingException("Error: .tran requires at least 2 parameters", "INCOMPLETE_TRAN");
        }

        parser.double_init_h = evaluateNumeric(tranParams[0], parser.params);
        parser.double_t_end = evaluateNumeric(tranParams[1], parser.params);
        parser.is_transient = true;
    }
    else if (type == ".dc" || type == ".DC")
    {
        std::string dcRest;
        std::getline(iss, dcRest);
        auto dcParams = tokenizeBraceAware(dcRest);
        if (dcParams.size() < 4) {
            throw ParsingException("Error: .dc requires 4 parameters", "INCOMPLETE_DC");
        }

        dc::DCSweepSpec spec;
        spec.sourceName = dcParams[0];
        spec.vstart = evaluateNumeric(dcParams[1], parser.params);
        spec.vend   = evaluateNumeric(dcParams[2], parser.params);
        spec.vstep  = evaluateNumeric(dcParams[3], parser.params);
        parser.dcSweep_parser = spec;

        parser.is_dc = true;
    }
    else if (type == ".ac" || type == ".AC")
    {
        std::string acRest;
        std::getline(iss, acRest);
        auto acParams = tokenizeBraceAware(acRest);
        if (acParams.size() < 4) {
            throw ParsingException("Error: .ac requires 4 parameters", "INCOMPLETE_AC");
        }

        std::string string_interval = acParams[0];

        ac::ACSweepSpec spec;
        if (string_interval == "dec" || string_interval == "DEC")
        {
            spec.interval = ac::ACSweepSpec::ACSweepType::DEC;
        }
        else if (string_interval == "oct" || string_interval == "OCT")
        {
            spec.interval = ac::ACSweepSpec::ACSweepType::OCT;
        }
        else if (string_interval == "lin" || string_interval == "LIN")
        {
            spec.interval = ac::ACSweepSpec::ACSweepType::LIN;
        }
        else
        {
            throw ParsingException("Error: Unknown AC sweep interval: " + string_interval, "UNKNOWN_AC_SWEEP_INTERVAL");
        }
        spec.numpts = evaluateNumeric(acParams[1], parser.params);
        spec.fstart = evaluateNumeric(acParams[2], parser.params);
        spec.fstop  = evaluateNumeric(acParams[3], parser.params);
        parser.acSweep_parser = spec;
        parser.is_ac = true;
    }
    else if (type == ".sparse" || type == ".SPARSE")
    {
        parser.use_sparse = true;  // Enable sparse matrix mode for LHS
    }

    else
    {

        throw ParsingException("Error: Unknown element type: " + type, "UNKNOWN_ELEMENT_TYPE");
    }
}

// Modified parser_netlist function - now just calls the recursive parser
void parser_netlist(CircuitParser &parser, Circuitmap &cktmap, Modelmap &modmap) {
    parser.parseTimer.start();
    std::unordered_set<std::string> includeStack;
    
    // Get canonical path of the main file
    std::filesystem::path mainPath(parser.filename);
    try {
        mainPath = std::filesystem::canonical(mainPath);
    } catch (const std::filesystem::filesystem_error&) {
        mainPath = std::filesystem::absolute(mainPath);
    }
    
    parseNetlistFile(mainPath.string(), parser, cktmap, modmap, includeStack, true);
    parser.parseTimer.stop();
}

// New recursive file parser that handles .INCLUDE and continuation lines
void parseNetlistFile(const std::string& filename, CircuitParser& parser, 
                     Circuitmap& cktmap, Modelmap& modmap, 
                     std::unordered_set<std::string>& includeStack,
                     bool isMainFile) {
    
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw ParsingException("Error opening file: " + filename, "FILE_OPEN_ERROR");
    }
    
    // Add current file to include stack to detect circular includes
    includeStack.insert(filename);
    
    std::string line;
    bool isFirstLine = true;
    
    while (std::getline(file, line)) {
        // Skip the title line (first line) only for the main file
        if (isFirstLine && isMainFile) {
            isFirstLine = false;
            continue;
        }
        isFirstLine = false;
        
        // Remove whole-line comments (* as first non-whitespace character)
        {
            size_t firstNonSpace = line.find_first_not_of(" \t");
            if (firstNonSpace != std::string::npos && line[firstNonSpace] == '*') {
                line.clear();
            }
        }
        
        // Trim and skip empty lines
        line = trim(line);
        if (line.empty()) {
            continue;
        }
        
        // Check first token to determine line type
        std::istringstream iss(line);
        std::string firstToken;
        iss >> firstToken;
        
        if (!firstToken.empty()) {
            std::string lowerToken = firstToken;
            std::transform(lowerToken.begin(), lowerToken.end(), lowerToken.begin(), ::tolower);
            
            // Handle .END directive
            if (lowerToken == ".end") {
                // .END terminates netlist parsing
                break;
            }
            
            // Handle .INCLUDE directive
            if (lowerToken == ".include") {
                std::string includeFile = extractIncludeFilename(line);
                
                if (includeFile.empty()) {
                    throw ParsingException("Error: .INCLUDE directive missing filename", "MISSING_INCLUDE_FILENAME");
                }
                
                // Resolve relative paths
                std::filesystem::path includePath(includeFile);
                if (!includePath.is_absolute()) {
                    std::filesystem::path currentDir = std::filesystem::path(filename).parent_path();
                    includePath = currentDir / includePath;
                }
                
                // Canonicalize the path
                try {
                    includePath = std::filesystem::canonical(includePath);
                } catch (const std::filesystem::filesystem_error&) {
                    includePath = std::filesystem::absolute(includePath);
                }
                
                std::string canonicalPath = includePath.string();
                
                // Check for circular includes
                if (includeStack.find(canonicalPath) != includeStack.end()) {
                    std::cerr << "Error: Circular include detected for file: " << canonicalPath << std::endl;
                    std::cerr << "Include chain: ";
                    for (const auto& f : includeStack) {
                        std::cerr << f << " -> ";
                    }
                    throw ParsingException("Circular include detected or file not found: " + canonicalPath, "CIRCULAR_INCLUDE_OR_FILE_NOT_FOUND");
                }
                
                // Recursively parse the included file
                parseNetlistFile(canonicalPath, parser, cktmap, modmap, includeStack, false);
                continue; // Skip further processing of this line
            }

            // Handle .OUTPUT directive (removed to favor command-line flags)
            if (lowerToken == ".output") {
                continue;
            }

            // TODO: When .title and .lib are implemented, add them here:
            // else if (lowerToken == ".title" || lowerToken == ".lib") {
            //     // These directives don't support continuation lines
            //     parseLine(line, parser, cktmap, modmap);
            //     continue;
            // }
            
            // For all other lines, check for continuation lines
            std::string completeLine = line;
            
            // Process continuation lines
            std::streampos lastPos = file.tellg();
            std::string nextLine;
            
            while (std::getline(file, nextLine)) {
                // Remove whole-line comments from continuation line
                {
                    size_t contFirstNonSpace = nextLine.find_first_not_of(" \t");
                    if (contFirstNonSpace != std::string::npos && nextLine[contFirstNonSpace] == '*') {
                        nextLine.clear();
                    }
                }
                
                // Check if this is a continuation line
                size_t firstNonSpace = nextLine.find_first_not_of(" \t");
                if (firstNonSpace != std::string::npos && nextLine[firstNonSpace] == '+') {
                    // Extract content after + and first whitespace
                    size_t contentStart = nextLine.find_first_not_of(" \t", firstNonSpace + 1);
                    if (contentStart != std::string::npos) {
                        completeLine += " " + nextLine.substr(contentStart);
                    }
                    lastPos = file.tellg(); // Update position
                } else {
                    // Not a continuation line, rewind
                    file.seekg(lastPos);
                    break;
                }
            }
            
            // Check for .param directive
            {
                std::istringstream tokenCheck(completeLine);
                std::string firstTok;
                tokenCheck >> firstTok;
                std::string lowerFirst = firstTok;
                std::transform(lowerFirst.begin(), lowerFirst.end(), lowerFirst.begin(), ::tolower);
                if (lowerFirst == ".param") {
                    parseParamLine(completeLine, parser.params);
                    continue;
                }
            }

            // Parse the complete line
            parseLine(completeLine, parser, cktmap, modmap);
        }
    }
    
    // Remove current file from include stack
    includeStack.erase(filename);
    file.close();
}
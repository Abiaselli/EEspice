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

#include <algorithm>
#include <deque>
#include <iomanip>
#include <typeinfo>

#include "global.hpp"
#include "map.hpp"
#include "model_parser.hpp"
#include "DC_calcs.hpp"
#include "AC_calcs.hpp"

struct CircuitParser
{
    std::string filename;
    CircuitElements elements;
    // Simulation tasks
    bool is_transient = false;
    bool is_dc = false;
    bool is_ac = false;
    bool timestep_control = true;
    // CKT parameters
    // Transient simulation parameters
    double double_t_end; // This double_t_end can be passed to the CKTcircuit class
    double double_init_h;
    // DC simulation parameters
    DCSweepSpec dcSweep_parser;
    // AC simulation parameters
    AC::ACSweepSpec acSweep_parser;

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
        if (nmos.modelType == MosfetModelType::LEVEL1)
        {
            internal_nodes += 3; // For LEVEL1 model
        }
    }
    for (const auto &pmos : elements.pmos)
    {
        if (pmos.modelType == MosfetModelType::LEVEL1)
        {
            internal_nodes += 3; // For LEVEL1 model
        }
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


void parseLine(const std::string &line, CircuitParser &parser, Circuitmap &cktmap, Modelmap &modmap)
{

    std::istringstream iss(line);
    std::string type, id_str, valueStr;
    std::string v_nodePos_str, v_nodeNeg_str;
    std::string v_type, t1_pulse;

    iss >> type; // Automatically skips leading whitespace before reading type

    if (type.empty())
    {
        return; // Skip empty lines or lines with only whitespaces
    }
    if(type == ".model" || type == ".MODEL")
    {
        parseModel(iss, line, modmap);
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

                // Split the pulseParamsString into individual parameters
                std::istringstream pulseParamsStream(pulseParamsString);
                std::string v1, v2, td, tr, tf, pw, per;

                pv.id_str = id_str;
                pv.id = v_id;
                pv.nodePos_str = v_nodePos_str;
                pv.nodeNeg_str = v_nodeNeg_str;
                pv.nodePos = convertToNode(v_nodePos_str, cktmap.map_nodes);
                pv.nodeNeg = convertToNode(v_nodeNeg_str, cktmap.map_nodes);

                pulseParamsStream >> v1 >> v2 >> td >> tr >> tf >> pw >> per;

                pv.V1 = convertToValue(v1);
                pv.V2 = convertToValue(v2);
                pv.td = convertToValue(td);
                pv.tr = convertToValue(tr);
                pv.tf = convertToValue(tf);
                pv.pw = convertToValue(pw);
                pv.per = convertToValue(per);

                parser.elements.pulseVoltages.emplace_back(pv);
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
                    vs.value = convertToValue(dc_val_str);

                    // Check for an optional AC part
                    if (iss >> ac_keyword && (ac_keyword == "ac" || ac_keyword == "AC"))
                    {
                        iss >> ac_amp_str;
                        vs.amplitude = convertToValue(ac_amp_str);
                        if (iss >> ac_phase_str)
                        {
                            vs.phase = convertToValue(ac_phase_str);
                        }
                    }
                }
                // Case 2: V... AC <amp> [<phase>]
                else if (first_param == "ac" || first_param == "AC")
                {
                    iss >> ac_amp_str;
                    vs.amplitude = convertToValue(ac_amp_str);
                    if (iss >> ac_phase_str)
                    {
                        vs.phase = convertToValue(ac_phase_str);
                    }
                }
                // Case 3: V... <dc_val>
                else
                {
                    vs.value = convertToValue(first_param);
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
            std::cerr << "Error: Incomplete voltage source definition for " << id_str << std::endl;
            exit(1);
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
        r.value = convertToValue(valueStr);

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
        c.value = convertToValue(valueStr);

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
        
        // Normal Current source
        cs.value = convertToValue(valueStr);
        
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
        d.Is = convertToValue(valueStr);

        iss >> valueStr;
        d.VT = convertToValue(valueStr);

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
        g.value = convertToValue(valueStr);

        parser.elements.vccs.emplace_back(g);
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
            mn.modelType = MosfetModelType::LEVEL1;

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
                        mn.W = convertToValue(valueStr);
                    }
                    else if (key == "L" || key == "l")
                    {
                        valueStr = value;
                        mn.L = convertToValue(valueStr);
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
            mp.modelType = MosfetModelType::LEVEL1;

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
                        mp.W = convertToValue(valueStr);
                    }
                    else if (key == "L" || key == "l")
                    {
                        valueStr = value;
                        mp.L = convertToValue(valueStr);
                    }
                }
            }

            parser.elements.pmos.emplace_back(mp);
        }
        // BSIM4v82
        else if (iter_bsim4 != modmap.bsim4Models.end())
        {
            // Parse BSIM4 instance
            // Create a new NMOS instance
            if (iter_bsim4->second->BSIM4type == bsim4::BSIM4_NMOS){
                NMOS mn;
                mn.id_str = id_str;
                mn.id = M_id;
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
                            mn.W = convertToValue(valueStr);
                        }
                        else if (key == "L" || key == "l")
                        {
                            valueStr = value;
                            mn.L = convertToValue(valueStr);
                        }
                    }
                }
                // Setup modelType
                mn.modelType = MosfetModelType::BSIM4V82;
                mn.bsim4v82Instance = bsim4::paserBSIM4instance(id_str, iter_bsim4->second, mn.node_vd, mn.node_vg, mn.node_vs, mn.node_vb, mn.W, mn.L);
                parser.elements.nmos.emplace_back(mn);
            }
            // Create a new PMOS instance
            else if (iter_bsim4->second->BSIM4type == bsim4::BSIM4_PMOS){
                PMOS mp;
                mp.id_str = id_str;
                mp.id = M_id;

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
                            mp.W = convertToValue(valueStr);
                        }
                        else if (key == "L" || key == "l")
                        {
                            valueStr = value;
                            mp.L = convertToValue(valueStr);
                        }
                    }
                }
                // Paser modelType
                mp.modelType = MosfetModelType::BSIM4V82;
                mp.bsim4v82Instance = bsim4::paserBSIM4instance(id_str, iter_bsim4->second, mp.node_vd, mp.node_vg, mp.node_vs, mp.node_vb, mp.W, mp.L);
                parser.elements.pmos.emplace_back(mp);
            }
            else{
                std::cerr << "Error: Unknown BSIM4 model type for: " << M_modelName << std::endl;
                exit(1);
            }
        }
        else
        {
            std::cerr << "Error: Unknown MOSFET model: " << M_modelName << std::endl;
            exit(1);
        }
    }
    else if (type == ".tran" || type == ".TRAN")
    {
        std::string string_h, string_t_end;

        iss >> string_h >> string_t_end;

        parser.double_init_h = convertToValue(string_h);
        parser.double_t_end = convertToValue(string_t_end);
        parser.is_transient = true;
    }
    else if (type == ".dc" || type == ".DC")
    {
        std::string srcnam, vstart, vend, vincr;
        iss >> srcnam >> vstart >> vend >> vincr;

        DCSweepSpec spec;
        spec.sourceName = srcnam;
        spec.vstart = convertToValue(vstart);
        spec.vend   = convertToValue(vend);
        spec.vstep  = convertToValue(vincr);
        parser.dcSweep_parser = spec;

        parser.is_dc = true;
    }
    else if (type == ".ac" || type == ".AC")
    {
        std::string string_interval, string_numpts, string_fstart, string_fstop;
        iss >> string_interval >> string_numpts >> string_fstart >> string_fstop;

        AC::ACSweepSpec spec;
        if (string_interval == "dec" || string_interval == "DEC")
        {
            spec.interval = AC::ACSweepSpec::ACSweepType::DEC;
        }
        else if (string_interval == "oct" || string_interval == "OCT")
        {
            spec.interval = AC::ACSweepSpec::ACSweepType::OCT;
        }
        else if (string_interval == "lin" || string_interval == "LIN")
        {
            spec.interval = AC::ACSweepSpec::ACSweepType::LIN;
        }
        else
        {
            std::cerr << "Error: Unknown AC sweep interval: " << string_interval << std::endl;
            exit(1);
        }
        spec.numpts = convertToValue(string_numpts);
        spec.fstart = convertToValue(string_fstart);
        spec.fstop  = convertToValue(string_fstop);
        parser.acSweep_parser = spec;
        parser.is_ac = true;
    }
    
    else
    {

        std::cerr << "Error: Unknown element type: " << type << std::endl;
        exit(1);
    }
}

void parser_netlist(CircuitParser &parser, Circuitmap &cktmap, Modelmap &modmap)
{

    std::ifstream file(parser.filename);
    if (!file.is_open())
    {
        std::cerr << "Error opening file: " << parser.filename << std::endl;
        exit(1);
        return;
    }

    std::string line;
    std::getline(file, line); // Read and discard the first line

    while (std::getline(file, line))
    {
        size_t commentPos = line.find('*');

        if (commentPos != std::string::npos)
        {
            line = line.substr(0, commentPos); // Remove comment
        }

        if (!line.empty())
        {
            parseLine(line, parser,cktmap, modmap);
        }
    }
}
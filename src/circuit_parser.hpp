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

struct CircuitParser
{
    std::string filename;
    std::vector<CircuitElement> elements;
    // Simulation tasks
    bool is_transient = false;
    bool is_dc = false;
    // CKT parameters
    int num_mosfets{};
    // Transient simulation parameters
    double double_t_end; // This double_t_end can be passed to the CKTcircuit class
    double double_init_h;
    // DC simulation parameters
    std::vector<DCSweepSpec> dcSweeps;

    CircuitParser(const std::string &filename) : filename(filename) {}

    const std::vector<CircuitElement> &getCircuitElements() const
    {
        return elements;
    }

};

int getMaxNode(const std::vector<CircuitElement> &elements)
{
    int maxNode = 0;
    for (const auto &element : elements)
    {
        std::visit([&maxNode](auto &&arg)
                    {
                        if constexpr (std::is_same_v<std::decay_t<decltype(arg)>, NMOS>)
                        {
                            maxNode = std::max({maxNode, arg.node_vd, arg.node_vg, arg.node_vs, arg.node_vb});
                        }
                        else if constexpr (std::is_same_v<std::decay_t<decltype(arg)>, PMOS>)
                        {
                            maxNode = std::max({maxNode, arg.node_vd, arg.node_vg, arg.node_vs, arg.node_vb});
                        }
                        else if constexpr (std::is_same_v<std::decay_t<decltype(arg)>, VCCS>)
                        {
                            maxNode = std::max({maxNode, arg.node_x, arg.node_y, arg.node_cx, arg.node_cy});
                        }
                        else
                        {
                            maxNode = std::max(maxNode, arg.nodePos);
                            maxNode = std::max(maxNode, arg.nodeNeg);
                        } },
                    element.element);
    }
    return maxNode;
}

// Parse [start:step:end] or (value1 value2 value3 ...) for batch simulation
std::vector<double> batchVector(std::string valueStr, const std::string &line){
    std::vector<double> vec;

    // Parse bracket [start:step:end]
    if(valueStr.front() == '[' && valueStr.back() == ']')
    {   
        valueStr.erase(0, 1);  // remove the first bracket [
        valueStr.pop_back();   // remove the last bracket ]
        // Now split on ':'
        std::vector<std::string> parts;
        {
            std::stringstream ss(valueStr);
            std::string piece;
            while (std::getline(ss, piece, ':')) {
                parts.push_back(piece);
            }
        }
        if (parts.size() == 3) {
            double bracketStart = convertToValue(parts[0]);
            double bracketStep  = convertToValue(parts[1]);
            double bracketEnd   = convertToValue(parts[2]);
            while(bracketStart <= bracketEnd){
                vec.push_back(bracketStart);
                bracketStart += bracketStep;
            }
            if(bracketEnd - vec.back() > LargeEpsilon){
                vec.push_back(bracketEnd);
            }
            
        } else {
            std::cerr << "Error parsing bracket in batchVector: " << line << std::endl;
        }

    }
    // Parse bracket (value1 value2 value3 ...)
    else if(valueStr.front() == '(' && valueStr.back() == ')')
    {
        valueStr.erase(0, 1);  // remove the first bracket (
        valueStr.pop_back();   // remove the last bracket )
        // Now split on ' '
        std::vector<std::string> parts;
        {
            std::stringstream ss(valueStr);
            std::string piece;
            while (std::getline(ss, piece, ' ')) {
                parts.push_back(piece);
            }
        }
        for(auto &part : parts){
            vec.push_back(convertToValue(part));
        }
    }
    else{
       std::cerr << "Error parsing bracket in batchVector: " << line << std::endl;
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

        iss >> v_nodePos_str >> v_nodeNeg_str >> v_type;
        if (v_type == "pulse" || v_type == "PULSE")
        { // Pulse voltage settings
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

            parser.elements.push_back(CircuitElement{pv});
        }
        else
        {
            VoltageSource vs;
            vs.id_str = id_str;
            vs.id = v_id;
            vs.nodePos_str = v_nodePos_str;
            vs.nodeNeg_str = v_nodeNeg_str;
            vs.nodePos = convertToNode(v_nodePos_str, cktmap.map_nodes);
            vs.nodeNeg = convertToNode(v_nodeNeg_str, cktmap.map_nodes);

            // Parse bracket [start:step:end]
            if(v_type.front() == '[' && v_type.back() == ']')
            {   
                v_type.erase(0, 1);  // remove the first bracket [
                v_type.pop_back();   // remove the last bracket ]
                // Now split on ':'
                std::vector<std::string> parts;
                {
                    std::stringstream ss(v_type);
                    std::string piece;
                    while (std::getline(ss, piece, ':')) {
                        parts.push_back(piece);
                    }
                }
                if (parts.size() == 3) {
                    vs.hasBracket = true;
                    vs.bracketStart = convertToValue(parts[0]);
                    vs.bracketStep  = convertToValue(parts[1]);
                    vs.bracketEnd   = convertToValue(parts[2]);
                    
                    // Also store a DCSweepSpec in parser's vector
                    DCSweepSpec spec;
                    spec.sourceName = id_str; // e.g. "Vg"
                    spec.vstart = vs.bracketStart;
                    spec.vend   =  vs.bracketEnd;
                    spec.vstep  = vs.bracketStep;
                    parser.dcSweeps.push_back(spec);
                } else {
                    std::cerr << "Error parsing bracket in voltage source: " << line << std::endl;
                }

            }
            // Normal DC voltage source
            else{
                vs.value = convertToValue(v_type);
            }
            parser.elements.push_back(CircuitElement{vs});
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

        parser.elements.push_back(CircuitElement{r});
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

        parser.elements.push_back(CircuitElement{c});
    }
    else if (type[0] == 'I' || type[0] == 'i')
    {
        CurrentSource cs;
        cs.id_str = type;
        cs.id = convertToDevice(cs.id_str, cktmap.map_currents);

        iss >> cs.nodePos_str >> cs.nodeNeg_str >> valueStr;

        cs.nodePos = convertToNode(cs.nodePos_str, cktmap.map_nodes);
        cs.nodeNeg = convertToNode(cs.nodeNeg_str, cktmap.map_nodes);
        // Parse bracket [start:step:end]
        if(valueStr.front() == '[' && valueStr.back() == ']')
        {   
            valueStr.erase(0, 1);  // remove the first bracket [
            valueStr.pop_back();   // remove the last bracket ]
            // Now split on ':'
            std::vector<std::string> parts;
            {
                std::stringstream ss(valueStr);
                std::string piece;
                while (std::getline(ss, piece, ':')) {
                    parts.push_back(piece);
                }
            }
            if (parts.size() == 3) {
                cs.hasBracket = true;
                cs.bracketStart = convertToValue(parts[0]);
                cs.bracketStep  = convertToValue(parts[1]);
                cs.bracketEnd   = convertToValue(parts[2]);
                
                // Also store a DCSweepSpec in parser's vector
                DCSweepSpec spec;
                spec.sourceName = id_str; 
                spec.vstart = cs.bracketStart;
                spec.vend   = cs.bracketEnd;
                spec.vstep  = cs.bracketStep;
                parser.dcSweeps.push_back(spec);
            } else {
                std::cerr << "Error parsing bracket in current source: " << line << std::endl;
            }

        }
        // Normal Current source
        else{
            cs.value = convertToValue(valueStr);
        }
        parser.elements.push_back(CircuitElement{cs});
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

        parser.elements.push_back(CircuitElement{d});
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

        parser.elements.push_back(CircuitElement{g});
    }

    else if (type[0] == 'M' || type[0] == 'm')
    {
        id_str = type;
        int M_id = convertToDevice(id_str, cktmap.map_mosfets);

        std::string M_node_vd_str, M_node_vg_str, M_node_vs_str, M_node_vb_str;
        std::string M_modelName, parameter;

        parser.num_mosfets = parser.num_mosfets + 1;

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

            parser.elements.push_back(CircuitElement{mn});
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

            parser.elements.push_back((CircuitElement{mp}));
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
                mn.bsim4v82Instance = bsim4::paserBSIM4instance(id_str, iter_bsim4->second, mn.node_vd, mn.node_vg, mn.node_vs, mn.node_vb);
                parser.elements.push_back(CircuitElement{mn});
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
                mp.bsim4v82Instance = bsim4::paserBSIM4instance(id_str, iter_bsim4->second, mp.node_vd, mp.node_vg, mp.node_vs, mp.node_vb);
                parser.elements.push_back((CircuitElement{mp}));
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
        parser.dcSweeps.push_back(spec);

        parser.is_dc = true;
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
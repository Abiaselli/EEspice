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

#include <algorithm>
#include <deque>
#include <iomanip>
#include <typeinfo>

#include "global.hpp"

struct CircuitParser
{
    std::string filename;
    std::vector<CircuitElement> elements;
    double double_t_end; // This double_t_end can be passed to the CKTcircuit class
    double double_init_h;
    int num_mosfets{};
    std::map<std::string, int> map_nodes;

    CircuitParser(const std::string &filename) : filename(filename) {}

    void parser()
    {

        std::ifstream file(filename);

        if (!file.is_open())
        {

            std::cerr << "Error opening file: " << filename << std::endl;

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

                parseLine(line);
            }
        }
    }

    const std::vector<CircuitElement> &getCircuitElements() const
    {

        return elements;
    }

    int getMaxNode() const
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

    void parseLine(const std::string &line)
    {

        std::istringstream iss(line);
        std::string type, valueStr;
        // int v_nodePos, v_nodeNeg;
        std::string v_nodePos_str, v_nodeNeg_str;
        std::string v_type, t1_pulse;

        iss >> type; // Automatically skips leading whitespace before reading type
        // std::cout<<"type is "<<type<<std::endl;

        if (type.empty())
        {
            return; // Skip empty lines or lines with only whitespaces
        }

        if (type[0] == 'V' || type[0] == 'v')
        {
            int v_id = std::stoi(type.substr(1));

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

                pv.id = v_id;
                pv.nodePos_str = v_nodePos_str;
                pv.nodeNeg_str = v_nodeNeg_str;

                pulseParamsStream >> v1 >> v2 >> td >> tr >> tf >> pw >> per;

                pv.V1 = convertToValue(v1);
                pv.V2 = convertToValue(v2);
                pv.td = convertToValue(td);
                pv.tr = convertToValue(tr);
                pv.tf = convertToValue(tf);
                pv.pw = convertToValue(pw);
                pv.per = convertToValue(per);

                elements.push_back(CircuitElement{pv});
            }
            else
            {

                VoltageSource vs;

                vs.id = v_id;
                vs.nodePos_str = v_nodePos_str;
                vs.nodeNeg_str = v_nodeNeg_str;

                vs.value = convertToValue(v_type);

                elements.push_back(CircuitElement{vs});
            }
        }
        else if (type[0] == 'R' || type[0] == 'r')
        {

            Resistor r;

            r.id = std::stoi(type.substr(1));

            iss >> r.nodePos_str >> r.nodeNeg_str >> valueStr;

            r.value = convertToValue(valueStr);

            elements.push_back(CircuitElement{r});
        }
        else if (type[0] == 'C' || type[0] == 'c')
        {

            Capacitor c;

            c.id = std::stoi(type.substr(1));

            iss >> c.nodePos_str >> c.nodeNeg_str >> valueStr;

            c.value = convertToValue(valueStr);

            elements.push_back(CircuitElement{c});
        }
        else if (type[0] == 'I' || type[0] == 'i')
        {

            CurrentSource cs;

            cs.id = std::stoi(type.substr(1));

            iss >> cs.nodePos_str >> cs.nodeNeg_str >> valueStr;

            cs.value = convertToValue(valueStr);

            elements.push_back(CircuitElement{cs});
        }
        else if (type[0] == 'D' || type[0] == 'd')
        {

            Diode d;

            d.id = std::stoi(type.substr(1));

            iss >> d.nodePos_str >> d.nodeNeg_str >> valueStr;

            d.Is = convertToValue(valueStr);

            iss >> valueStr;

            d.VT = convertToValue(valueStr);

            elements.push_back(CircuitElement{d});
        }
        else if (type[0] == 'G' || type[0] == 'g')
        {
            VCCS g;
            g.id = std::stoi(type.substr(1));

            iss >> g.node_x_str >> g.node_y_str >> g.node_cx_str >> g.node_cy_str >> valueStr;

            g.value = convertToValue(valueStr);

            elements.push_back(CircuitElement{g});
        }

        else if (type[0] == 'M' || type[0] == 'm')
        {

            int M_id = std::stoi(type.substr(1));

            std::string M_node_vd_str, M_node_vg_str, M_node_vs_str, M_node_vb_str;

            std::string M_model, parameter;

            num_mosfets = num_mosfets + 1;

            iss >> M_node_vd_str >> M_node_vg_str >> M_node_vs_str >> M_node_vb_str >> M_model;

            if (M_model == "NMOS")
            {

                NMOS mn;

                mn.id = M_id;

                mn.node_vd_str = M_node_vd_str;

                mn.node_vg_str = M_node_vg_str;

                mn.node_vs_str = M_node_vs_str;

                mn.node_vb_str = M_node_vb_str;

                // Read and parse the W and L parameters with their prefixes
                while (iss >> parameter)
                {
                    size_t pos = parameter.find('=');
                    if (pos != std::string::npos)
                    {
                        std::string key = parameter.substr(0, pos);
                        std::string value = parameter.substr(pos + 1);

                        if (key == "W")
                        {
                            valueStr = value;
                            mn.W = convertToValue(valueStr);
                        }
                        else if (key == "L")
                        {
                            valueStr = value;
                            mn.L = convertToValue(valueStr);

                            // std::cout<<"L is "<<mn.L<<std::endl;
                            // std::cout<<"type of L is "<<typeid(mn.L).name() << std::endl;
                        }
                    }
                }

                elements.push_back(CircuitElement{mn});
            }
            else if (M_model == "PMOS")
            {

                PMOS mp;

                mp.id = M_id;

                mp.node_vd_str = M_node_vd_str;

                mp.node_vg_str = M_node_vg_str;

                mp.node_vs_str = M_node_vs_str;

                mp.node_vb_str = M_node_vb_str;

                while (iss >> parameter)
                {
                    size_t pos = parameter.find('=');
                    if (pos != std::string::npos)
                    {
                        std::string key = parameter.substr(0, pos);
                        std::string value = parameter.substr(pos + 1);

                        if (key == "W")
                        {
                            valueStr = value;
                            mp.W = convertToValue(valueStr);
                        }
                        else if (key == "L")
                        {
                            valueStr = value;
                            mp.L = convertToValue(valueStr);
                        }
                    }
                }

                elements.push_back(CircuitElement{mp});
            }
            else
            {
                std::cerr << "Error: Unknown MOSFET model: " << M_model << std::endl;
                exit(1);
            }
        }
        else if (type == ".tran")
        {
            std::string string_h, string_t_end;

            iss >> string_h >> string_t_end;

            double_init_h = convertToValue(string_h);
            double_t_end = convertToValue(string_t_end);
        }
        else
        {

            std::cerr << "Error: Unknown element type: " << type << std::endl;
            exit(1);
        }
    }
};

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

class CircuitParser
{
public:
    std::string filename;
    std::vector<CircuitElement> elements;
    double double_t_end; // This double_t_end can be passed to the CKTcircuit class
    double double_init_h;
    int num_mosfets{};

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
        int v_nodePos, v_nodeNeg;
        std::string v_type, t1_pulse, v1, v2, td, tr, tf, pw, per;

        iss >> type; // Automatically skips leading whitespace before reading type
        // std::cout<<"type is "<<type<<std::endl;

        if (type.empty())
        {
            return; // Skip empty lines or lines with only whitespaces
        }

        if (type[0] == 'V')
        {
            int v_id = std::stoi(type.substr(1));

            iss >> v_nodePos >> v_nodeNeg >> v_type;
            if (v_type == "pulse")
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
                pv.nodePos = v_nodePos;
                pv.nodeNeg = v_nodeNeg;

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
                vs.nodePos = v_nodePos;
                vs.nodeNeg = v_nodeNeg;

                vs.value = convertToValue(v_type);

                elements.push_back(CircuitElement{vs});
            }
        }
        else if (type[0] == 'R')
        {

            Resistor r;

            r.id = std::stoi(type.substr(1));

            iss >> r.nodePos >> r.nodeNeg >> valueStr;

            r.value = convertToValue(valueStr);

            elements.push_back(CircuitElement{r});
        }
        else if (type[0] == 'C')
        {

            Capacitor c;

            c.id = std::stoi(type.substr(1));

            iss >> c.nodePos >> c.nodeNeg >> valueStr;

            c.value = convertToValue(valueStr);

            elements.push_back(CircuitElement{c});
        }
        else if (type[0] == 'I')
        {

            CurrentSource cs;

            cs.id = std::stoi(type.substr(1));

            iss >> cs.nodePos >> cs.nodeNeg >> valueStr;

            cs.value = convertToValue(valueStr);

            elements.push_back(CircuitElement{cs});
        }
        else if (type[0] == 'D')
        {

            Diode d;

            d.id = std::stoi(type.substr(1));

            iss >> d.nodePos >> d.nodeNeg >> valueStr;

            d.Is = convertToValue(valueStr);

            iss >> valueStr;

            d.VT = convertToValue(valueStr);

            elements.push_back(CircuitElement{d});
        }
        else if (type[0] == 'G')
        {
            VCCS g;
            g.id = std::stoi(type.substr(1));

            iss >> g.node_x >> g.node_y >> g.node_cx >> g.node_cy >> valueStr;

            g.value = convertToValue(valueStr);

            elements.push_back(CircuitElement{g});
        }

        else if (type[0] == 'M')
        {

            int M_id = std::stoi(type.substr(1));

            int M_node_vd, M_node_vg, M_node_vs, M_node_vb;

            std::string M_model, parameter;

            num_mosfets = num_mosfets + 1;

            iss >> M_node_vd >> M_node_vg >> M_node_vs >> M_node_vb >> M_model;

            if (M_model == "NMOS")
            {

                NMOS mn;

                mn.id = M_id;

                mn.node_vd = M_node_vd;

                mn.node_vg = M_node_vg;

                mn.node_vs = M_node_vs;

                mn.node_vb = M_node_vb;

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

                mp.node_vd = M_node_vd;

                mp.node_vg = M_node_vg;

                mp.node_vs = M_node_vs;

                mp.node_vb = M_node_vb;

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
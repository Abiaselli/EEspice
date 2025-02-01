#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <variant>

#include "device.hpp"
#include "matrix.hpp"
#include "CKT.hpp"
#include "Transient.hpp"

void save_csv(const CKTcircuit &ckt, const std::vector<Transient> &vec_trans, const Circuitmap &map)
{
    std::ofstream file("final_solution.csv");

    // All std::map elements are <std::string, int>

    // 1) Build a helper vector for nodeIndex -> nodeName
    //    We'll index from 1..ckt.external_nodes
    std::vector<std::string> nodeIndexToName(ckt.external_nodes + 1);
    for (const auto &pair : map.map_nodes)
    {
        int nodeIndex = pair.second;   // e.g. 1,2,3 ...
        if (nodeIndex >= 1 && nodeIndex <= ckt.external_nodes)
        {
            nodeIndexToName[nodeIndex] = pair.first;
        }
    }

    std::vector<std::string> volIndexToName(ckt.no_of_V_sources + 1);
    for (const auto &pair : map.map_voltages)
    {
        int volIndex = pair.second;   // e.g. 1,2,3 ...
        if (volIndex >= 1 && volIndex <= ckt.no_of_V_sources)
        {
            volIndexToName[volIndex] = pair.first;
        }
    }

    // 2) Write the header
    file << "Time, Time Step";
    // Output node voltages in ascending node index
    for (int j = 1; j <= ckt.external_nodes; ++j)
    {
        file << ", Voltage " << nodeIndexToName[j];
    }

    for (int j = 1; j <= ckt.no_of_V_sources; ++j)
    {
        file << ", Current " << volIndexToName[j];
    }

    file << std::endl;

    // 3) Write the data rows
    for (size_t i = 0; i < vec_trans.size(); ++i)
    {
        file << std::scientific << std::setprecision(20);                // Set precision to 20 decimal places
        file << vec_trans.at(i).time_trans << ", " << vec_trans.at(i).h; // Time and Timestep

        for (size_t j = 0; j < ckt.external_nodes; ++j)
        {
            file << ", " << vec_trans.at(i).solution(j, 0); // Voltages
        }

        for (size_t z = ckt.no_of_V_sources; z > 0; --z)
        {
            file << ", " << vec_trans.at(i).solution(ckt.cktdematrix->n_rows - z, 0); // Currents
        }
        // file << ", " << vec_trans.at(i).C_current(0,0);
        file << std::endl;
    }

    file.close();
}
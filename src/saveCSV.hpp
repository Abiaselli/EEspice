#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <variant>
#include <string>
#include <filesystem>
#include <algorithm>
#include <sstream>

#include "device.hpp"
#include "matrix.hpp"
#include "CKT.hpp"
#include "Transient.hpp"
#include "DC.hpp"
#include "batch.hpp"

void save_csv(const std::string &filename, const CKTcircuit &ckt, const std::vector<Transient> &vec_trans, const Circuitmap &map)
{
    std::ofstream file(filename);

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

void save_csv_dc(const std::string &filename, const CKTcircuit &ckt, const std::vector<DC> &vec_dc, const Circuitmap &map){
    std::ofstream file(filename);

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

    // Sanity check: must have at least one DC solution:
    if (vec_dc.empty())
    {
        std::cerr << "Warning: vec_dc is empty, no DC solutions to write.\n";
        file.close();
        return;
    }

    // 2) Write the header
    // Single sweep (as per the DC struct implementation)
    file << "Sweep(" << vec_dc.front().sweepName << ")";
    file << ", ";
    
    for (int j = 1; j <= ckt.external_nodes; ++j)
    {
        file << "Voltage " << nodeIndexToName[j] << ", ";
    }
    for (int j = 1; j <= ckt.no_of_V_sources; ++j)
    {
        file << "Current " << volIndexToName[j] << ", ";
    }

    file << std::endl;

    // 3) Write the data rows
    for (const auto &dc : vec_dc) {
        // a) Print the sweep value
        file << dc.sweepValue;
        file << ", ";
        
        // b) Print node voltages
        for (size_t j = 0; j < ckt.external_nodes; ++j)
        {
            file << dc.solution(j, 0) << ", "; // Voltages
        }
        // c) Print current sources
        for (size_t z = ckt.no_of_V_sources; z > 0; --z)
        {
            file << dc.solution(ckt.cktdematrix->n_rows - z, 0)  << ", "; // Currents
        }

        file << std::endl;
    }

    file.close();
}

namespace batch{
void save_csv_batch(const std::vector<BatchRunResult> &batch_results) {
    const std::string output_dir = "batch_results";
    std::filesystem::create_directory(output_dir);

    for (const auto& run_result : batch_results) {
        // Create a unique filename for this run
        std::string filename_suffix;
        for (const auto& [param, value] : run_result.config) {
            // Sanitize param name for filename
            std::string sanitized_param = param;
            std::replace(sanitized_param.begin(), sanitized_param.end(), '.', '_');
            
            std::ostringstream val_stream;
            val_stream << std::scientific << std::setprecision(4) << value;

            filename_suffix += "_" + sanitized_param + "_" + val_stream.str();
        }
        auto &ckt = run_result.ckt; // Get the circuit state from the run result
        if (run_result.simulation_type == "dc") {
            std::string filename = output_dir + "/dc_solution" + filename_suffix + ".csv";
            const auto& vec_dc_result = std::get<std::vector<DC>>(run_result.results);
            save_csv_dc(filename, ckt, vec_dc_result, ckt.map);
        } else if (run_result.simulation_type == "tran") {
            std::string filename = output_dir + "/tran_solution" + filename_suffix + ".csv";
            const auto& vec_trans_result = std::get<std::vector<Transient>>(run_result.results);
            save_csv(filename, ckt, vec_trans_result, ckt.map);
        } // TODO: Add AC results saving
    }
}

}// namespace batch
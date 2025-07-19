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

// We'll index from 1..ckt.external_nodes
std::vector<std::string> buildNodeIndexToNameMap(const CKTcircuit& ckt, const Circuitmap& map) {
    std::vector<std::string> nodeIndexToName(ckt.external_nodes + 1); 
    for (const auto& [name, index] : map.map_nodes) {
        if (index >= 1 && index <= ckt.external_nodes) {
            nodeIndexToName[index] = name;
        }
    }
    return nodeIndexToName;
}

// Overloaded function to pass std::ofstream directly
void save_csv(std::ofstream &file, const CKTcircuit &ckt, const std::vector<Transient> &vec_trans, const Circuitmap &map)
{
    // All std::map elements are <std::string, int>

    // 1) Build a helper vector for nodeIndex -> nodeName
    //    This vector is created to ensure the nodes are arranged in order from 1 to n.
    std::vector<std::string> nodeIndexToName = buildNodeIndexToNameMap(ckt, map);

    // 2) Write the header
    file << "Time,Time Step";
    for (int j = 1; j <= ckt.external_nodes; ++j) {
        file << ",Voltage " << nodeIndexToName[j];
    }
    for (const auto& [name, index] : map.map_branch_currents) {
        file << ",I(" << name << ")";
    }
    file << std::endl;

    // 3) Write the data rows
    for (const auto& trans_result : vec_trans) {
        file << std::scientific << std::setprecision(20);
        file << trans_result.time_trans << "," << trans_result.h;
        for (size_t j = 0; j < ckt.external_nodes; ++j) {
            file << "," << trans_result.solution(j);
        }
        for (const auto& [name, index] : map.map_branch_currents) {
            file << "," << trans_result.solution(index);
        }
        file << std::endl;
    }
}
void save_csv(const std::string &filename, const CKTcircuit &ckt, const std::vector<Transient> &vec_trans, const Circuitmap &map)
{
    std::ofstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Error: Could not open file " << filename << " for writing." << std::endl;
        return;
    }

    // Call the overloaded function that takes std::ofstream directly
    save_csv(file, ckt, vec_trans, map);

    file.close();
}

// Overloaded function to pass std::ofstream directly
void save_csv_dc(std::ofstream &file, const CKTcircuit &ckt, const std::vector<DC> &vec_dc, const Circuitmap &map){
    // Sanity check: must have at least one DC solution:
    if (vec_dc.empty())
    {
        std::cerr << "Warning: vec_dc is empty, no DC solutions to write.\n";
        file.close();
        return;
    }

    // 1) Build a helper vector for nodeIndex -> nodeName
    //    This vector is created to ensure the nodes are arranged in order from 1 to n.
    std::vector<std::string> nodeIndexToName = buildNodeIndexToNameMap(ckt, map);

    // 2) Write the header
    // Single sweep (as per the DC struct implementation)
    file << "Sweep(" << vec_dc.front().sweepName << ")";
    for (int j = 1; j <= ckt.external_nodes; ++j) {
        file << ",Voltage " << nodeIndexToName[j];
    }
    for (const auto& [name, index] : map.map_branch_currents) {
        file << ",I(" << name << ")";
    }
    file << std::endl;

    // 3) Write the data rows
   for (const auto& dc_result : vec_dc) {
        file << dc_result.sweepValue;
        for (size_t j = 0; j < ckt.external_nodes; ++j) {
            file << "," << dc_result.solution(j);
        }
        for (const auto& [name, index] : map.map_branch_currents) {
            file << "," << dc_result.solution(index);
        }
        file << std::endl;
    }
}
void save_csv_dc(const std::string &filename, const CKTcircuit &ckt, const std::vector<DC> &vec_dc, const Circuitmap &map){
    std::ofstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Error: Could not open file " << filename << " for writing." << std::endl;
        return;
    }

    // Call the overloaded function that takes std::ofstream directly
    save_csv_dc(file, ckt, vec_dc, map);

    file.close();
}

namespace batch{
void save_csv_batch(const std::vector<BatchRunResult> &batch_results) {
    const std::string output_dir = "batch_results";
    std::filesystem::create_directory(output_dir);

    int dc_counter = 1;
    int tran_counter = 1;
    int ac_counter = 1;

    for (const auto& run_result : batch_results) {
        auto &ckt = run_result.ckt; // Get the circuit state from the run result
        
        if (run_result.simulation_type == "dc") {
            std::string filename = output_dir + "/dc" + std::to_string(dc_counter) + ".csv";
            std::ofstream file(filename);
            
            // Write circuit information header
            file << "# Circuit Configuration for DC Analysis " << dc_counter << std::endl;
            file << "# Parameters:" << std::endl;
            for (const auto& [param, value] : run_result.config) {
                file << "# " << param << " = " << std::scientific << std::setprecision(6) << value << std::endl;
            }
            file << "# End of Configuration" << std::endl;
            file << std::endl;
            
            // Call the overloaded save_csv_dc function that takes std::ofstream
            const auto& vec_dc_result = std::get<std::vector<DC>>(run_result.results);
            save_csv_dc(file, ckt, vec_dc_result, ckt.map);
            file.close();
            dc_counter++;
            
        } else if (run_result.simulation_type == "tran") {
            std::string filename = output_dir + "/tran" + std::to_string(tran_counter) + ".csv";
            std::ofstream file(filename);
            
            // Write circuit information header
            file << "# Circuit Configuration for Transient Analysis " << tran_counter << std::endl;
            file << "# Parameters:" << std::endl;
            for (const auto& [param, value] : run_result.config) {
                file << "# " << param << " = " << std::scientific << std::setprecision(6) << value << std::endl;
            }
            file << "# End of Configuration" << std::endl;
            file << std::endl;
            
            // Call the overloaded save_csv function that takes std::ofstream
            const auto& vec_trans_result = std::get<std::vector<Transient>>(run_result.results);
            save_csv(file, ckt, vec_trans_result, ckt.map);
            file.close();
            tran_counter++;
            
        } else if (run_result.simulation_type == "ac") {
            std::string filename = output_dir + "/ac" + std::to_string(ac_counter) + ".csv";
            std::ofstream file(filename);
            
            // Write circuit information header
            file << "# Circuit Configuration for AC Analysis " << ac_counter << std::endl;
            file << "# Parameters:" << std::endl;
            for (const auto& [param, value] : run_result.config) {
                file << "# " << param << " = " << std::scientific << std::setprecision(6) << value << std::endl;
            }
            file << "# End of Configuration" << std::endl;
            file << std::endl;
            
            // TODO: Add AC results saving function call here
            // const auto& vec_ac_result = std::get<std::vector<AC::AC>>(run_result.results);
            // save_csv_ac(file, ckt, vec_ac_result, ckt.map);
            file.close();
            ac_counter++;
        }
    }
}

}// namespace batch
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
#include "OP_calcs.hpp"

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
void save_csv_dc(std::ofstream &file, const CKTcircuit &ckt, const std::vector<dc::DCResult> &vec_dc, const Circuitmap &map){
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
void save_csv_dc(const std::string &filename, const CKTcircuit &ckt, const std::vector<dc::DCResult> &vec_dc, const Circuitmap &map){
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

// Overloaded function to pass std::ofstream directly
void save_txt_op(std::ofstream &file, const OPResult &op_result, const Circuitmap &map){

    const auto mosfet_data = op_result.mosfet_data;
    const auto solution = op_result.solution;

    // Write to file with the specified formatting
    file << "Semiconductor Device Operating Points:\n";
    file << "                      --- BSIM4 MOSFETS ---\n";
    if (mosfet_data.empty()) {
        file << "No BSIM4 devices found in the circuit.\n";
        return;
    }

    const size_t devices_per_chunk = 5; // Max devices to print per block

    for (size_t i = 0; i < mosfet_data.size(); i += devices_per_chunk) {
        size_t end = std::min(i + devices_per_chunk, mosfet_data.size());

        // Helper lambda to print a row of data for the current chunk
        auto print_row = [&](const std::string& label, auto value_extractor) {
            file << std::left << std::setw(8) << label;
            for (size_t j = i; j < end; ++j) {
                file << std::right << std::setw(11) << value_extractor(mosfet_data[j]);
            }
            file << "\n";
        };

        // Print string properties
        print_row("Name:",  [](const auto& d) { return d.name; });
        print_row("Model:", [](const auto& d) { return d.model; });

        // Set formatting for all numeric values
        file << std::scientific << std::setprecision(2);
        
        // Print numeric properties
        print_row("Id:",    [](const auto& d) { return d.id; });
        print_row("Vgs:",   [](const auto& d) { return d.vgs; });
        print_row("Vds:",   [](const auto& d) { return d.vds; });
        print_row("Vbs:",   [](const auto& d) { return d.vbs; });
        print_row("Vth:",   [](const auto& d) { return d.vth; });
        print_row("Vdsat:", [](const auto& d) { return d.vdsat; });
        print_row("Gm:",    [](const auto& d) { return d.gm; });
        print_row("Gds:",   [](const auto& d) { return d.gds; });
        print_row("Gmb",    [](const auto& d) { return d.gmb; }); // Note: Gmb label from example
        print_row("Cbd:",   [](const auto& d) { return d.cbd; });
        print_row("Cbs:",   [](const auto& d) { return d.cbs; });
        
        // Reset formatting and add a newline if there are more chunks to print
        file << std::defaultfloat; 
        if (end < mosfet_data.size()) {
            file << "\n";
        }
    }

    file << std::endl;

    // Print the operating point solution to file
    // Temporarily redirect std::cout to the file stream
    auto cout_buffer = std::cout.rdbuf();
    std::cout.rdbuf(file.rdbuf());
    
    printOperatingPointWithNames(solution, map);
    
    // Restore std::cout to its original buffer
    std::cout.rdbuf(cout_buffer);
    file << "\n";
}

void save_txt_op(const std::string &filename, const OPResult &op_result, const Circuitmap &map){
    std::ofstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Error: Could not open file " << filename << " for writing." << std::endl;
        return;
    }
    save_txt_op(file, op_result, map);
    file.close();
}

void save_csv_ac(std::ofstream &file, const CKTcircuit &ckt, const std::vector<ac::ACResult> &vec_ac, const Circuitmap &map, const ac::ACResultType type){
    // Sanity check: must have at least one AC solution:
    if (vec_ac.empty())
    {
        std::cerr << "Warning: vec_ac is empty, no AC solutions to write.\n";
        file.close();
        return;
    }

    // 1) Build a helper vector for nodeIndex -> nodeName
    //    This vector is created to ensure the nodes are arranged in order from 1 to n.
    std::vector<std::string> nodeIndexToName = buildNodeIndexToNameMap(ckt, map);

    // 2) Write the header
    file << "frequency";
    
    // Add voltage node headers
    for (int j = 1; j <= ckt.external_nodes; ++j) {
        if (type == ac::ACResultType::RealImag) {
            file << ",v(" << nodeIndexToName[j] << ") (real),v(" << nodeIndexToName[j] << ") (img)";
        } else { // MagPhase
            file << ",v(" << nodeIndexToName[j] << ") (mag),v(" << nodeIndexToName[j] << ") (phase)";
        }
    }
    
    // Add branch current headers
    for (const auto& [name, index] : map.map_branch_currents) {
        if (type == ac::ACResultType::RealImag) {
            file << ",i(" << name << ") (real),i(" << name << ") (img)";
        } else { // MagPhase
            file << ",i(" << name << ") (mag),i(" << name << ") (phase)";
        }
    }
    file << std::endl;

    // 3) Write the data rows
    if (type == ac::ACResultType::RealImag) {
        // Use complex solution directly
        for (const auto& ac_result : vec_ac) {
            file << std::scientific << std::setprecision(20);
            file << ac_result.freq;
            
            // Write node voltages (real and imaginary parts)
            for (size_t j = 0; j < ckt.external_nodes; ++j) {
                file << "," << std::real(ac_result.solution(j)) << "," << std::imag(ac_result.solution(j));
            }
            
            // Write branch currents (real and imaginary parts)
            for (const auto& [name, index] : map.map_branch_currents) {
                file << "," << std::real(ac_result.solution(index)) << "," << std::imag(ac_result.solution(index));
            }
            file << std::endl;
        }
    } else { // MagPhase
        // Convert to polar form first
        std::vector<ac::ACResultPolar> polarResults = ac::convertToPolar(vec_ac);
        
        for (const auto& polar_result : polarResults) {
            file << std::scientific << std::setprecision(20);
            file << polar_result.freq;
            
            // Write node voltages (magnitude and phase)
            for (size_t j = 0; j < ckt.external_nodes; ++j) {
                file << "," << polar_result.magnitudes(j) << "," << polar_result.phases_rad(j);
            }
            
            // Write branch currents (magnitude and phase)
            for (const auto& [name, index] : map.map_branch_currents) {
                file << "," << polar_result.magnitudes(index) << "," << polar_result.phases_rad(index);
            }
            file << std::endl;
        }
    }
}
void save_csv_ac(const std::string &filename, const CKTcircuit &ckt, const std::vector<ac::ACResult> &vec_ac, const Circuitmap &map, const ac::ACResultType type){
    std::ofstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Error: Could not open file " << filename << " for writing." << std::endl;
        return;
    }
    save_csv_ac(file, ckt, vec_ac, map, type);
    file.close();
}

namespace batch{
// Helper function to write error information to CSV file
void write_error_csv(std::ofstream &file, const BatchRunResult &run_result, const std::string &analysis_type, int counter) {
    file << "# Circuit Configuration for " << analysis_type << " Analysis " << counter << std::endl;
    file << "# Parameters:" << std::endl;
    for (const auto& [param, value] : run_result.config) {
        file << "# " << param << " = " << std::scientific << std::setprecision(6) << value << std::endl;
    }
    file << "# SIMULATION FAILED" << std::endl;
    file << "# Error Type: " << run_result.error_type << std::endl;
    file << "# Error Message: " << run_result.error_message << std::endl;
    file << "# End of Configuration" << std::endl;
    file << std::endl;
    file << "# No simulation data available due to failure" << std::endl;
}

void save_csv_batch(const std::vector<BatchRunResult> &batch_results) {
    const std::string output_dir = "batch_results";
    std::filesystem::create_directory(output_dir);

    int dc_counter = 1;
    int tran_counter = 1;
    int ac_counter = 1;
    int op_counter = 1;
    
    int successful_count = 0;
    int failed_count = 0;

    for (const auto& run_result : batch_results) {
        if (run_result.simulation_type == "op"){
            std::string filename = output_dir + "/op" + std::to_string(op_counter) + ".txt";
            std::ofstream file(filename);
            
            if (!run_result.success) {
                // Write error information for OP analysis
                file << "# Circuit Configuration for OP Analysis " << op_counter << std::endl;
                file << "# Parameters:" << std::endl;
                for (const auto& [param, value] : run_result.config) {
                    file << "# " << param << " = " << std::scientific << std::setprecision(6) << value << std::endl;
                }
                file << "# SIMULATION FAILED" << std::endl;
                file << "# Error Type: " << run_result.error_type << std::endl;
                file << "# Error Message: " << run_result.error_message << std::endl;
                file << "# End of Configuration" << std::endl;
                file << std::endl;
                file << "# No simulation data available due to failure" << std::endl;
                failed_count++;
            } else {
                // Write circuit information header
                file << "# Circuit Configuration for OP Analysis " << op_counter << std::endl;
                file << "# Parameters:" << std::endl;
                for (const auto& [param, value] : run_result.config) {
                    file << "# " << param << " = " << std::scientific << std::setprecision(6) << value << std::endl;
                }
                file << "# End of Configuration" << std::endl;
                file << std::endl;
                
                // Extract OPResult and call save_txt_op with mosfet_data
                const auto& op_result = std::get<OPResult>(run_result.results);
                save_txt_op(file, op_result, run_result.ckt.map);
                successful_count++;
            }
            file.close();
            op_counter++;
            
        } else if (run_result.simulation_type == "dc") {
            std::string filename = output_dir + "/dc" + std::to_string(dc_counter) + ".csv";
            std::ofstream file(filename);
            
            if (!run_result.success) {
                write_error_csv(file, run_result, "DC", dc_counter);
                failed_count++;
            } else {
                // Write circuit information header
                file << "# Circuit Configuration for DC Analysis " << dc_counter << std::endl;
                file << "# Parameters:" << std::endl;
                for (const auto& [param, value] : run_result.config) {
                    file << "# " << param << " = " << std::scientific << std::setprecision(6) << value << std::endl;
                }
                file << "# End of Configuration" << std::endl;
                file << std::endl;
                
                // Call the overloaded save_csv_dc function that takes std::ofstream
                const auto& vec_dc_result = std::get<std::vector<dc::DCResult>>(run_result.results);
                save_csv_dc(file, run_result.ckt, vec_dc_result, run_result.ckt.map);
                successful_count++;
            }
            file.close();
            dc_counter++;
            
        } else if (run_result.simulation_type == "tran") {
            std::string filename = output_dir + "/tran" + std::to_string(tran_counter) + ".csv";
            std::ofstream file(filename);
            
            if (!run_result.success) {
                write_error_csv(file, run_result, "Transient", tran_counter);
                failed_count++;
            } else {
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
                save_csv(file, run_result.ckt, vec_trans_result, run_result.ckt.map);
                successful_count++;
            }
            file.close();
            tran_counter++;
            
        } else if (run_result.simulation_type == "ac") {
            std::string filename = output_dir + "/ac" + std::to_string(ac_counter) + ".csv";
            std::ofstream file(filename);
            
            if (!run_result.success) {
                write_error_csv(file, run_result, "AC", ac_counter);
                failed_count++;
            } else {
                // Write circuit information header
                file << "# Circuit Configuration for AC Analysis " << ac_counter << std::endl;
                file << "# Parameters:" << std::endl;
                for (const auto& [param, value] : run_result.config) {
                    file << "# " << param << " = " << std::scientific << std::setprecision(6) << value << std::endl;
                }
                file << "# End of Configuration" << std::endl;
                file << std::endl;
                
                // Call the overloaded save_csv_ac function that takes std::ofstream
                const auto& vec_ac_result = std::get<std::vector<ac::ACResult>>(run_result.results);
                save_csv_ac(file, run_result.ckt, vec_ac_result, run_result.ckt.map, ac::ACResultType::MagPhase);
                successful_count++;
            }
            file.close();
            ac_counter++;
        }
    }
    
    // Print batch simulation summary
    std::cout << std::endl << "Batch Simulation Summary:" << std::endl;
    std::cout << "Total simulations: " << batch_results.size() << std::endl;
    std::cout << "Successful: " << successful_count << std::endl;
    std::cout << "Failed: " << failed_count << std::endl;
}

}// namespace batch
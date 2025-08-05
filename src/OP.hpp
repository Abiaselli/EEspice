#pragma once
#include <iostream>
#include <iomanip>
#include <armadillo>
#include "Newton.hpp"
#include "CKT.hpp"
#include "map.hpp"
#include "SPICEcompatible.hpp"
#include "OP_calcs.hpp"
#include "bsim4v82/bsim4v82ask.hpp"


// Function to perform operating point analysis
arma::vec OperatingPointAnalysis(CKTcircuit &ckt, const Modelmap &modmap, bool non_linear){
    // Set the flags for SPICE compatibility
    ckt.spiceCompatible.setFlagsOP();

    // Get the initial LHS and RHS matrices
    arma::mat init_LHS = ckt.cktdematrix->get_init_LHS();
    arma::vec init_RHS = ckt.cktdematrix->get_init_RHS();

    // Solve the linear system to get the operating point
    if(non_linear){
        // If the circuit is non-linear, use the Newton-Raphson method
        return NewtonRaphson_system(ckt, init_LHS, init_RHS, modmap);
    } else {
        // If the circuit is linear, use the direct solver
        return arma::solve(init_LHS, init_RHS);
    }

    // Todo: Add gmin stepping and stepping sources...
}

std::vector<MosfetOpData> extractMosfetData(const CKTcircuit &ckt) {
    std::vector<MosfetOpData> mosfet_data;
    mosfet_data.reserve(ckt.CKTelements.nmos.size() + ckt.CKTelements.pmos.size());

    // A helper lambda to populate the data for any BSIM4 device
    auto extract_data = [](const auto& mosfet) -> MosfetOpData {
        MosfetOpData data;
        const auto& inst = mosfet.bsim4v82Instance;

        // Use a default value of 0.0 if the optional is empty or conversion fails.
        constexpr double default_val = 0.0;

        data.name    = mosfet.id_str;
        data.model   = mosfet.modelName;
        data.id      = bsim4::BSIM4ask(inst, bsim4::BSIM4_CD).template get_as<double>().value_or(default_val);
        data.vgs     = bsim4::BSIM4ask(inst, bsim4::BSIM4_VGS).template get_as<double>().value_or(default_val);
        data.vds     = bsim4::BSIM4ask(inst, bsim4::BSIM4_VDS).template get_as<double>().value_or(default_val);
        data.vbs     = bsim4::BSIM4ask(inst, bsim4::BSIM4_VBS).template get_as<double>().value_or(default_val);
        data.vth     = bsim4::BSIM4ask(inst, bsim4::BSIM4_VON).template get_as<double>().value_or(default_val); // VON is Vth
        data.vdsat   = bsim4::BSIM4ask(inst, bsim4::BSIM4_VDSAT).template get_as<double>().value_or(default_val);
        data.gm      = bsim4::BSIM4ask(inst, bsim4::BSIM4_GM).template get_as<double>().value_or(default_val);
        data.gds     = bsim4::BSIM4ask(inst, bsim4::BSIM4_GDS).template get_as<double>().value_or(default_val);
        data.gmb     = bsim4::BSIM4ask(inst, bsim4::BSIM4_GMBS).template get_as<double>().value_or(default_val); // Gmb is GMBS
        data.cbd     = bsim4::BSIM4ask(inst, bsim4::BSIM4_CBD).template get_as<double>().value_or(default_val);
        data.cbs     = bsim4::BSIM4ask(inst, bsim4::BSIM4_CBS).template get_as<double>().value_or(default_val);
        return data;
    };

    for (const auto& nmos : ckt.CKTelements.nmos){
        if (nmos.modelType == MosfetModelType::BSIM4V82) {
            mosfet_data.push_back(extract_data(nmos));
        }
    }

    for (const auto& pmos : ckt.CKTelements.pmos){
        if (pmos.modelType == MosfetModelType::BSIM4V82) {
            mosfet_data.push_back(extract_data(pmos));
        }
    }

    return mosfet_data;
}

void printOperatingPoint(const arma::vec &op_solution, const CKTcircuit &ckt){
    // Print the operating point solution
    arma::vec node_volt_print = op_solution.submat(0, 0, ckt.external_nodes - 1, 0);
    arma::vec current_print = op_solution.submat(op_solution.n_rows - ckt.no_of_V_sources, 0, op_solution.n_rows - 1, 0);
    arma::vec op_solution_print = arma::join_vert(node_volt_print, current_print);
    op_solution_print.print("Operating Point Solution: ");
}

void printOperatingPointWithNames(const arma::vec &op_solution, const Circuitmap &map){
    std::cout << "Operating Bias Point Solution:" << std::endl;
    
    // Find maximum label width for alignment
    size_t max_width = 0;
    for (const auto& node : map.map_nodes) {
        size_t label_width = 3 + node.first.length(); // "V(" + node_name + ")"
        max_width = std::max(max_width, label_width);
    }
    for (const auto& current : map.map_branch_currents) {
        size_t label_width = 3 + current.first.length(); // "I(" + source_name + ")"
        max_width = std::max(max_width, label_width);
    }
    
    // Print node voltages
    for (const auto& node : map.map_nodes) {
        std::string label = "V(" + node.first + ")";
        double voltage = op_solution[node.second - 1]; // node_id - 1 for matrix index
        std::cout << std::left << std::setw(max_width) << label << " " << voltage << std::endl;
    }
    
    // Print branch currents
    for (const auto& current : map.map_branch_currents) {
        std::string label = "I(" + current.first + ")";
        double current_val = op_solution[current.second]; // matrix index directly
        std::cout << std::left << std::setw(max_width) << label << " " << current_val << std::endl;
    }
}

// Only for .op command
OPResult OP_ops(CKTcircuit &ckt, const Modelmap &modmap, bool non_linear) {
    OPResult result;

    // Run the operating point analysis
    result.solution = OperatingPointAnalysis(ckt, modmap, non_linear);

    // Extract MOSFET data
    result.mosfet_data = extractMosfetData(ckt);

    return result;
}
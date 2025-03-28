#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <variant>
#include <map>
#include <memory>

#include "global.hpp"
#include "matrix.hpp"
#include "circuit_parser.hpp"
#include "device.hpp"
#include "map.hpp"

struct CKTcircuit
{
    int pulse_num{};                                // Number of pulse voltages
    // uword is a typedef for an unsigned integer type; it is used for matrix indices as well as all internal counters and loops

    std::vector<CircuitElement> CKTelements;        // Vector of circuit elements
    int external_nodes{};                           // Number of external nodes excluding ground and nodes inside the MOSFETs
    // int external_mosfets{};                      // Number of standalone mosfets (excluding mosfets from ring oscillator)
    int no_of_mosfets{};                            // Total number of MOSFETs
    int no_of_V_sources{};                          // Total number of voltage sources
    int T_nodes{};                                  // Total number of nodes excluding ground (Used to create the initial matrix in CKTsetup)

    std::shared_ptr<DenseMatrix> cktdematrix;       // Dense matrix struct
    std::vector<Capacitor> C_list;                  // Vector of capacitors (including the parasitic capacitance of the MOSFETs)
    Circuitmap map;                                 // Circuit and Model Map struct
    bool ckt_loaded{};                              // To check if the circuit is loaded or not

    double CKTtemp{};                               // Actual temperature of CKT, initialzed to 300.15 K 
    const double CKTnomTemp = 300.15;               // Reference temperature 300.15 K
};

void CKTinstanceSetup(CKTcircuit &ckt, const Modelmap &modmap){
    // 1. Setup the model shared pointer for the instances
    // 2. Setup the instance parameters
    for (auto &element : ckt.CKTelements)
    {
        std::visit([&](auto &&arg) {
            if constexpr (std::is_same_v<std::decay_t<decltype(arg)>, NMOS> || std::is_same_v<std::decay_t<decltype(arg)>, PMOS>) {

                if(arg.modelType == MosfetModelType::BSIM4V82){
                    auto bsim_iter = modmap.bsim4Models.find(arg.modelName);
                    if(bsim_iter != modmap.bsim4Models.end()){
                        // Don't need to set the model pointer again, it's already set in the parser
                        // arg.bsim4v82Instance.BSIM4modPtr = bsim_iter->second;
                        bsim4::instanceSetup(*arg.bsim4v82Instance.BSIM4modPtr, arg.bsim4v82Instance);
                    }
                    else{
                        std::cerr << "Error: BSIM4 model not found in model map for: " << arg.modelName << std::endl;
                        exit(1);
                    }
                }
            }
        }, element.element);
    }
}

void CKTsetup(CKTcircuit &ckt, const CircuitParser &parser, std::shared_ptr<DenseMatrix> denseMatrixPtr, const Modelmap &modmap)
{
    // Careful! getCircuitElements function is const, so it can't be used to modify the elements vector
    // ckt.elements = parser.getCircuitElements();
    ckt.CKTelements = parser.elements;
    ckt.external_nodes = getMaxNode(ckt.CKTelements);
    ckt.no_of_mosfets = parser.num_mosfets;
    ckt.T_nodes = ckt.external_nodes + 3 * ckt.no_of_mosfets;
    // ckt.T_nodes = ckt.external_nodes;
    ckt.CKTtemp = 300.15;   // Initial temperature of the circuit

    // Setup the instances in the circuit (only bsim4)
    if(!modmap.bsim4Models.empty()){
        CKTinstanceSetup(ckt, modmap);
    }

    // Size of matrix
    ckt.cktdematrix = denseMatrixPtr;
    ckt.cktdematrix->Maxi = ckt.T_nodes;
    ckt.cktdematrix->Maxj = ckt.cktdematrix->Maxi;
    ckt.cktdematrix->LHS = arma::zeros(ckt.cktdematrix->Maxi, ckt.cktdematrix->Maxj);    // LHS matrix
    ckt.cktdematrix->RHS = arma::zeros(ckt.cktdematrix->Maxi, 1);                        // RHS matrix
}

void CKTload(CKTcircuit &ckt)
{
    // ASSIGNING THE STAMPS TO THE LHS AND RHS MATRICES

    for (auto &element : ckt.CKTelements)
    {
        std::visit([&](auto &&arg)
                   {
                       if constexpr (std::is_same_v<std::decay_t<decltype(arg)>, Resistor>)
                       {
                            // std::cout << "R Element ID: " << arg.id << ", Node Pos: " << arg.nodePos << ", Node Neg: " << arg.nodeNeg << ", value: "<< arg.value << std::endl;

                           if (arg.value == 0)
                           {
                               arg.value = 1e-3;
                           }
                           R_assigner(arg.nodePos, arg.nodeNeg, 1 / arg.value, ckt.cktdematrix->LHS);
                       }
                       else if constexpr (std::is_same_v<std::decay_t<decltype(arg)>, VoltageSource>)
                       {
                           // std::cout << "VS Element ID: " << arg.id << ", Node Pos: " << arg.nodePos << ", Node Neg: " << arg.nodeNeg << ", value: "<< arg.value << std::endl;

                           Vs_assigner(arg.nodePos, arg.nodeNeg, arg.value, ckt.cktdematrix->LHS, ckt.cktdematrix->RHS);

                           ckt.no_of_V_sources++;
                       }
                       else if constexpr (std::is_same_v<std::decay_t<decltype(arg)>, Pulsevoltage>)
                       {
                           // std::cout << "VPulse Element ID: " << arg.id << ", Node Pos: " << arg.nodePos << ", Node Neg: " << arg.nodeNeg << std::endl;

                           ckt.no_of_V_sources++;

                           ckt.pulse_num++;

                           arg.RHS_locate = V_pulse_assigner(arg.nodePos, arg.nodeNeg, arg.V1, ckt.cktdematrix->LHS, ckt.cktdematrix->RHS);
                       }
                       else if constexpr (std::is_same_v<std::decay_t<decltype(arg)>, CurrentSource>)
                       {
                           // std::cout << "I Element ID: " << arg.id << ", Node Pos: " << arg.nodePos << ", Node Neg: " << arg.nodeNeg << ", value: "<< arg.value << std::endl;

                           Is_assigner(arg.nodePos, arg.nodeNeg, arg.value, ckt.cktdematrix->RHS);
                       }
                       else if constexpr (std::is_same_v<std::decay_t<decltype(arg)>, Capacitor>)
                       {
                           // std::cout << "C Element ID: " << arg.id << ", Node Pos: " << arg.nodePos << ", Node Neg: " << arg.nodeNeg << ", value: "<< arg.value << std::endl;

                           ckt.C_list.push_back(arg);
                       }
                       else if constexpr (std::is_same_v<std::decay_t<decltype(arg)>, Diode>)
                       {
                           // std::cout << "D Element ID: " << arg.id << ", Node Pos: " << arg.nodePos << ", Node Neg: " << arg.nodeNeg << ", value: "<< arg.Is << std::endl;
                       }
                       else if constexpr (std::is_same_v<std::decay_t<decltype(arg)>, NMOS>)
                       {

                       }
                       else if constexpr (std::is_same_v<std::decay_t<decltype(arg)>, PMOS>)
                       {
                        //    std::cout << "PMOS Element ID: " << arg.id << ", Node VD: " << arg.node_vd << ", Node VG: " << arg.node_vg << ", Node VS: " << arg.node_vs << ", Node VB: " << arg.node_vb << std::endl;
                       }
                       else if constexpr (std::is_same_v<std::decay_t<decltype(arg)>, VCCS>)
                       {
                           // std::cout << "G Element ID: " << arg.id << ", Node X: " << arg.node_x << ", Node Y: " << arg.node_y << ", Node CX: " << arg.node_cx << ", Node CY: " << arg.node_cy << ", value: " << arg.value << std::endl;

                           VCCS_assigner(arg.node_x, arg.node_y, arg.node_cx, arg.node_cy, arg.value, ckt.cktdematrix->LHS);
                       } },
                   element.element);
    }
}
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
#include "bsim4v82/bsim4v82temp.hpp"
#include "SPICEcompatible.hpp"

struct CKTcircuit
{
    int pulse_num{};                                // Number of pulse voltages
    // uword is a typedef for an unsigned integer type; it is used for matrix indices as well as all internal counters and loops

    CircuitElements CKTelements;                    // Vector of circuit elements
    int external_nodes{};                           // Number of external nodes excluding ground and nodes inside the MOSFETs
    int internal_nodes{};                           // Number of internal nodes (nodes inside the MOSFETs)
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
    double CKTgmin = 1.0e-20;                               // Gmin value

    SPICECompatible spiceCompatible;                // SPICE-compatible (cktmode)
};

void CKTinstanceSetup(CKTcircuit &ckt, const Modelmap &modmap){
    // Setup the instance parameters
    // Setup the mosfet instance parameters
    // Don't need to set the model pointer again, it's already set in the parser
    // arg.bsim4v82Instance.BSIM4modPtr = bsim_iter->second;

    for (auto &nmos : ckt.CKTelements.nmos){
        if(nmos.modelType == MosfetModelType::BSIM4V82){
            bsim4::instanceSetup(*nmos.bsim4v82Instance.BSIM4modPtr, nmos.bsim4v82Instance);
            bsim4::instanceTemp(nmos.bsim4v82Instance,*nmos.bsim4v82Instance.BSIM4modPtr);
        }
    }
    for (auto &pmos : ckt.CKTelements.pmos){
        if(pmos.modelType == MosfetModelType::BSIM4V82){
            bsim4::instanceSetup(*pmos.bsim4v82Instance.BSIM4modPtr, pmos.bsim4v82Instance);
            bsim4::instanceTemp(pmos.bsim4v82Instance,*pmos.bsim4v82Instance.BSIM4modPtr);
        }
    }
}

void CKTsetup(CKTcircuit &ckt, const CircuitParser &parser, std::shared_ptr<DenseMatrix> denseMatrixPtr, const Modelmap &modmap)
{
    // Careful! getCircuitElements function is const, so it can't be used to modify the elements vector
    // ckt.elements = parser.getCircuitElements();
    ckt.CKTelements = parser.elements;
    ckt.external_nodes = getMaxNode(ckt.CKTelements);
    ckt.internal_nodes = getInternalMosfetNodes(ckt.CKTelements);
    ckt.no_of_mosfets = parser.num_mosfets;
    // ckt.T_nodes = ckt.external_nodes + 3 * ckt.no_of_mosfets;
    ckt.T_nodes = ckt.external_nodes + ckt.internal_nodes;  // Total number of nodes excluding ground
    ckt.CKTtemp = 300.15;                                   // Initial temperature of the circuit
    ckt.spiceCompatible.setMode(0);                         // Initialize the CKTmode to 0

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

    for (const auto &vol : ckt.CKTelements.voltageSources)
    {
        Vs_assigner(vol.nodePos, vol.nodeNeg, vol.value, ckt.cktdematrix->LHS, ckt.cktdematrix->RHS);
        ckt.no_of_V_sources++;
    }
    for (const auto &cur : ckt.CKTelements.currentSources)
    {
        Is_assigner(cur.nodePos, cur.nodeNeg, cur.value, ckt.cktdematrix->RHS);
    }
    for (auto &res : ckt.CKTelements.resistors)
    {
        if (res.value < 1.0e-3)
        {
            res.value = 1.0e-3;
        }
        R_assigner(res.nodePos, res.nodeNeg, 1.0 / res.value, ckt.cktdematrix->LHS);
    }
    for (const auto &cap : ckt.CKTelements.capacitors)
    {
        ckt.C_list.push_back(cap);
    }
    for (auto &pulse : ckt.CKTelements.pulseVoltages)
    {
        ckt.no_of_V_sources++;
        ckt.pulse_num++;
        pulse.RHS_locate = V_pulse_assigner(pulse.nodePos, pulse.nodeNeg, pulse.V1, ckt.cktdematrix->LHS, ckt.cktdematrix->RHS);
    }
    for (auto &vccs : ckt.CKTelements.vccs)
    {
        VCCS_assigner(vccs.node_x, vccs.node_y, vccs.node_cx, vccs.node_cy, vccs.value, ckt.cktdematrix->LHS);
    }
    
}

void updateDeviceState(CKTcircuit &ckt){
    for (auto &element : ckt.CKTelements)
    {
        std::visit([&](auto &&arg)
        {
            if constexpr (std::is_same_v<std::decay_t<decltype(arg)>, NMOS> ||
                                std::is_same_v<std::decay_t<decltype(arg)>, PMOS>)
            {    
                if (arg.modelType == MosfetModelType::BSIM4V82)
                {
                    updateState1(arg.bsim4v82Instance);
                }
            }
        },
        element.element);
    }
}
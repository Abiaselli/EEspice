#pragma once
#include "CKT.hpp"

#include "global.hpp"
#include "matrix.hpp"
#include "circuit_parser.hpp"
#include "device.hpp"
#include "map.hpp"
#include "bsim4v82/bsim4v82setup.hpp"
#include "bsim4v82/bsim4v82temp.hpp"
#include "bsim4v82/bsim4v82stamp.hpp"
#include "bsim4v82/bsim4v82paser.hpp"

bool CKTisNonLinear(const CircuitElements &elements)
{
    bool non_linear = false;
    if(!elements.nmos.empty() || !elements.pmos.empty() || !elements.diodes.empty()){
        non_linear = true;
    }
    return non_linear;
}

// Check if the circuit has any AC sources
bool CKTcheckAC(const CircuitElements &elements){
    bool ACsource = false;
    for (const auto &src : elements.voltageSources) {
        if (src.amplitude >= 0) {
            ACsource = true;
            break;
        }
    }
    return ACsource;
}

void CKTinstanceSetup(CKTcircuit &ckt, const Modelmap &modmap){
    // Create and Setup the instance parameters
    // Setup the mosfet instance parameters

    for (auto &nmos : ckt.CKTelements.nmos){
        if(nmos.modelType == MosfetModelType::BSIM4V82){
            const auto &iter_bsim4 = modmap.bsim4Models.find(nmos.modelName);
            nmos.bsim4v82Instance = bsim4::paserBSIM4instance(nmos.id_str, iter_bsim4->second, nmos.node_vd, nmos.node_vg, nmos.node_vs, nmos.node_vb, nmos.W, nmos.L);
            bsim4::instanceSetup(*nmos.bsim4v82Instance.BSIM4modPtr, nmos.bsim4v82Instance);
            bsim4::instanceTemp(nmos.bsim4v82Instance,*nmos.bsim4v82Instance.BSIM4modPtr);
            if (nmos.bsim4v82Instance.BSIM4trnqsMod){
                ckt.num_of_states++; // BSIM4qcdump
            }
            if(nmos.bsim4v82Instance.BSIM4rbodyMod){
                ckt.num_of_states++; // BSIM4qbs
            }
            if(nmos.bsim4v82Instance.BSIM4rgateMod == 3){
                ckt.num_of_states++; // BSIM4qgmid
            }
            ckt.num_of_states += 3; // BSIM4qb, BSIM4qg, BSIM4qd
        }
    }
    for (auto &pmos : ckt.CKTelements.pmos){
        if(pmos.modelType == MosfetModelType::BSIM4V82){
            const auto &iter_bsim4 = modmap.bsim4Models.find(pmos.modelName);
            pmos.bsim4v82Instance = bsim4::paserBSIM4instance(pmos.id_str, iter_bsim4->second, pmos.node_vd, pmos.node_vg, pmos.node_vs, pmos.node_vb, pmos.W, pmos.L);
            bsim4::instanceSetup(*pmos.bsim4v82Instance.BSIM4modPtr, pmos.bsim4v82Instance);
            bsim4::instanceTemp(pmos.bsim4v82Instance,*pmos.bsim4v82Instance.BSIM4modPtr);
            if (pmos.bsim4v82Instance.BSIM4trnqsMod){
                ckt.num_of_states++; // BSIM4qcdump
            }
            if(pmos.bsim4v82Instance.BSIM4rbodyMod){
                ckt.num_of_states++; // BSIM4qbs
            }
            if(pmos.bsim4v82Instance.BSIM4rgateMod == 3){
                ckt.num_of_states++; // BSIM4qgmid
            }
            ckt.num_of_states += 3; // BSIM4qb, BSIM4qg, BSIM4qd
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
    // ckt.T_nodes = ckt.external_nodes + 3 * ckt.no_of_mosfets;
    ckt.T_nodes = ckt.external_nodes + ckt.internal_nodes;  // Total number of nodes excluding ground
    ckt.CKTtemp = 300.15;                                   // Initial temperature of the circuit
    ckt.spiceCompatible.setMode(0);                         // Initialize the CKTmode to 0
    ckt.CKTintegrateMethod = BACKWARD_EULER;                // Set the integration method to Backward Euler
    ckt.CKTorder = 1;                                       // Set the order of the integration method to 1
    ckt.CKTag.fill(0.0);


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
    ckt.cktdematrix->LHS_cx = arma::cx_mat(ckt.cktdematrix->Maxi, ckt.cktdematrix->Maxj, arma::fill::zeros); // Complex LHS matrix for AC analysis
    ckt.cktdematrix->RHS_cx = arma::cx_mat(ckt.cktdematrix->Maxi, 1, arma::fill::zeros);                     // Complex RHS matrix for AC analysis
}

void CKTload(CKTcircuit &ckt)
{
    // ASSIGNING THE STAMPS TO THE LHS AND RHS MATRICES

    for (const auto &vol : ckt.CKTelements.voltageSources)
    {
        Vs_assigner(vol.nodePos, vol.nodeNeg, vol.value, ckt.cktdematrix->LHS, ckt.cktdematrix->RHS);
        ckt.map.map_branch_currents.insert({vol.id_str, ckt.cktdematrix->RHS.n_rows - 1}); // Store the branch current index in the map
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
        ckt.num_of_states++;
    }
    for (auto &pulse : ckt.CKTelements.pulseVoltages)
    {
        ckt.pulse_num++;
        pulse.RHS_locate = V_pulse_assigner(pulse.nodePos, pulse.nodeNeg, pulse.V1, ckt.cktdematrix->LHS, ckt.cktdematrix->RHS);
        ckt.map.map_branch_currents.insert({pulse.id_str, ckt.cktdematrix->RHS.n_rows - 1}); // Store the branch current index in the map
    }
    for (auto &vccs : ckt.CKTelements.vccs)
    {
        VCCS_assigner(vccs.node_x, vccs.node_y, vccs.node_cx, vccs.node_cy, vccs.value, ckt.cktdematrix->LHS);
    }
    for (auto &vcvs : ckt.CKTelements.vcvs){
        VCVS_assigner(vcvs.node_x, vcvs.node_y, vcvs.node_cx, vcvs.node_cy, vcvs.value, ckt.cktdematrix->LHS, ckt.cktdematrix->RHS);
        ckt.map.map_branch_currents.insert({vcvs.id_str, ckt.cktdematrix->RHS.n_rows - 1}); // Store the branch current index in the map
    }
}
void CKTloadAC(CKTcircuit &ckt){
    // For safety, clear the branch current map before load any linear devices
    ckt.map.map_branch_currents.clear(); 
    // ASSIGNING THE STAMPS TO THE LHS AND RHS MATRICES FOR AC ANALYSIS
    // Only the components that change the size of MNA have their own ACassigner function.
    for (const auto &vol : ckt.CKTelements.voltageSources)
    {
        Vs_ACassigner(vol.nodePos, vol.nodeNeg, vol.acReal, vol.acImag, ckt.cktdematrix->LHS_cx, ckt.cktdematrix->RHS_cx);
        ckt.map.map_branch_currents.insert({vol.id_str, ckt.cktdematrix->RHS_cx.n_rows - 1});
    }
    for (const auto &cur : ckt.CKTelements.currentSources)
    {
       // TODO: AC current source
    }
    for (auto &res : ckt.CKTelements.resistors)
    {
        if (res.value < 1.0e-3)
        {
            res.value = 1.0e-3;
        }
        arma::mat LHS_real = arma::real(ckt.cktdematrix->LHS_cx);
        R_assigner(res.nodePos, res.nodeNeg, 1.0 / res.value, LHS_real);
        ckt.cktdematrix->LHS_cx.set_real(LHS_real);
    }
    // All sources with no acparameter specified are disabled!!!
    // i.e., voltage sources are shorted and current sources removed from the circuit
    constexpr double Shorted = 0.0;
    for (auto &pulse : ckt.CKTelements.pulseVoltages){
        Vs_ACassigner(pulse.nodePos, pulse.nodeNeg, Shorted, Shorted, ckt.cktdematrix->LHS_cx, ckt.cktdematrix->RHS_cx);
        pulse.RHS_locate = ckt.cktdematrix->RHS_cx.n_rows;
        ckt.map.map_branch_currents.insert({pulse.id_str, ckt.cktdematrix->RHS_cx.n_rows - 1});
    }
    for (auto &vccs : ckt.CKTelements.vccs)
    {
        arma::mat LHS_real = arma::real(ckt.cktdematrix->LHS_cx);
        VCCS_assigner(vccs.node_x, vccs.node_y, vccs.node_cx, vccs.node_cy, vccs.value, LHS_real);  // Doesn't change the size of MNA
        ckt.cktdematrix->LHS_cx.set_real(LHS_real);
    }
    for (auto &vcvs : ckt.CKTelements.vcvs){
        VCVS_ACassigner(vcvs.node_x, vcvs.node_y, vcvs.node_cx, vcvs.node_cy, vcvs.value, ckt.cktdematrix->LHS_cx, ckt.cktdematrix->RHS_cx);
        ckt.map.map_branch_currents.insert({vcvs.id_str, ckt.cktdematrix->RHS_cx.n_rows - 1});
    }
}

void updateDeviceState(CKTcircuit &ckt){

    for (auto &nmos : ckt.CKTelements.nmos)
    {
        if (nmos.modelType == MosfetModelType::BSIM4V82)
        {
            bsim4::updateState1(nmos.bsim4v82Instance);
        }
    }
    for (auto &pmos : ckt.CKTelements.pmos)
    {
        if (pmos.modelType == MosfetModelType::BSIM4V82)
        {
            bsim4::updateState1(pmos.bsim4v82Instance);
        }
    }
}

// double CKTterr(double qcap, double ccap, double timestep){

// }
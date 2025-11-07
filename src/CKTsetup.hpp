#pragma once
#include "CKT.hpp"

#include "global.hpp"
#include "matrix.hpp"
#include "circuit_parser.hpp"
#include "device.hpp"
#include "map.hpp"
#include "CKTmkVolt.hpp"
#include "bsim4v82/bsim4v82setup.hpp"
#include "bsim4v82/bsim4v82temp.hpp"
#include "bsim4v82/bsim4v82load/bsim4v82load.hpp"
#include "bsim4v82/bsim4v82paser.hpp"

bool CKTisNonLinear(const CircuitElements &elements)
{
    bool non_linear = false;
    if(!elements.nmos.empty() || 
        !elements.pmos.empty() || 
        !elements.bsim4.empty() ||
        !elements.diodes.empty())
    {
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
        CKTmkVolt(ckt, nmos.id_str + "#drain");
        CKTmkVolt(ckt, nmos.id_str + "#gate");
        CKTmkVolt(ckt, nmos.id_str + "#source"); 
    }
    for (auto &pmos : ckt.CKTelements.pmos){
        CKTmkVolt(ckt, pmos.id_str + "#drain");
        CKTmkVolt(ckt, pmos.id_str + "#gate");
        CKTmkVolt(ckt, pmos.id_str + "#source");
    }
    for (auto &bsim4 : ckt.CKTelements.bsim4){
        const auto &iter_bsim4 = modmap.bsim4Models.find(bsim4.modelName);
        bsim4.bsim4v82Instance = bsim4::paserBSIM4instance(bsim4.id_str, iter_bsim4->second, bsim4.node_vd, bsim4.node_vg, bsim4.node_vs, bsim4.node_vb, bsim4.W, bsim4.L);
        bsim4::instanceSetup(*bsim4.bsim4v82Instance.BSIM4modPtr, bsim4.bsim4v82Instance, ckt);
        bsim4::instanceTemp(bsim4.bsim4v82Instance,*bsim4.bsim4v82Instance.BSIM4modPtr);
        if (bsim4.bsim4v82Instance.BSIM4trnqsMod){
            ckt.num_of_states++; // BSIM4qcdump
        }
        if(bsim4.bsim4v82Instance.BSIM4rbodyMod){
            ckt.num_of_states++; // BSIM4qbs
        }
        if(bsim4.bsim4v82Instance.BSIM4rgateMod == 3){
            ckt.num_of_states++; // BSIM4qgmid
        }
        ckt.num_of_states += 3; // BSIM4qb, BSIM4qg, BSIM4qd
    }
    for (auto &sin : ckt.CKTelements.sinVoltages){
        sin.phase_rad = sin.phase * M_PI / 180.0;   // Convert phase to radians
        if (sin.freq == 0.0){
            sin.freq = 1.0 / ckt.CKTfinalTime;      // Default frequency
        }
    }
}

void CKTloadStatistics(CKTcircuit &ckt){
    ckt.sim_stats.num_threads = ckt.num_threads;
    ckt.sim_stats.num_colors = ckt.b4coloring.getNumColors();
}

void CKTsetup(CKTcircuit &ckt, const CircuitParser &parser, std::shared_ptr<DenseMatrix> denseMatrixPtr, const Modelmap &modmap)
{
    // Careful! getCircuitElements function is const, so it can't be used to modify the elements vector
    // ckt.elements = parser.getCircuitElements();
    ckt.CKTelements = parser.elements;
    ckt.CKTtemp = 300.15;                                   // Initial temperature of the circuit
    ckt.CKTfinalTime = parser.double_t_end;                 // Final time for simulation
    ckt.spiceCompatible.setMode(0);                         // Initialize the CKTmode to 0
    ckt.CKTintegrateMethod = BACKWARD_EULER;                // Set the integration method to Backward Euler
    ckt.CKTorder = 1;                                       // Set the order of the integration method to 1
    ckt.CKTag.fill(0.0);
    ckt.is_batch = parser.is_batch;                         // Set the batch mode
    ckt.CKTmultithreaded = parser.multithreaded;            // Set multithreading mode
    ckt.num_threads = parser.num_threads;                   // Set number of threads for multithreading
    omp_set_num_threads(ckt.num_threads);                   // Set number of OpenMP threads
    ckt.b4coloring.computeColoring(ckt.CKTelements.bsim4);  // Compute coloring for BSIM4 instances

    // Setup the instances in the circuit and create internal nodes
    if(!modmap.bsim4Models.empty()){
        CKTinstanceSetup(ckt, modmap);
    }

    // Get the external and internal nodes
    // ckt.external_nodes = getMaxNode(ckt.CKTelements);
    // ckt.internal_nodes = getInternalMosfetNodes(ckt.CKTelements);
    ckt.external_nodes = ckt.map.map_nodes.size();
    ckt.internal_nodes = ckt.map.map_internal_nodes.size();
    ckt.T_nodes = ckt.external_nodes + ckt.internal_nodes;  // Total number of nodes excluding ground

    // Size of matrix
    ckt.cktdematrix = denseMatrixPtr;
    ckt.cktdematrix->Maxi = ckt.T_nodes;
    ckt.cktdematrix->Maxj = ckt.cktdematrix->Maxi;
    ckt.cktdematrix->LHS = arma::zeros(ckt.cktdematrix->Maxi, ckt.cktdematrix->Maxj);    // LHS matrix
    ckt.cktdematrix->RHS = arma::zeros(ckt.cktdematrix->Maxi, 1);                        // RHS matrix
    ckt.cktdematrix->LHS_cx = arma::cx_mat(ckt.cktdematrix->Maxi, ckt.cktdematrix->Maxj, arma::fill::zeros); // Complex LHS matrix for AC analysis
    ckt.cktdematrix->RHS_cx = arma::cx_mat(ckt.cktdematrix->Maxi, 1, arma::fill::zeros);                     // Complex RHS matrix for AC analysis

    // Load simulation statistics
    CKTloadStatistics(ckt);
}

void CKTload(CKTcircuit &ckt)
{
    // ASSIGNING THE STAMPS TO THE LHS AND RHS MATRICES

    for (const auto &vol : ckt.CKTelements.voltageSources)
    {
        Vs_assigner(vol.nodePos, vol.nodeNeg, vol.value, ckt.cktdematrix->LHS, ckt.cktdematrix->RHS);
        ckt.map.map_branch_currents.emplace(vol.id_str, static_cast<int>(ckt.cktdematrix->RHS.n_rows - 1)); // Store the branch current index in the map
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
        ckt.map.map_branch_currents.emplace(pulse.id_str, static_cast<int>(ckt.cktdematrix->RHS.n_rows - 1)); // Store the branch current index in the map
    }
    for (auto &sin : ckt.CKTelements.sinVoltages){
        sin.RHS_locate = V_sin_assigner(sin.nodePos, sin.nodeNeg, sin.vo, ckt.cktdematrix->LHS, ckt.cktdematrix->RHS);
        ckt.map.map_branch_currents.emplace(sin.id_str, sin.RHS_locate);
    }
    for (auto &vccs : ckt.CKTelements.vccs)
    {
        VCCS_assigner(vccs.node_x, vccs.node_y, vccs.node_cx, vccs.node_cy, vccs.value, ckt.cktdematrix->LHS);
    }
    for (auto &vcvs : ckt.CKTelements.vcvs){
        VCVS_assigner(vcvs.node_x, vcvs.node_y, vcvs.node_cx, vcvs.node_cy, vcvs.value, ckt.cktdematrix->LHS, ckt.cktdematrix->RHS);
        ckt.map.map_branch_currents.emplace(vcvs.id_str, static_cast<int>(ckt.cktdematrix->RHS.n_rows - 1)); // Store the branch current index in the map
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
        ckt.map.map_branch_currents.emplace(vol.id_str, static_cast<int>(ckt.cktdematrix->RHS_cx.n_rows - 1));
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
        pulse.RHS_locate = ckt.cktdematrix->RHS_cx.n_rows - 1;
        ckt.map.map_branch_currents.emplace(pulse.id_str, static_cast<int>(ckt.cktdematrix->RHS_cx.n_rows - 1));
    }
    for (auto &vccs : ckt.CKTelements.vccs)
    {
        arma::mat LHS_real = arma::real(ckt.cktdematrix->LHS_cx);
        VCCS_assigner(vccs.node_x, vccs.node_y, vccs.node_cx, vccs.node_cy, vccs.value, LHS_real);  // Doesn't change the size of MNA
        ckt.cktdematrix->LHS_cx.set_real(LHS_real);
    }
    for (auto &vcvs : ckt.CKTelements.vcvs){
        VCVS_ACassigner(vcvs.node_x, vcvs.node_y, vcvs.node_cx, vcvs.node_cy, vcvs.value, ckt.cktdematrix->LHS_cx, ckt.cktdematrix->RHS_cx);
        ckt.map.map_branch_currents.emplace(vcvs.id_str, static_cast<int>(ckt.cktdematrix->RHS_cx.n_rows - 1));
    }
}

void updateDeviceState(CKTcircuit &ckt){
    for (auto &bsim4 : ckt.CKTelements.bsim4)
    {
        bsim4::updateState1(bsim4.bsim4v82Instance);
    }
}

// double CKTterr(double qcap, double ccap, double timestep){

// }
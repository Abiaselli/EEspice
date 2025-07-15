#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <variant>
#include <armadillo>
#include <cmath>

#include "sim_variables.hpp"
#include "global.hpp"

#define DEBUG_PRINT(x)               \
    if (isDebugMode())               \
    {                                \
        std::cout << x << std::endl; \
    }

// branch extender function
arma::mat branch_ext(arma::mat M, int node_x, int node_y)
{
    int M_cols = M.n_cols;
    int M_rows = M.n_rows;
    arma::mat va = arma::zeros(M_rows, 1);
    if (node_x == 0)
    {
        va.row(node_y - 1).col(0) = 1;
    }
    else if (node_y == 0)
    {
        va.row(node_x - 1).col(0) = 1;
    }
    else
    {
        va.row(node_x - 1).col(0) = -1;
        va.row(node_y - 1).col(0) = 1;
    }

    arma::mat zero_ext = arma::zeros(1, 1);
    arma::mat ha = va.as_row();
    arma::mat haz = arma::join_rows(ha, zero_ext);

    arma::mat M1 = arma::join_rows(M, va);
    arma::mat M2 = arma::join_cols(M1, haz);

    return M2;
}

// create resistor matrix stamp
void R_assigner(int node_x, int node_y, double G, arma::mat &LHS)
{

    // if(G == 0){
    //     std::cerr << "Error: Resistor value cannot be zero" << std::endl;
    //     exit(1);
    // }

    if ((node_x == 0) && (node_y == 0))
    {
        // x = 0;
        return;
    }
    else
    {
        if (node_x == 0)
        {
            LHS.row(node_y - 1).col(node_y - 1) += G;
        }
        else if (node_y == 0)
        {
            LHS.row(node_x - 1).col(node_x - 1) += G;
        }
        else
        {
            LHS.row(node_x - 1).col(node_x - 1) += G;
            LHS.row(node_x - 1).col(node_y - 1) += -G;
            LHS.row(node_y - 1).col(node_x - 1) += -G;
            LHS.row(node_y - 1).col(node_y - 1) += G;
        }
    }

    // RR = RR +1;
}

// Voltage source stamp assigner
void Vs_assigner(int node_x, int node_y, double V_value, arma::mat &LHS, arma::vec &RHS)
{

    arma::vec value(1);
    value = V_value;
    // Extending the branch at the LHS matrix
    LHS = branch_ext(LHS, node_x, node_y);

    // Assigning the value at RHS
    RHS = arma::join_cols(RHS, value);

    // location of voltage value for transient simulation
    // double V_locate1 = RHS.n_rows;

    // return V_locate1;
}
// AC voltage source stamp assigner
void Vs_ACassigner(int node_x, int node_y, double acReal, double acImag, arma::cx_dmat &LHS, arma::cx_dvec &RHS){
    arma::mat real_LHS = arma::real(LHS);
    arma::vec real_RHS = arma::real(RHS);
    arma::mat imag_LHS = arma::imag(LHS);
    arma::vec imag_RHS = arma::imag(RHS);

    Vs_assigner(node_x, node_y, acReal, real_LHS, real_RHS);
    int size = imag_LHS.n_cols;
    imag_LHS.resize(size + 1, size + 1);            // In Armadillo ≥ 10.5, constructors default to zero‐initialization
    imag_RHS.resize(imag_RHS.n_elem + 1);
    imag_RHS(imag_RHS.n_elem - 1) = acImag;
    LHS = arma::cx_dmat(real_LHS, imag_LHS);
    RHS = arma::cx_dvec(real_RHS, imag_RHS);
}

// create a current matrix stamp
void Is_assigner(double node_x, double node_y, double I, arma::vec &RHS)
{

    if ((node_x == 0) && (node_y == 0))
    {
        I = 0;
    }
    else
    {
        if (node_x == 0)
        {
            RHS.row(node_y - 1).col(0) += I;
        }
        else if (node_y == 0)
        {
            RHS.row(node_x - 1).col(0) += -I;
        }
        else
        {
            RHS.row(node_x - 1).col(0) += -I;
            RHS.row(node_y - 1).col(0) += I;
        }
    }
}
// delete the current matrix stamp
void Is_assigner_reverse(double node_x, double node_y, double I, arma::vec &RHS){
    if ((node_x == 0) && (node_y == 0))
    {
        I = 0;
    }
    else
    {
        if (node_x == 0)
        {
            RHS.row(node_y - 1).col(0) -= I;
        }
        else if (node_y == 0)
        {
            RHS.row(node_x - 1).col(0) -= -I;
        }
        else
        {
            RHS.row(node_x - 1).col(0) -= -I;
            RHS.row(node_y - 1).col(0) -= I;
        }
    }
}

// Capacitor stamp assigner with backward euler method
void C_assigner_BE(int node_x, int node_y, double C, double h, arma::mat &LHS, arma::vec &RHS, const arma::vec &pre_solution, int mode)
{
    if (node_x == 0 && node_y == 0)
    {
        return;
    }
    double x{}; // x = C/h
    double vol{};

    if (mode > 0)
    {
        x = C / h;

        if (node_x == 0)
        {
            vol = pre_solution(node_y - 1, 0);
        }
        else if (node_y == 0)
        {
            vol = pre_solution(node_x - 1, 0);
        }
        else
        {
            vol = pre_solution(node_x - 1, 0) - pre_solution(node_y - 1, 0);
        }
    }
    else
    {
        x = 0;
        vol = 0;
    }
    // Matrix stamp for a capacitor on LHS
    R_assigner(node_x, node_y, x, LHS);
    // Matrix stamp for a capacitor on RHS
    Is_assigner(node_x, node_y, -(x * vol), RHS);
}

// Correct Diode stamp assigner (Original Diode_assigner is wrong with the node voltage assignment)
double Diode_assigner(int node_x, int node_y, double Is, double VT, arma::mat &LHS, arma::vec &RHS,
                      const arma::vec &solution)
{
    int maxi = LHS.n_cols;
    int maxj = 1;

    double G_eq = 0;
    double I_eq = 0;
    double val_nodex = 0;
    double val_nodey = 0;

    double Id = 0;
    double vol = 0; // Vd(k)

    // The Companion Model of a Diode is Resistor and Current Source in parallel.
    // Id(k+1) = G_ed * Vd(k+1) + I_eq(k)
    // G_eq = Id' or g'(vd) = Isat/Vt * e^(vd(k)/Vt)
    // I_eq = Id - G_eq * Vd(k)

    if ((node_x == 0 && node_y == 0))
    {
        G_eq = 0;
        I_eq = 0;
        vol = 0;
        Id = Is * (exp(vol / VT) - 1);
        return Id;
    }
    else if (node_x == 0)
    {
        val_nodey = solution(node_y - 1, 0);
        vol = -val_nodey;
        Id = Is * (exp(vol / VT) - 1);
        G_eq = (Is / VT) * (exp(vol / VT));
        I_eq = Id - (G_eq * vol);
    }
    else if (node_y == 0)
    {
        val_nodex = solution(node_x - 1, 0);
        vol = val_nodex;
        Id = Is * (exp(vol / VT) - 1);
        G_eq = (Is / VT) * (exp(vol / VT));
        I_eq = Id - (G_eq * vol);
    }
    else
    {
        val_nodex = solution(node_x - 1, 0);
        val_nodey = solution(node_y - 1, 0);
        vol = val_nodex - val_nodey;
        Id = Is * (exp(vol / VT) - 1);
        G_eq = (Is / VT) * (exp(vol / VT));
        I_eq = Id - (G_eq * vol);
    }

    // Matrix stamp for a diode on RHS
    Is_assigner(node_x, node_y, I_eq, RHS);
    // Matrix stamp for a diode on LHS
    R_assigner(node_x, node_y, G_eq, LHS);

    return Id;
}

int V_pulse_assigner(int node_x, int node_y, double V_value, arma::mat &LHS, arma::vec &RHS)
{

    arma::vec value(1);
    value = V_value;
    // Extending the branch at the LHS matrix
    LHS = branch_ext(LHS, node_x, node_y);

    // Assigning the value at RHS
    RHS = arma::join_cols(RHS, value);

    // location of voltage value for transient simulation
    int Pulse_RHS_locate1 = RHS.n_rows;

    return Pulse_RHS_locate1;
}

double V_pulse_value(double V1, double V2, double t1, double td, double tr, double tf, double tpw, double tper)
{
    double v{};
    t1 = fmod(t1, tper);
    if (tper < tr + tpw + tf)
    {
        std::cerr << "Period is incorrect" << std::endl;
        exit(1);
    }

    if (t1 < 0)
    {
        std::cerr << "Simulation time in pulse voltage it wrong" << std::endl;
        exit(1);
    }
    else if (t1 >= 0 && t1 < td)
    {
        v = V1;
    }
    else if (t1 >= td && t1 < td + tr)
    {
        v = V1 + (V2 - V1) * (t1 - td) / tr;
    }
    else if (t1 >= td + tr && t1 < td + tr + tpw)
    {
        v = V2;
    }
    else if (t1 >= td + tr + tpw && t1 < td + tr + tpw + tf)
    {
        v = V2 + (V1 - V2) * (t1 - td - tr - tpw) / tf;
    }
    else if (t1 >= td + tr + tpw + tf && t1 < tper)
    {
        v = V1;
    }
    else
    {
        std::cerr << "Pulse voltage Error" << std::endl;
        exit(1);
    }
    return v;
}

// assigning the matrix stamps for the VCCS
void VCCS_assigner(int node_x, int node_y, int node_cx, int node_cy, double R, arma::mat &LHS)
{
    int maxi = LHS.n_cols;
    int maxj = LHS.n_rows;

    if (node_x == 0)
    {
        if (node_cx == 0)
        {
            if (node_cy > 0)
            {
                LHS.row(node_y - 1).col(node_cy - 1) += R;
            }
        }
        else if (node_cy == 0)
        {
            if (node_cx > 0)
            {
                LHS.row(node_y - 1).col(node_cx - 1) += -R;
            }
        }
        else
        {
            LHS.row(node_y - 1).col(node_cx - 1) += -R;
            LHS.row(node_y - 1).col(node_cy - 1) += R;
        }
    }
    else if (node_y == 0)
    {
        if (node_cx == 0)
        {
            if (node_cy > 0)
            {
                LHS.row(node_x - 1).col(node_cy - 1) += -R;
            }
        }
        else if (node_cy == 0)
        {
            if (node_cx > 0)
            {
                LHS.row(node_x - 1).col(node_cx - 1) += R;
            }
        }
        else
        {
            LHS.row(node_x - 1).col(node_cx - 1) += R;
            LHS.row(node_x - 1).col(node_cy - 1) += -R;
        }
    }
    else
    {
        if (node_cx == 0)
        {
            if (node_cy > 0)
            {
                LHS.row(node_x - 1).col(node_cy - 1) += -R;
                LHS.row(node_y - 1).col(node_cy - 1) += R;
            }
        }
        else if (node_cy == 0)
        {
            if (node_cx > 0)
            {
                LHS.row(node_x - 1).col(node_cx - 1) += R;
                LHS.row(node_y - 1).col(node_cx - 1) += -R;
            }
        }
        else
        {
            LHS.row(node_x - 1).col(node_cx - 1) += R;
            LHS.row(node_x - 1).col(node_cy - 1) += -R;
            LHS.row(node_y - 1).col(node_cx - 1) += -R;
            LHS.row(node_y - 1).col(node_cy - 1) += R;
        }
    }
}

void VCVS_assigner(int node_x, int node_y, int node_cx, int node_cy, double E, arma::mat &LHS, arma::vec &RHS){
    if (node_x == 0 && node_y == 0){
        std::cerr << "Error: VCVS cannot have both nodes as ground" << std::endl;
        exit(1);
    }
    LHS = LHS.resize(LHS.n_rows + 1, LHS.n_cols + 1); // Resize LHS matrix to accommodate the new row and column
    RHS = RHS.resize(RHS.n_elem + 1); // Resize RHS vector to accommodate the new value
    int row = LHS.n_rows - 1; // New row index for the VCVS
    int col = LHS.n_cols - 1; // New column index for the VCVS

    if (node_x > 0){
        LHS(row, node_x - 1) += 1;  // Last row, column for node_x
        LHS(node_x - 1, col) += 1;  // Last column, row for node_x
    }
    if (node_y > 0){
        LHS(row, node_y - 1) += -1; // Last row, column for node_y
        LHS(node_y - 1, col) += -1; // Last column, row for node_y
    }
    if (node_cx > 0){
        LHS(row, node_cx - 1) += -E; // Last row, column for node_cx
    }
    if (node_cy > 0){
        LHS(row, node_cy - 1) += E;  // Last row, column for node_cy
    }
}

double NMOS_assigner(int number, int node_vd, int node_vg, int node_vs, int node_vb, double W, double L,
                     const arma::vec &solution, int T_nodes, arma::mat &LHS, arma::vec &RHS, const NMOSModel &nmosModel)
{   
    // Get the level 1 parameters
    const auto& params = std::get<NMOSParamLV1>(nmosModel.params);
    
    double vd = 0.0;
    double vg = 0.0;
    double vs = 0.0;
    double vb = 0.0;

    // The settings for the large signal analysis model
    if (node_vd == 0)
    {
        vd = 0.0;
    }
    else
        vd = solution(node_vd - 1, 0);

    if (node_vg == 0)
    {
        vg = 0.0;
    }
    else
        vg = solution(node_vg - 1, 0);

    if (node_vs == 0)
    {
        vs = 0.0;
    }
    else
        vs = solution(node_vs - 1, 0);

    if (node_vb == 0)
    {
        vb = 0.0;
    }
    else
        vb = solution(node_vb - 1, 0);
    // enhancement mode voltages
    double vgs = vg - vs;
    double vds = vd - vs;
    double vbs = vb - vs;
    // depletion mode voltages
    double vbd = vb - vd;
    double vgd = vg - vd;

    double vt = 0.0;
    double id = 0.0;
    double gds = 0.0;
    double gm = 0.0;
    double gmb = 0.0;
    double Leff = L;
    double Beta = (params.mCox) * (W / Leff);
    double I_DSeq = 0.0;

    // Capacitances settings
    double CGS = params.CGSO * W;
    double CGD = params.CGDO * W;
    double CGB = params.CGBO * L;
    // if (vbd <= FC * PB)
    // {
    //     CBD = (CJ * AD) / (pow((1.0 - vbd / PB), MJ)) + (CJSW * PD) / (pow((1.0 - vbd / PB), MJSW));
    // }
    // else
    // {
    //     CBD = ((CJ * AD) * (1.0 - (1.0 + MJ) * FC + MJ * vbd / PB)) / (pow((1.0 - FC), (1.0 + MJ))) + (CJSW * PD) * (1.0 - (1.0 + MJSW) * FC + MJSW * vbd / PB) / (pow((1.0 - FC), (1.0 + MJSW)));
    // }
    // if (vbs <= FC * PB)
    // {
    //     CBS = (CJ * AS) / (pow((1.0 - vbs / PB), MJ)) + (CJSW * PS) / (pow((1.0 - vbs / PB), MJSW));
    // }
    // else
    // {
    //     CBS = ((CJ * AS) * (1.0 - (1.0 + MJ) * FC + MJ * vbs / PB)) / (pow((1.0 - FC), (1.0 + MJ))) + (CJSW * PS) * (1.0 - (1.0 + MJSW) * FC + MJSW * vbs / PB) / (pow((1.0 - FC), (1.0 + MJSW)));
    // }

    R_assigner(node_vd, T_nodes - (3 * number) + 1, 1.0 / params.RD, LHS); // # RD
    R_assigner(node_vg, T_nodes - (3 * number) + 2, 1.0 / params.RG, LHS); // # RG
    R_assigner(T_nodes - (3 * number) + 3, node_vs, 1.0 / params.RS, LHS); // # RS

    vt = params.vt0 + params.gamma * ((sqrt(params.phi - vbs) - sqrt(params.phi))); // # already taking into account the body effect of MOSFETs
    double n_vt = std::abs(vt);

    if ((vds <= (vgs - vt)) && (vgs >= vt))
    { // # the transistor is in linear
        id = Beta * (vgs - vt - vds / 2.0) * vds * (1.0 + params.LAMBDA * vds);
        gds = Beta * (1.0 + params.LAMBDA * vds) * (vgs - vt - vds) + Beta * params.LAMBDA * vds * (vgs - vt - vds / 2.0);
        gm = Beta * vds * (1.0 + params.LAMBDA * vds);
        gmb = gm * params.gamma / (2.0 * sqrt(params.phi - vbs));

        // CGS = 2.0 / 3.0 * params.mCox * (1.0 - pow(vgs - vds - vt, 2) / pow(2.0 * (vgs - vt) - vds, 2)) + params.CGSO * W;
        // CGD = 2.0 / 3.0 * params.mCox * (1.0 - pow(vgs - vt, 2) / pow(2.0 * (vgs - vt) - vds, 2)) + params.CGDO * W;
        // CGB = 0.0 + params.CGBO * L;
    }
    else if ((vds >= (vgs - vt)) && (vgs >= vt))
    { // # the transistor is in saturation
        id = (Beta / 2.0) * pow((vgs - vt), 2) * (1.0 + params.LAMBDA * vds);
        gds = (Beta / 2.0) * params.LAMBDA * pow((vgs - vt), 2);
        gm = Beta * (1 + params.LAMBDA * vds) * (vgs - vt);
        gmb = gm * params.gamma / (2.0 * sqrt(params.phi - vbs));

        // CGS = 2.0 / 3.0 * params.mCox + params.CGSO * W;
        // CGD = params.CGDO * W;
        // CGB = 0.0 + params.CGBO * L;
    }
    else
    {
        id = 0.0;
        gds = 0.0;
        gm = 0.0;
        gmb = 0.0;
    }

    // gds cannot be zero.
    if (gds == 0)
    {
        gds = 1.0e-12;
    }
    if (gm == 0)
    {
        gm = 1.0e-12;
    }

    I_DSeq = id - gds * vds - gm * vgs; // # 10.190 equation

    Is_assigner(T_nodes - (3 * number) + 1, T_nodes - (3 * number) + 3, I_DSeq, RHS);
    R_assigner(T_nodes - (3 * number) + 1, T_nodes - (3 * number) + 3, gds, LHS);                                                           // # assigning gds
    VCCS_assigner(T_nodes - (3 * number) + 1, T_nodes - (3 * number) + 3, T_nodes - (3 * number) + 2, T_nodes - (3 * number) + 3, gm, LHS); // # assigning gm

    return id;
}

double PMOS_assigner(int number, int node_vd, int node_vg, int node_vs, int node_vb, double W, double L,
                     const arma::vec &solution, int T_nodes, arma::mat &LHS, arma::vec &RHS, const PMOSModel &pmosModel)
{   
    // Get the level 1 parameters
    const auto& params = std::get<PMOSParamLV1>(pmosModel.params);

    double vd = 0.0;
    double vg = 0.0;
    double vs = 0.0;
    double vb = 0.0;

    // The settings for the large signal analysis model
    if (node_vd == 0)
    {
        vd = 0.0;
    }
    else
        vd = solution(node_vd - 1, 0);

    if (node_vg == 0)
    {
        vg = 0.0;
    }
    else
        vg = solution(node_vg - 1, 0);

    if (node_vs == 0)
    {
        vs = 0.0;
    }
    else
        vs = solution(node_vs - 1, 0);

    if (node_vb == 0)
    {
        vb = 0.0;
    }
    else
        vb = solution(node_vb - 1, 0);

    // enhancement mode voltages
    double vgs = vg - vs;
    double vds = vd - vs;
    double vbs = vb - vs;
    // enhancement mode voltages
    double vsg = vs - vg;
    double vsd = vs - vd;
    double vsb = vs - vb;
    // depletion mode voltages
    double vdb = vd - vb;
    double vdg = vd - vg;
    double vbd = vb - vd;
    double vgd = vg - vd;

    double id = 0.0;
    double gds = 0.0;
    double gm = 0.0;
    double gmb = 0.0;
    double Leff = L;
    double Beta = (params.mCox) * (W / Leff);
    double vt = 0.0;
    double I_DSeq = 0.0;

    // Capacitances settings
    double CGS = params.CGSO * W;
    double CGD = params.CGDO * W;
    double CGB = params.CGBO * L;
    // if (vdb <= FC * PB)
    // {
    //     CBD = (CJ * AD) / (pow((1 - vdb / PB), MJ)) + (CJSW * PD) / (pow((1 - vdb / PB), MJSW));
    // }
    // else
    // {
    //     CBD = ((CJ * AD) * (1 - (1 + MJ) * FC + MJ * vdb / PB)) / (pow((1 - FC), (1 + MJ))) + (CJSW * PD) * (1 - (1 + MJSW) * FC + MJSW * vdb / PB) / (pow((1 - FC), (1 + MJSW)));
    // }
    // if (vsb <= FC * PB)
    // {
    //     CBS = (CJ * AS) / (pow((1 - vsb / PB), MJ)) + (CJSW * PS) / (pow((1 - vsb / PB), MJSW));
    // }
    // else
    // {
    //     CBS = ((CJ * AS) * (1 - (1 + MJ) * FC + MJ * vsb / PB)) / (pow((1 - FC), (1 + MJ))) + (CJSW * PS) * (1 - (1 + MJSW) * FC + MJSW * vsb / PB) / (pow((1 - FC), (1 + MJSW)));
    // }

    // # the settings for fet model based on the large signal analysis

    R_assigner(T_nodes - (3 * number) + 1, node_vd, 1.0 / params.RD, LHS); // # RD
    R_assigner(node_vg, T_nodes - (3 * number) + 2, 1.0 / params.RG, LHS); // # RG
    R_assigner(node_vs, T_nodes - (3 * number) + 3, 1.0 / params.RS, LHS); // # RS

    vt = params.vt0 - params.gamma * ((sqrt(params.phi - vsb) - sqrt(params.phi))); // already taking into account the body effect of MOSFETs

    double n_vt = std::abs(vt);

    if ((vds >= (vgs - vt)) && (vgs <= vt))
    { // # the transistor is in linear
        id = Beta * (vsg - n_vt - vsd / 2.0) * vsd * (1.0 + params.LAMBDA * vsd);
        gds = Beta * (1.0 + params.LAMBDA * vsd) * (vsg - n_vt - vsd) + Beta * params.LAMBDA * vsd * (vsg - n_vt - vsd / 2.0);
        gm = Beta * vsd * (1.0 + params.LAMBDA * vsd);
        gmb = gm * params.gamma / (2.0 * sqrt(params.phi - vsb));

        // CGS = 2.0 / 3.0 * params.mCox * (1.0 - pow(vsg - vsd - n_vt, 2) / pow(2 * (vsg - n_vt) - vsd, 2)) + params.CGSO * W;
        // CGD = 2.0 / 3.0 * params.mCox * (1.0 - pow(vsg - n_vt, 2) / pow(2 * (vsg - n_vt) - vsd, 2)) + params.CGDO * W;
        // CGB = 0 + params.CGBO * L;
    }
    else if ((vds <= (vgs - vt)) && (vgs <= vt))
    { // # the transistor is in saturation
        id = (Beta / 2.0) * pow((vsg - n_vt), 2) * (1 + params.LAMBDA * vsd);
        gds = (Beta / 2.0) * params.LAMBDA * pow((vsg - n_vt), 2);
        gm = Beta * (1.0 + params.LAMBDA * vsd) * (vsg - n_vt);
        gmb = gm * params.gamma / (2.0 * sqrt(params.phi - vsb));

        // CGS = 2.0 / 3.0 * params.mCox + params.CGSO * W;
        // CGD = 0.0 + params.CGDO * W;
        // CGB = 0.0 + params.CGBO * L;
    }
    else
    { // # the transistor is in cutoff
        id = 0.0;
        gds = 0.0;
        gm = 0.0;
        gmb = 0.0;

        // CGS = 0 + CGSO * W;
        // CGD = 0 + CGDO * W;
        // CGB = mCox + CGBO * L;

        // if(vsg - n_vt <= -phi){
        //     CGS = 0 + CGSO * W;
        //     CGD = 0 + CGDO * W;
        //     CGB = mCox + CGBO * L;
        // }
        // else if(vsg - n_vt > -phi && vsg - n_vt <= -phi/2){
        //     CGS = 0 + CGSO * W;
        //     CGD = 0 + CGDO * W;
        //     CGB = -mCox * ((vsg - n_vt)/phi) + CGBO * L;
        // }
        // else{
        //     CGS = 2.0/3.0 * mCox + 4.0/3.0 * mCox * (vsg - n_vt)/phi + CGSO * W;
        //     CGD = 0.0 + CGDO * W;
        //     CGB = -mCox * ((vsg - n_vt)/phi) + CGBO * L;
        // }
    }

    // gds cannot be zero.
    if (gds == 0)
    {
        gds = 1.0e-12;
    }
    if (gm == 0)
    {
        gm = 1.0e-12;
    }
    I_DSeq = id - gds * vsd - gm * vsg;

    Is_assigner(T_nodes - (3 * number) + 3, T_nodes - (3 * number) + 1, I_DSeq, RHS);
    R_assigner(T_nodes - (3 * number) + 3, T_nodes - (3 * number) + 1, gds, LHS);                                                           // # assigning gds
    VCCS_assigner(T_nodes - (3 * number) + 3, T_nodes - (3 * number) + 1, T_nodes - (3 * number) + 3, T_nodes - (3 * number) + 2, gm, LHS); // # assigning gm

    return id;
}
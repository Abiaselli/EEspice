#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <variant>
#include <armadillo>

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
                      const arma::vec &solution, int mode)
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

double NMOS_assigner(int number, int node_vd, int node_vg, int node_vs, int node_vb, double W, double L, double h,
                     const arma::vec &solution, int T_nodes, arma::mat &LHS, arma::vec &RHS, int mode, std::vector<Capacitor> &C_list, int NR_iteration_counter)
{

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

    double id = 0.0;
    double gds = 0.0;
    double gm = 0.0;
    double gmb = 0.0;
    double Ld = 0.08e-6;
    double Leff = L;
    double kp = 2.0e-5; // default is 2e-5
    double mCox = kp;
    double LAMBDA = 0.1;

    double Beta = (mCox) * (W / Leff);
    double gamma = 0.37;
    double phi = 0.65;
    double vt0 = 0.7; // default is 0
    double vt = 0.0;
    double I_DSeq = 0.0;

    // Capacitances settings
    double CGSO = 4.0e-11;
    double CGDO = 4.0e-11;
    double CGBO = 2.0e-11;
    double CGS = CGSO * W;
    double CGD = CGDO * W;
    double CGB = CGBO * L;
    double CBD = 20.0e-15; // typical value for CBD
    double CBS = 20.0e-15; // typical value for CBS

    double FC = 0.5;        // Coefficient for forward-bias depletion capacitance formula
    double PB = 0.9;        // Bulk junction potential
    double CJ = 0.56e-3;    // Zero bias bulk junction capacitance per unit area
    double MJ = 0.45;       // Bulk junction bottom grading coefficient
    double AD = 200.0e-12;  // Drain area
    double AS = 200.0e-12;  // Source area
    double CJSW = 0.35e-11; // Zero bias bulk junction sidewall capacitance per unit periphery
    double PD = 20.0e-6;    // Drain junction potential
    double PS = 20.0e-6;    // Source junction potential
    double TT = 1.0e-9;     // Transit time
    double MJSW = 0.2;      // Bulk junction sidewall grading coefficient level 1
    double JSSW = 1.0e-9;   // Bulk junction saturation current per meter of sidewall
    double JS = 1.0e-8;     // Bulk junction saturation current per meter of junction perimeter

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

    // # the settings for fet model based on the large signal analysis

    R_assigner(node_vd, T_nodes - (3 * number) + 1, 1.0 / RD, LHS); // # RD
    R_assigner(node_vg, T_nodes - (3 * number) + 2, 1.0 / RG, LHS); // # RG
    R_assigner(T_nodes - (3 * number) + 3, node_vs, 1.0 / RS, LHS); // # RS
    // R_assigner(node_vb ,T_nodes - (4 * number) + 4, 1, LHS, RHS);   // # RB

    // // New model for Diode
    // Diode_assigner(T_nodes - (4 * number) + 3, node_vd, 1e-14, 0.05, LHS, RHS, solution, mode); // # Diode BD
    // Diode_assigner(T_nodes - (4 * number) + 4, node_vs, 1e-14, 0.05, LHS, RHS, solution, mode); // # Diode BS
    // Diode_assigner(node_vb, node_vd, 1e-14, 0.05, LHS, RHS, solution, mode); // # Diode BD
    // Diode_assigner(node_vb, node_vs, 1e-14, 0.05, LHS, RHS, solution, mode); // # Diode BS

    vt = vt0 + gamma * ((sqrt(phi - vbs) - sqrt(phi))); // # already taking into account the body effect of MOSFETs
    double n_vt = std::abs(vt);

    if ((vds <= (vgs - vt)) && (vgs >= vt))
    { // # the transistor is in linear
        id = Beta * (vgs - vt - vds / 2.0) * vds * (1.0 + LAMBDA * vds);
        // gds = Beta/2.0 * LAMBDA * pow(vds,2) + Beta * (vgs - n_vt - vds) * (1 + 2 * LAMBDA * vds);  // From the book: Circuit simulation
        gds = Beta * (1.0 + LAMBDA * vds) * (vgs - vt - vds) + Beta * LAMBDA * vds * (vgs - vt - vds / 2.0);
        gm = Beta * vds * (1.0 + LAMBDA * vds);
        gmb = gm * gamma / (2.0 * sqrt(phi - vbs));

        CGS = 2.0 / 3.0 * mCox * (1.0 - pow(vgs - vds - vt, 2) / pow(2.0 * (vgs - vt) - vds, 2)) + CGSO * W;
        CGD = 2.0 / 3.0 * mCox * (1.0 - pow(vgs - vt, 2) / pow(2.0 * (vgs - vt) - vds, 2)) + CGDO * W;
        CGB = 0.0 + CGBO * L;
    }
    else if ((vds >= (vgs - vt)) && (vgs >= vt))
    { // # the transistor is in saturation
        id = (Beta / 2.0) * pow((vgs - vt), 2) * (1.0 + LAMBDA * vds);
        gds = (Beta / 2.0) * LAMBDA * pow((vgs - vt), 2);
        gm = Beta * (1 + LAMBDA * vds) * (vgs - vt);
        gmb = gm * gamma / (2.0 * sqrt(phi - vbs));

        CGS = 2.0 / 3.0 * mCox + CGSO * W;
        CGD = CGDO * W;
        CGB = 0.0 + CGBO * L;
    }
    else
    {
        id = 0.0;
        gds = 0.0;
        gm = 0.0;
        gmb = 0.0;
        // CGS = 0 + CGSO * W;
        // CGD = 0 + CGDO * W;
        // CGB = mCox + CGBO * L;

        // if(vgs - vt <= -phi){
        //     CGS = 0.0 + CGSO * W;
        //     CGD = 0.0 + CGDO * W;
        //     CGB = mCox + CGBO * L;
        // }else if(vgs - vt > -phi && vgs - vt <= -phi/2){
        //     CGS = 0.0 + CGSO * W;
        //     CGD = 0.0 + CGDO * W;
        //     CGB = -mCox * ((vgs - vt)/phi) + CGBO * L;
        // }else{
        //     CGS = 2.0/3.0 * mCox + 4.0/3.0 * mCox * (vgs - vt)/phi + CGSO * W;
        //     CGD = 0.0 + CGDO * W;
        //     CGB = -mCox * ((vgs - vt)/phi) + CGBO * L;
        // }
    }

    // gds cannot be zero.
    if (gds == 0)
    {
        gds = 1.0e-12;
    }
    // if(id == 0){
    //     id = 1.0e-12;
    // }
    if (gm == 0)
    {
        gm = 1.0e-12;
    }

    I_DSeq = id - gds * vds - gm * vgs; // # 10.190 equation

    // C_assigner_BE(node_vg, T_nodes - (4 * number) + 2, CGD, h, LHS, RHS, solution, mode); // # Capacitor CGD
    // C_assigner_BE(node_vg, T_nodes - (4 * number) + 1, CGS, h, LHS, RHS, solution, mode); // # Capacitor CGS
    // C_assigner_BE(T_nodes - (4 * number) + 3, node_vd, CBD, h, LHS, RHS, solution, mode); // # Capacitor CBD
    // C_assigner_BE(T_nodes - (4 * number) + 4, node_vs, CBS, h, LHS, RHS, solution, mode); // # Capacitor CBS

    // C_assigner_BE(node_vg, node_vd, CGD, h, LHS, RHS, solution, mode); // # Capacitor CGD
    // C_assigner_BE(node_vg, node_vs, CGS, h, LHS, RHS, solution, mode); // # Capacitor CGS
    // if(node_vb != 0 || node_vs != 0){
    //     C_assigner_BE(node_vb, node_vs, CBS, h, LHS, RHS, solution, mode); // # Capacitor CBS
    // }
    // C_assigner_BE(node_vg, node_vb, CGB, h, LHS, RHS, solution, mode); // # Capacitor CGB
    // C_assigner_BE(node_vb, node_vd, CBD, h, LHS, RHS, solution, mode); // # Capacitor CBD

    Is_assigner(T_nodes - (3 * number) + 1, T_nodes - (3 * number) + 3, I_DSeq, RHS);
    R_assigner(T_nodes - (3 * number) + 1, T_nodes - (3 * number) + 3, gds, LHS);                                                           // # assigning gds
    VCCS_assigner(T_nodes - (3 * number) + 1, T_nodes - (3 * number) + 3, T_nodes - (3 * number) + 2, T_nodes - (3 * number) + 3, gm, LHS); // # assigning gm

    // VCCS_assigner(node_vd, node_vs, node_vb, node_vs, gmb, LHS); // assigning gmb
    // std::cout << "gds_NMOS: " << gds << std::endl;
    // std::cout << "gm_NMOS: " << gm << std::endl;
    // std::cout << "id_eq_NMOS: " << I_DSeq << std::endl;
    // std::cout << "id_NMOS: " << id << std::endl;

    DEBUG_PRINT("NMOS gds: " << gds);
    DEBUG_PRINT("NMOS gm: " << gm);
    DEBUG_PRINT("NMOS id_eq: " << I_DSeq);
    DEBUG_PRINT("NMOS id: " << id);

    // if(mode == 0 && NR_iteration_counter == 0){

    // Capacitor c1{0, "M"+std::to_string(number)+".1", node_vg, node_vd, CGD}; // M1.1
    // Capacitor c2{0, "M"+std::to_string(number)+".2", node_vg, node_vs, CGS}; // M1.2
    // if(node_vb != 0 || node_vs != 0){
    //     Capacitor c3{0, "M"+std::to_string(number)+".3", node_vb, node_vs, CBS}; // M1.3
    //     C_list.push_back(c3);
    // }
    // Capacitor c4{0, "M"+std::to_string(number)+".4", node_vg, node_vb, CGB}; // M1.4
    // Capacitor c5{0, "M"+std::to_string(number)+".5", node_vb, node_vd, CBD}; // M1.5

    // C_list.push_back(c1);
    // C_list.push_back(c2);

    // C_list.push_back(c4);
    // C_list.push_back(c5);
    // }
    // else{

    //     for(auto& Cap: C_list){                                      // Update the Capacitance values
    //         if(Cap.name == "M"+std::to_string(number)+".1"){
    //             Cap.value = CGD;
    //         }
    //         else if(Cap.name == "M"+std::to_string(number)+".2"){
    //             Cap.value = CGS;
    //         }
    //         // else if(Cap.name == "M"+std::to_string(number)+".3"){
    //         //     Cap.value = CBS;
    //         // }
    //         else if(Cap.name == "M"+std::to_string(number)+".4"){
    //             Cap.value = CGB;
    //         }
    //         // else if(Cap.name == "M"+std::to_string(number)+".5"){
    //         //     Cap.value = CBD;
    //         // }

    //     }
    // }

    return id;
}

double PMOS_assigner(int number, int node_vd, int node_vg, int node_vs, int node_vb, double W, double L, double h,
                     const arma::vec &solution, int T_nodes, arma::mat &LHS, arma::vec &RHS, int mode, std::vector<Capacitor> &C_list, int NR_iteration_counter)
{
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
    double Ld = 0.0;
    double Leff = L;
    double kp = 2.0e-5; // default is 2e-5
    double mCox = kp;
    double LAMBDA = 0.1;

    double Beta = (mCox) * (W / Leff);
    double gamma = 0.4;
    double phi = 0.8;
    double vt0 = -0.8; // default is 0
    double vt = 0.0;
    double I_DSeq = 0.0;

    // Capacitances settings
    double CGSO = 4.0e-10;
    double CGDO = 4.0e-10;
    double CGBO = 2.0e-10;
    double CGS = CGSO * W;
    double CGD = CGDO * W;
    double CGB = CGBO * L;
    double CBD = 20.0e-15; // typical value for CBD
    double CBS = 20.0e-15; // typical value for CBS

    double FC = 0.5;       // Coefficient for forward-bias depletion capacitance formula
    double PB = 0.8;       // Bulk junction potential
    double CJ = 2.0e-4;    // Zero bias bulk junction capacitance per unit area
    double MJ = 0.5;       // Bulk junction bottom grading coefficient
    double AD = 200.0e-12; // Drain area
    double AS = 200.0e-12; // Source area
    double CJSW = 1.0e-12; // Zero bias bulk junction sidewall capacitance per unit periphery
    double PD = 20.0e-6;   // Drain junction potential
    double PS = 20.0e-6;   // Source junction potential
    double TT = 1.0e-9;    // Transit time
    double MJSW = 0.5;     // Bulk junction sidewall grading coefficient level 1
    double JSSW = 1.0e-9;  // Bulk junction saturation current per meter of sidewall
    double JS = 1.0e-6;    // Bulk junction saturation current per meter of junction perimeter

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

    R_assigner(T_nodes - (3 * number) + 1, node_vd, 1.0 / RD, LHS); // # RD
    R_assigner(node_vg, T_nodes - (3 * number) + 2, 1.0 / RG, LHS); // # RG
    R_assigner(node_vs, T_nodes - (3 * number) + 3, 1.0 / RS, LHS); // # RS
    // R_assigner(T_nodes - (4 * number) + 4, node_vb, 1, LHS, RHS);   // # RB

    // Diode_assigner(T_nodes - (4 * number) + 2, T_nodes - (4 * number) + 3, 1e-14, 0.05, LHS, RHS,  solution, mode); // # Diode DB
    // Diode_assigner(T_nodes - (4 * number) + 4, T_nodes - (4 * number) + 3, 1e-14, 0.05, LHS, RHS,  solution, mode); // # Diode SB
    // Diode_assigner(node_vd, node_vb, 1e-14, 0.05, LHS, RHS, solution, mode); // # Diode BD
    // Diode_assigner(node_vs, node_vb, 1e-14, 0.05, LHS, RHS, solution, mode); // # Diode BS

    vt = vt0 - gamma * ((sqrt(phi - vsb) - sqrt(phi))); // already taking into account the body effect of MOSFETs

    double n_vt = std::abs(vt);

    if ((vds >= (vgs - vt)) && (vgs <= vt))
    { // # the transistor is in linear
        id = Beta * (vsg - n_vt - vsd / 2.0) * vsd * (1.0 + LAMBDA * vsd);
        gds = Beta * (1.0 + LAMBDA * vsd) * (vsg - n_vt - vsd) + Beta * LAMBDA * vsd * (vsg - n_vt - vsd / 2.0);
        gm = Beta * vsd * (1.0 + LAMBDA * vsd);
        gmb = gm * gamma / (2.0 * sqrt(phi - vsb));

        CGS = 2.0 / 3.0 * mCox * (1.0 - pow(vsg - vsd - n_vt, 2) / pow(2 * (vsg - n_vt) - vsd, 2)) + CGSO * W;
        CGD = 2.0 / 3.0 * mCox * (1.0 - pow(vsg - n_vt, 2) / pow(2 * (vsg - n_vt) - vsd, 2)) + CGDO * W;
        CGB = 0 + CGBO * L;
    }
    else if ((vds <= (vgs - vt)) && (vgs <= vt))
    { // # the transistor is in saturation
        id = (Beta / 2.0) * pow((vsg - n_vt), 2) * (1 + LAMBDA * vsd);
        gds = (Beta / 2.0) * LAMBDA * pow((vsg - n_vt), 2);
        gm = Beta * (1.0 + LAMBDA * vsd) * (vsg - n_vt);
        gmb = gm * gamma / (2.0 * sqrt(phi - vsb));

        CGS = 2.0 / 3.0 * mCox + CGSO * W;
        CGD = 0.0 + CGDO * W;
        CGB = 0.0 + CGBO * L;
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

    // C_assigner_BE(T_nodes - (4 * number) + 2, T_nodes - (4 * number) + 3, CBD, h, LHS, RHS, solution, mode); // # Capacitor CBD
    // C_assigner_BE(T_nodes - (4 * number) + 4, T_nodes - (4 * number) + 3, CBS, h, LHS, RHS, solution, mode); // # Capacitor CBS
    // C_assigner_BE(T_nodes - (4 * number) + 1, T_nodes - (4 * number) + 4, CGS, h, LHS, RHS, solution, mode); // # Capacitor GSO
    // C_assigner_BE(T_nodes - (4 * number) + 1, T_nodes - (4 * number) + 3, CGB, h, LHS, RHS, solution, mode); // # Capacitor GBO
    // C_assigner_BE(T_nodes - (4 * number) + 4, T_nodes - (4 * number) + 3, CGD, h, LHS, RHS, solution, mode); // # Capacitor GDO

    // C_assigner_BE(node_vd, node_vg, CGD, h, LHS, RHS, solution, mode); // # Capacitor CGD
    // C_assigner_BE(node_vs, node_vg, CGS, h, LHS, RHS, solution, mode); // # Capacitor CGS
    // C_assigner_BE(node_vb, node_vg, CGB, h, LHS, RHS, solution, mode); // # Capacitor CGB
    // C_assigner_BE(node_vd, node_vb, CBD, h, LHS, RHS, solution, mode); // # Capacitor CBD
    // C_assigner_BE(node_vs, node_vb, CBS, h, LHS, RHS, solution, mode); // # Capacitor CBS

    // gds cannot be zero.
    if (gds == 0)
    {
        gds = 1.0e-12;
    }
    // if(id == 0){
    //     id = 1.0e-12;
    // }
    if (gm == 0)
    {
        gm = 1.0e-12;
    }

    I_DSeq = id - gds * vsd - gm * vsg;

    Is_assigner(T_nodes - (3 * number) + 3, T_nodes - (3 * number) + 1, I_DSeq, RHS);
    // VCCS_assigner(node_vs, node_vd, node_vs, node_vb, gmb, LHS); // assigning gmb
    R_assigner(T_nodes - (3 * number) + 3, T_nodes - (3 * number) + 1, gds, LHS);                                                           // # assigning gds
    VCCS_assigner(T_nodes - (3 * number) + 3, T_nodes - (3 * number) + 1, T_nodes - (3 * number) + 3, T_nodes - (3 * number) + 2, gm, LHS); // # assigning gm

    // std::cout << "gds_PMOS: " << gds << std::endl;
    // std::cout << "gm_PMOS: " << gm << std::endl;
    // std::cout << "id_eq_PMOS: " << I_DSeq << std::endl;
    // std::cout << "id_PMOS: " << id << std::endl;

    DEBUG_PRINT("PMOS gds: " << gds);
    DEBUG_PRINT("PMOS gm: " << gm);
    DEBUG_PRINT("PMOS id_eq: " << I_DSeq);
    DEBUG_PRINT("PMOS id: " << id);

    // if(mode == 0 && NR_iteration_counter == 0){

    //     Capacitor c1{0, "M"+std::to_string(number)+".1", node_vd, node_vg, CGD}; // M1.1
    //     Capacitor c2{0, "M"+std::to_string(number)+".2", node_vs, node_vg, CGS}; // M1.2
    //     Capacitor c3{0, "M"+std::to_string(number)+".3", node_vb, node_vg, CGB}; // M1.3
    //     Capacitor c4{0, "M"+std::to_string(number)+".4", node_vd, node_vb, CBD}; // M1.4
    //     Capacitor c5{0, "M"+std::to_string(number)+".5", node_vs, node_vb, CBS}; // M1.5

    //     C_list.push_back(c1);
    //     C_list.push_back(c2);
    //     C_list.push_back(c3);
    //     C_list.push_back(c4);
    //     C_list.push_back(c5);
    // }
    // else{

    //     for(auto& Cap: C_list){                                      // Update the Capacitance values
    //         if(Cap.name == "M"+std::to_string(number)+".1"){
    //             Cap.value = CGD;
    //         }
    //         else if(Cap.name == "M"+std::to_string(number)+".2"){
    //             Cap.value = CGS;
    //         }
    //         else if(Cap.name == "M"+std::to_string(number)+".3"){
    //             Cap.value = CGB;
    //         }
    //         else if(Cap.name == "M"+std::to_string(number)+".4"){
    //             Cap.value = CBD;
    //         }
    //         else if(Cap.name == "M"+std::to_string(number)+".5"){
    //             Cap.value = CBS;
    //         }

    //     }
    // }

    return id;
}
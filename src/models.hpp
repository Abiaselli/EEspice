#pragma once
#include <string>
#include <variant>
#include <map>
struct NMOSModel;
struct PMOSModel;
struct NMOSParamLV1;
struct PMOSParamLV1;

struct NMOSParamLV1{
    double Ld = 0.08e-6;
    double kp = 2.0e-5; // default is 2e-5
    double mCox = kp;
    double LAMBDA = 0.1;

    // double Beta = (mCox) * (W / Leff);   Leff = L;
    double gamma = 0.37;
    double phi = 0.65;
    double vt0 = 0.7; // default is 0

    // Capacitances settings
    double CGSO = 4.0e-11;
    double CGDO = 4.0e-11;
    double CGBO = 2.0e-11;
    // double CGS = CGSO * W;
    // double CGD = CGDO * W;
    // double CGB = CGBO * L;
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

    double RD = 1.0;
    double RG = 1.0;
    double RS = 1.0;

    void setFromMap(const std::map<std::string, std::string> &kvMap);
};

struct PMOSParamLV1{
    double Ld = 0.0;
    // double Leff = L;
    double kp = 2.0e-5; // default is 2e-5
    double mCox = kp;
    double LAMBDA = 0.1;

    // double Beta = (mCox) * (W / Leff);
    double gamma = 0.4;
    double phi = 0.8;
    double vt0 = -0.8; // default is 0

    // Capacitances settings
    double CGSO = 4.0e-10;
    double CGDO = 4.0e-10;
    double CGBO = 2.0e-10;
    // double CGS = CGSO * W;
    // double CGD = CGDO * W;
    // double CGB = CGBO * L;
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

    double RD = 1.0;
    double RG = 1.0;
    double RS = 1.0;

    void setFromMap(const std::map<std::string, std::string> &kvMap);
};


struct NMOSModel{
    int level = 1;
    std::string version;

    // NMOS parameters
    std::variant<NMOSParamLV1> params = NMOSParamLV1{};
};

struct PMOSModel {
    int level = 1;
    std::string version;

    // PMOS parameters
    std::variant<PMOSParamLV1> params = PMOSParamLV1{};
};
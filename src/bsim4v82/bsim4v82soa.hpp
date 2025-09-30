#pragma once
#include "bsim4v82.hpp"
#include <vector>

struct BSIM4V82soa{

    std::string modelName;

    std::vector<std::string> BSIM4name;
    std::array<std::vector<double>, 29> BSIM4states0;
    std::array<std::vector<double>, 29> BSIM4states1;

    // Node
    enum NodeType{
        D_NODE = 0, G_NODE_EXT, S_NODE, B_NODE, D_NODE_PRIME, G_NODE_PRIME, G_NODE_MID, S_NODE_PRIME,
        B_NODE_PRIME, DB_NODE, SB_NODE, Q_NODE
    };
    std::array<std::vector<uint8_t>, 12> BSIM4nodeValid; // // 0=false, 1=true

    std::vector<int> BSIM4dNode;
    std::vector<int> BSIM4gNodeExt;
    std::vector<int> BSIM4sNode;
    std::vector<int> BSIM4bNode;
    std::vector<int> BSIM4dNodePrime;
    std::vector<int> BSIM4gNodePrime;
    std::vector<int> BSIM4gNodeMid;
    std::vector<int> BSIM4sNodePrime;
    std::vector<int> BSIM4bNodePrime;
    std::vector<int> BSIM4dbNode;
    std::vector<int> BSIM4sbNode;
    std::vector<int> BSIM4qNode;

    std::vector<double> BSIM4ueff;
    std::vector<double> BSIM4thetavth;
    std::vector<double> BSIM4von;
    std::vector<double> BSIM4vdsat;
    std::vector<double> BSIM4cgdo;
    std::vector<double> BSIM4qgdo;
    std::vector<double> BSIM4cgso;
    std::vector<double> BSIM4qgso;
    std::vector<double> BSIM4grbsb;
    std::vector<double> BSIM4grbdb;
    std::vector<double> BSIM4grbpb;
    std::vector<double> BSIM4grbps;
    std::vector<double> BSIM4grbpd;

    std::vector<double> BSIM4vjsmFwd;
    std::vector<double> BSIM4vjsmRev;
    std::vector<double> BSIM4vjdmFwd;
    std::vector<double> BSIM4vjdmRev;
    std::vector<double> BSIM4XExpBVS;
    std::vector<double> BSIM4XExpBVD;
    std::vector<double> BSIM4SslpFwd;
    std::vector<double> BSIM4SslpRev;
    std::vector<double> BSIM4DslpFwd;
    std::vector<double> BSIM4DslpRev;
    std::vector<double> BSIM4IVjsmFwd;
    std::vector<double> BSIM4IVjsmRev;
    std::vector<double> BSIM4IVjdmFwd;
    std::vector<double> BSIM4IVjdmRev;

    std::vector<double> BSIM4grgeltd;
    std::vector<double> BSIM4Pseff;
    std::vector<double> BSIM4Pdeff;
    std::vector<double> BSIM4Aseff;
    std::vector<double> BSIM4Adeff;

    std::vector<double> BSIM4l;
    std::vector<double> BSIM4w;
    std::vector<double> BSIM4drainArea;
    std::vector<double> BSIM4sourceArea;
    std::vector<double> BSIM4drainSquares;
    std::vector<double> BSIM4sourceSquares;
    std::vector<double> BSIM4drainPerimeter;
    std::vector<double> BSIM4sourcePerimeter;
    std::vector<double> BSIM4sourceConductance;
    std::vector<double> BSIM4drainConductance;

    /* stress effect instance param */
    std::vector<double> BSIM4sa;
    std::vector<double> BSIM4sb;
    std::vector<double> BSIM4sd;
    std::vector<double> BSIM4sca;
    std::vector<double> BSIM4scb;
    std::vector<double> BSIM4scc;
    std::vector<double> BSIM4sc;

    std::vector<double> BSIM4rbdb;
    std::vector<double> BSIM4rbsb;
    std::vector<double> BSIM4rbpb;
    std::vector<double> BSIM4rbps;
    std::vector<double> BSIM4rbpd;

    std::vector<double> BSIM4delvto;
    std::vector<double> BSIM4xgw;
    std::vector<double> BSIM4ngcon;

    /* added here to account stress effect instance dependence */

    std::vector<double> BSIM4u0temp;
    std::vector<double> BSIM4vsattemp;
    std::vector<double> BSIM4vth0;
    std::vector<double> BSIM4vfb;
    std::vector<double> BSIM4vfbzb;
    std::vector<double> BSIM4vtfbphi1;
    std::vector<double> BSIM4vtfbphi2;
    std::vector<double> BSIM4k2;
    std::vector<double> BSIM4vbsc;
    std::vector<double> BSIM4k2ox;
    std::vector<double> BSIM4eta0;
    std::vector<double> BSIM4toxp;
    std::vector<double> BSIM4coxp;

    std::vector<double> BSIM4icVDS;
    std::vector<double> BSIM4icVGS;
    std::vector<double> BSIM4icVBS;
    std::vector<double> BSIM4nf;
    std::vector<int> BSIM4off;

    std::vector<int> BSIM4mode;
    std::vector<int> BSIM4trnqsMod;
    std::vector<int> BSIM4acnqsMod;
    std::vector<int> BSIM4rbodyMod;
    std::vector<int> BSIM4rgateMod;
    std::vector<int> BSIM4geoMod;
    std::vector<int> BSIM4rgeoMod;
    std::vector<int> BSIM4min;


    /* OP point */
    std::vector<double> BSIM4Vgsteff;
    std::vector<double> BSIM4vgs_eff;
    std::vector<double> BSIM4vgd_eff;
    std::vector<double> BSIM4dvgs_eff_dvg;
    std::vector<double> BSIM4dvgd_eff_dvg;
    std::vector<double> BSIM4Vdseff;
    std::vector<double> BSIM4nstar;
    std::vector<double> BSIM4Abulk;
    std::vector<double> BSIM4EsatL;
    std::vector<double> BSIM4AbovVgst2Vtm;
    std::vector<double> BSIM4qinv;
    std::vector<double> BSIM4cd;
    std::vector<double> BSIM4cbs;
    std::vector<double> BSIM4cbd;
    std::vector<double> BSIM4csub;
    std::vector<double> BSIM4Igidl;
    std::vector<double> BSIM4Igisl;
    std::vector<double> BSIM4gm;
    std::vector<double> BSIM4gds;
    std::vector<double> BSIM4gmbs;
    std::vector<double> BSIM4gbd;
    std::vector<double> BSIM4gbs;
    std::vector<double> BSIM4noiGd0;   /* tnoiMod=2 (v4.7) */
    std::vector<double> BSIM4Coxeff;

    std::vector<double> BSIM4gbbs;
    std::vector<double> BSIM4gbgs;
    std::vector<double> BSIM4gbds;
    std::vector<double> BSIM4ggidld;
    std::vector<double> BSIM4ggidlg;
    std::vector<double> BSIM4ggidls;
    std::vector<double> BSIM4ggidlb;
    std::vector<double> BSIM4ggisld;
    std::vector<double> BSIM4ggislg;
    std::vector<double> BSIM4ggisls;
    std::vector<double> BSIM4ggislb;

    std::vector<double> BSIM4Igcs;
    std::vector<double> BSIM4gIgcsg;
    std::vector<double> BSIM4gIgcsd;
    std::vector<double> BSIM4gIgcss;
    std::vector<double> BSIM4gIgcsb;
    std::vector<double> BSIM4Igcd;
    std::vector<double> BSIM4gIgcdg;
    std::vector<double> BSIM4gIgcdd;
    std::vector<double> BSIM4gIgcds;
    std::vector<double> BSIM4gIgcdb;

    std::vector<double> BSIM4Igs;
    std::vector<double> BSIM4gIgsg;
    std::vector<double> BSIM4gIgss;
    std::vector<double> BSIM4Igd;
    std::vector<double> BSIM4gIgdg;
    std::vector<double> BSIM4gIgdd;

    std::vector<double> BSIM4Igb;
    std::vector<double> BSIM4gIgbg;
    std::vector<double> BSIM4gIgbd;
    std::vector<double> BSIM4gIgbs;
    std::vector<double> BSIM4gIgbb;

    std::vector<double> BSIM4grdsw;
    std::vector<double> BSIM4IdovVds;
    std::vector<double> BSIM4gcrg;
    std::vector<double> BSIM4gcrgd;
    std::vector<double> BSIM4gcrgg;
    std::vector<double> BSIM4gcrgs;
    std::vector<double> BSIM4gcrgb;

    std::vector<double> BSIM4gstot;
    std::vector<double> BSIM4gstotd;
    std::vector<double> BSIM4gstotg;
    std::vector<double> BSIM4gstots;
    std::vector<double> BSIM4gstotb;

    std::vector<double> BSIM4gdtot;
    std::vector<double> BSIM4gdtotd;
    std::vector<double> BSIM4gdtotg;
    std::vector<double> BSIM4gdtots;
    std::vector<double> BSIM4gdtotb;

    std::vector<double> BSIM4cggb;
    std::vector<double> BSIM4cgdb;
    std::vector<double> BSIM4cgsb;
    std::vector<double> BSIM4cbgb;
    std::vector<double> BSIM4cbdb;
    std::vector<double> BSIM4cbsb;
    std::vector<double> BSIM4cdgb;
    std::vector<double> BSIM4cddb;
    std::vector<double> BSIM4cdsb;
    std::vector<double> BSIM4csgb;
    std::vector<double> BSIM4csdb;
    std::vector<double> BSIM4cssb;
    std::vector<double> BSIM4cgbb;
    std::vector<double> BSIM4cdbb;
    std::vector<double> BSIM4csbb;
    std::vector<double> BSIM4cbbb;
    std::vector<double> BSIM4capbd;
    std::vector<double> BSIM4capbs;

    std::vector<double> BSIM4cqgb;
    std::vector<double> BSIM4cqdb;
    std::vector<double> BSIM4cqsb;
    std::vector<double> BSIM4cqbb;

    std::vector<double> BSIM4qgate;
    std::vector<double> BSIM4qbulk;
    std::vector<double> BSIM4qdrn;
    std::vector<double> BSIM4qsrc;
    std::vector<double> BSIM4qdef;

    std::vector<double> BSIM4qchqs;
    std::vector<double> BSIM4taunet;
    std::vector<double> BSIM4gtau;
    std::vector<double> BSIM4gtg;
    std::vector<double> BSIM4gtd;
    std::vector<double> BSIM4gts;
    std::vector<double> BSIM4gtb;
    std::vector<double> BSIM4SjctTempRevSatCur;
    std::vector<double> BSIM4DjctTempRevSatCur;
    std::vector<double> BSIM4SswTempRevSatCur;
    std::vector<double> BSIM4DswTempRevSatCur;
    std::vector<double> BSIM4SswgTempRevSatCur;
    std::vector<double> BSIM4DswgTempRevSatCur;
};
/*
    bsim4def.h - BSIM4v4.8.2
*/
/* ******************************************************************************
   *  BSIM4 4.8.2 released by Chetan Kumar Dabhi 01/01/2020                     *
   *  BSIM4 Model Equations                                                     *
   ******************************************************************************

   ******************************************************************************
   *  Copyright (c) 2020 University of California                               *
   *                                                                            *
   *  Project Director: Prof. Chenming Hu.                                      *
   *  Current developers: Chetan Kumar Dabhi   (Ph.D. student, IIT Kanpur)      *
   *                      Prof. Yogesh Chauhan (IIT Kanpur)                     *
   *                      Dr. Pragya Kushwaha  (Postdoc, UC Berkeley)           *
   *                      Dr. Avirup Dasgupta  (Postdoc, UC Berkeley)           *
   *                      Ming-Yen Kao         (Ph.D. student, UC Berkeley)     *
   *  Authors: Gary W. Ng, Weidong Liu, Xuemei Xi, Mohan Dunga, Wenwei Yang     *
   *           Ali Niknejad, Chetan Kumar Dabhi, Yogesh Singh Chauhan,          *
   *           Sayeef Salahuddin, Chenming Hu                                   * 
   ******************************************************************************/
#pragma once
#include <string>
#include <memory>
#include <vector>
#include <array>

namespace bsim4{
struct bsim4SizeDependParam;
// Temporary parameters for bsim4v82temp.hpp
struct bsim4v82temp
{
    double tmp = 0.0, tmp1 = 0.0, tmp2 = 0.0, tmp3 = 0.0, Eg = 0.0, Eg0 = 0.0, ni = 0.0, epssub = 0.0;
    double T0 = 0.0, T1 = 0.0, T2 = 0.0, T3 = 0.0, T4 = 0.0, T5 = 0.0, T6 = 0.0, T7 = 0.0, T8 = 0.0, T9 = 0.0, Lnew = 0.0, Wnew = 0.0;
    double delTemp = 0.0, Temp = 0.0, TRatio = 0.0, Inv_L = 0.0, Inv_W = 0.0, Inv_LW = 0.0, Dw = 0.0, Dl = 0.0, Vtm0 = 0.0, Tnom = 0.0;
    double dumPs = 0.0, dumPd = 0.0, dumAs = 0.0, dumAd = 0.0, PowWeffWr = 0.0;
    double DMCGeff = 0.0, DMCIeff = 0.0, DMDGeff = 0.0;
    double Nvtms = 0.0, Nvtmd = 0.0, SourceSatCurrent = 0.0, DrainSatCurrent = 0.0;
    double T10 = 0.0, T11 = 0.0;
    double Inv_saref = 0.0, Inv_sbref = 0.0, Inv_sa = 0.0, Inv_sb = 0.0, rho = 0.0, Ldrn = 0.0, dvth0_lod = 0.0;
    double W_tmp = 0.0, Inv_ODeff = 0.0, OD_offset = 0.0, dk2_lod = 0.0, deta0_lod = 0.0;
    double lnl = 0.0, lnw = 0.0, lnnf = 0.0, rbpbx = 0.0, rbpby = 0.0, rbsbx = 0.0, rbsby = 0.0, rbdbx = 0.0, rbdby = 0.0, bodymode = 0.0;
    double kvsat = 0.0, wlod = 0.0, sceff = 0.0, Wdrn = 0.0;
    double V0 = 0.0, lt1 = 0.0, ltw = 0.0, Theta0 = 0.0, Delt_vth = 0.0, TempRatio = 0.0, Vth_NarrowW = 0.0, Lpe_Vb = 0.0, Vth = 0.0;
    double n = 0.0, n0 = 0.0, Vgsteff = 0.0, Vgs_eff = 0.0, niter = 0.0, toxpf = 0.0, toxpi = 0.0, Tcen = 0.0, toxe = 0.0, epsrox = 0.0, vddeot = 0.0;
    double vtfbphi2eot = 0.0, phieot = 0.0, TempRatioeot = 0.0, Vtm0eot = 0.0, Vtmeot = 0.0, vbieot = 0.0;
    std::shared_ptr<bsim4SizeDependParam> pParam = nullptr;
};

struct bsim4SizeDependParam
{
    double Width;
    double Length;
    double NFinger;

    double BSIM4cdsc;
    double BSIM4cdscb;
    double BSIM4cdscd;
    double BSIM4cit;
    double BSIM4nfactor;
    double BSIM4xj;
    double BSIM4vsat;
    double BSIM4at;
    double BSIM4a0;
    double BSIM4ags;
    double BSIM4a1;
    double BSIM4a2;
    double BSIM4keta;
    double BSIM4nsub;
    double BSIM4ndep;
    double BSIM4nsd;
    double BSIM4phin;
    double BSIM4ngate;
    double BSIM4gamma1;
    double BSIM4gamma2;
    double BSIM4vbx;
    double BSIM4vbi;
    double BSIM4vbm;
    double BSIM4xt;
    double BSIM4phi;
    double BSIM4litl;
    double BSIM4k1;
    double BSIM4kt1;
    double BSIM4kt1l;
    double BSIM4kt2;
    double BSIM4k2;
    double BSIM4k3;
    double BSIM4k3b;
    double BSIM4w0;
    double BSIM4dvtp0;
    double BSIM4dvtp1;
    double BSIM4dvtp2;  /* New DIBL/Rout */
    double BSIM4dvtp3;
    double BSIM4dvtp4;
    double BSIM4dvtp5;
    double BSIM4lpe0;
    double BSIM4lpeb;
    double BSIM4dvt0;
    double BSIM4dvt1;
    double BSIM4dvt2;
    double BSIM4dvt0w;
    double BSIM4dvt1w;
    double BSIM4dvt2w;
    double BSIM4drout;
    double BSIM4dsub;
    double BSIM4vth0;
    double BSIM4ua;
    double BSIM4ua1;
    double BSIM4ub;
    double BSIM4ub1;
    double BSIM4uc;
    double BSIM4uc1;
    double BSIM4ud;
    double BSIM4ud1;
    double BSIM4up;
    double BSIM4lp;
    double BSIM4u0;
    double BSIM4eu;
    double BSIM4ucs;
    double BSIM4ute;
    double BSIM4ucste;
    double BSIM4voff;
    double BSIM4tvoff;
    double BSIM4tnfactor;   /* v4.7 Temp dep of leakage current */
    double BSIM4teta0;      /* v4.7 temp dep of leakage current */
    double BSIM4tvoffcv;    /* v4.7 temp dep of leakage current */
    double BSIM4minv;
    double BSIM4minvcv;
    double BSIM4vfb;
    double BSIM4delta;
    double BSIM4rdsw;
    double BSIM4rds0;
    double BSIM4rs0;
    double BSIM4rd0;
    double BSIM4rsw;
    double BSIM4rdw;
    double BSIM4prwg;
    double BSIM4prwb;
    double BSIM4prt;
    double BSIM4eta0;
    double BSIM4etab;
    double BSIM4pclm;
    double BSIM4pdibl1;
    double BSIM4pdibl2;
    double BSIM4pdiblb;
    double BSIM4fprout;
    double BSIM4pdits;
    double BSIM4pditsd;
    double BSIM4pscbe1;
    double BSIM4pscbe2;
    double BSIM4pvag;
    double BSIM4wr;
    double BSIM4dwg;
    double BSIM4dwb;
    double BSIM4b0;
    double BSIM4b1;
    double BSIM4alpha0;
    double BSIM4alpha1;
    double BSIM4beta0;
    double BSIM4agidl;
    double BSIM4bgidl;
    double BSIM4cgidl;
    double BSIM4egidl;
    double BSIM4fgidl; /* v4.7 New GIDL/GISL */
    double BSIM4kgidl; /* v4.7 New GIDL/GISL */
    double BSIM4rgidl; /* v4.7 New GIDL/GISL */
    double BSIM4agisl;
    double BSIM4bgisl;
    double BSIM4cgisl;
    double BSIM4egisl;
    double BSIM4fgisl; /* v4.7 New GIDL/GISL */
    double BSIM4kgisl; /* v4.7 New GIDL/GISL */
    double BSIM4rgisl; /* v4.7 New GIDL/GISL */
    double BSIM4aigc;
    double BSIM4bigc;
    double BSIM4cigc;
    double BSIM4aigsd;
    double BSIM4bigsd;
    double BSIM4cigsd;
    double BSIM4aigs;
    double BSIM4bigs;
    double BSIM4cigs;
    double BSIM4aigd;
    double BSIM4bigd;
    double BSIM4cigd;
    double BSIM4aigbacc;
    double BSIM4bigbacc;
    double BSIM4cigbacc;
    double BSIM4aigbinv;
    double BSIM4bigbinv;
    double BSIM4cigbinv;
    double BSIM4nigc;
    double BSIM4nigbacc;
    double BSIM4nigbinv;
    double BSIM4ntox;
    double BSIM4eigbinv;
    double BSIM4pigcd;
    double BSIM4poxedge;
    double BSIM4xrcrg1;
    double BSIM4xrcrg2;
    double BSIM4lambda; /* overshoot */
    double BSIM4vtl; /* thermal velocity limit */
    double BSIM4xn; /* back scattering parameter */
    double BSIM4lc; /* back scattering parameter */
    double BSIM4tfactor;  /* ballistic transportation factor  */
    double BSIM4vfbsdoff;  /* S/D flatband offset voltage  */
    double BSIM4tvfbsdoff;

/* added for stress effect */
    double BSIM4ku0;
    double BSIM4kvth0;
    double BSIM4ku0temp;
    double BSIM4rho_ref;
    double BSIM4inv_od_ref;
/* added for well proximity effect */
    double BSIM4kvth0we;
    double BSIM4k2we;
    double BSIM4ku0we;

    /* CV model */
    double BSIM4cgsl;
    double BSIM4cgdl;
    double BSIM4ckappas;
    double BSIM4ckappad;
    double BSIM4cf;
    double BSIM4clc;
    double BSIM4cle;
    double BSIM4vfbcv;
    double BSIM4noff;
    double BSIM4voffcv;
    double BSIM4acde;
    double BSIM4moin;

/* Pre-calculated constants */

    double BSIM4dw;
    double BSIM4dl;
    double BSIM4leff;
    double BSIM4weff;

    double BSIM4dwc;
    double BSIM4dlc;
    double BSIM4dwj;
    double BSIM4leffCV;
    double BSIM4weffCV;
    double BSIM4weffCJ;
    double BSIM4abulkCVfactor;
    double BSIM4cgso;
    double BSIM4cgdo;
    double BSIM4cgbo;

    double BSIM4u0temp;
    double BSIM4vsattemp;
    double BSIM4sqrtPhi;
    double BSIM4phis3;
    double BSIM4Xdep0;
    double BSIM4sqrtXdep0;
    double BSIM4theta0vb0;
    double BSIM4thetaRout;
    double BSIM4mstar;
    double BSIM4VgsteffVth;
    double BSIM4mstarcv;
    double BSIM4voffcbn;
    double BSIM4voffcbncv;
    double BSIM4rdswmin;
    double BSIM4rdwmin;
    double BSIM4rswmin;
    double BSIM4vfbsd;

    double BSIM4cof1;
    double BSIM4cof2;
    double BSIM4cof3;
    double BSIM4cof4;
    double BSIM4cdep0;
    double BSIM4ToxRatio;
    double BSIM4Aechvb;
    double BSIM4Bechvb;
    double BSIM4ToxRatioEdge;
    double BSIM4AechvbEdgeS;
    double BSIM4AechvbEdgeD;
    double BSIM4BechvbEdge;
    double BSIM4ldeb;
    double BSIM4k1ox;
    double BSIM4k2ox;
    double BSIM4vfbzbfactor;
    double BSIM4dvtp2factor; /* v4.7 */
};

struct BSIM4model
{
    // int BSIM4modType;
    // struct sBSIM4model *BSIM4nextModel;
    bsim4v82temp BSIM4temp;
    std::string BSIM4modName;
    int BSIM4type;

    int    BSIM4mobMod;             // BSIM4 provides three different models of the effective mobility.
    int    BSIM4cvchargeMod;
    int    BSIM4capMod;
    int    BSIM4dioMod;
    int    BSIM4trnqsMod;
    int    BSIM4acnqsMod;
    int    BSIM4fnoiMod;
    int    BSIM4tnoiMod;
    int    BSIM4rdsMod;
    int    BSIM4rbodyMod;
    int    BSIM4rgateMod;
    int    BSIM4perMod;
    int    BSIM4geoMod;
    int    BSIM4mtrlMod;
    int    BSIM4mtrlCompatMod; /* v4.7 */
    int    BSIM4gidlMod; /* v4.7 New GIDL/GISL */
    int    BSIM4igcMod;
    int    BSIM4igbMod;
    int    BSIM4tempMod;
    int    BSIM4binUnit;
    int    BSIM4paramChk;
    std::string BSIM4version;
    double BSIM4eot;
    double BSIM4vddeot;
    double BSIM4tempeot;
    double BSIM4leffeot;
    double BSIM4weffeot;
    double BSIM4ados;
    double BSIM4bdos;
    double BSIM4toxe;
    double BSIM4toxp;
    double BSIM4toxm;
    double BSIM4dtox;
    double BSIM4epsrox;
    double BSIM4cdsc;
    double BSIM4cdscb;
    double BSIM4cdscd;
    double BSIM4cit;
    double BSIM4nfactor;
    double BSIM4xj;
    double BSIM4vsat;
    double BSIM4at;
    double BSIM4a0;
    double BSIM4ags;
    double BSIM4a1;
    double BSIM4a2;
    double BSIM4keta;
    double BSIM4nsub;
    double BSIM4phig;
    double BSIM4epsrgate;
    double BSIM4easub;
    double BSIM4epsrsub;
    double BSIM4ni0sub;
    double BSIM4bg0sub;
    double BSIM4tbgasub;
    double BSIM4tbgbsub;
    double BSIM4ndep;
    double BSIM4nsd;
    double BSIM4phin;
    double BSIM4ngate;
    double BSIM4gamma1;
    double BSIM4gamma2;
    double BSIM4vbx;
    double BSIM4vbm;
    double BSIM4xt;
    double BSIM4k1;
    double BSIM4kt1;
    double BSIM4kt1l;
    double BSIM4kt2;
    double BSIM4k2;
    double BSIM4k3;
    double BSIM4k3b;
    double BSIM4w0;
    double BSIM4dvtp0;
    double BSIM4dvtp1;
    double BSIM4dvtp2;  /* New DIBL/Rout */
    double BSIM4dvtp3;
    double BSIM4dvtp4;
    double BSIM4dvtp5;
    double BSIM4lpe0;
    double BSIM4lpeb;
    double BSIM4dvt0;
    double BSIM4dvt1;
    double BSIM4dvt2;
    double BSIM4dvt0w;
    double BSIM4dvt1w;
    double BSIM4dvt2w;
    double BSIM4drout;
    double BSIM4dsub;
    double BSIM4vth0;
    double BSIM4eu;
    double BSIM4ucs;
    double BSIM4ua;
    double BSIM4ua1;
    double BSIM4ub;
    double BSIM4ub1;
    double BSIM4uc;
    double BSIM4uc1;
    double BSIM4ud;
    double BSIM4ud1;
    double BSIM4up;
    double BSIM4lp;
    double BSIM4u0;
    double BSIM4ute;
    double BSIM4ucste;
    double BSIM4voff;
    double BSIM4tvoff;
    double BSIM4tnfactor;   /* v4.7 Temp dep of leakage current */
    double BSIM4teta0;      /* v4.7 temp dep of leakage current */

    double BSIM4tvoffcv;    /* v4.7 temp dep of leakage current */
    double BSIM4minv;
    double BSIM4minvcv;
    double BSIM4voffl;
    double BSIM4voffcvl;
    double BSIM4delta;
    double BSIM4rdsw;
    double BSIM4rdswmin;
    double BSIM4rdwmin;
    double BSIM4rswmin;
    double BSIM4rsw;
    double BSIM4rdw;
    double BSIM4prwg;
    double BSIM4prwb;
    double BSIM4prt;
    double BSIM4eta0;
    double BSIM4etab;
    double BSIM4pclm;
    double BSIM4pdibl1;
    double BSIM4pdibl2;
    double BSIM4pdiblb;
    double BSIM4fprout;
    double BSIM4pdits;
    double BSIM4pditsd;
    double BSIM4pditsl;
    double BSIM4pscbe1;
    double BSIM4pscbe2;
    double BSIM4pvag;
    double BSIM4wr;
    double BSIM4dwg;
    double BSIM4dwb;
    double BSIM4b0;
    double BSIM4b1;
    double BSIM4alpha0;
    double BSIM4alpha1;
    double BSIM4beta0;
    double BSIM4agidl;
    double BSIM4bgidl;
    double BSIM4cgidl;
    double BSIM4egidl;
    double BSIM4fgidl; /* v4.7 New GIDL/GISL */
    double BSIM4kgidl; /* v4.7 New GIDL/GISL */
    double BSIM4rgidl; /* v4.7 New GIDL/GISL */
    double BSIM4agisl;
    double BSIM4bgisl;
    double BSIM4cgisl;
    double BSIM4egisl;
    double BSIM4fgisl; /* v4.7 New GIDL/GISL */
    double BSIM4kgisl; /* v4.7 New GIDL/GISL */
    double BSIM4rgisl; /* v4.7 New GIDL/GISL */
    double BSIM4aigc;
    double BSIM4bigc;
    double BSIM4cigc;
    double BSIM4aigsd;
    double BSIM4bigsd;
    double BSIM4cigsd;
    double BSIM4aigs;
    double BSIM4bigs;
    double BSIM4cigs;
    double BSIM4aigd;
    double BSIM4bigd;
    double BSIM4cigd;
    double BSIM4aigbacc;
    double BSIM4bigbacc;
    double BSIM4cigbacc;
    double BSIM4aigbinv;
    double BSIM4bigbinv;
    double BSIM4cigbinv;
    double BSIM4nigc;
    double BSIM4nigbacc;
    double BSIM4nigbinv;
    double BSIM4ntox;
    double BSIM4eigbinv;
    double BSIM4pigcd;
    double BSIM4poxedge;
    double BSIM4toxref;
    double BSIM4ijthdfwd;
    double BSIM4ijthsfwd;
    double BSIM4ijthdrev;
    double BSIM4ijthsrev;
    double BSIM4xjbvd;
    double BSIM4xjbvs;
    double BSIM4bvd;
    double BSIM4bvs;

    double BSIM4jtss;
    double BSIM4jtsd;
    double BSIM4jtssws;
    double BSIM4jtsswd;
    double BSIM4jtsswgs;
    double BSIM4jtsswgd;
    double BSIM4jtweff;
    double BSIM4njts;
    double BSIM4njtssw;
    double BSIM4njtsswg;
    double BSIM4njtsd;
    double BSIM4njtsswd;
    double BSIM4njtsswgd;
    double BSIM4xtss;
    double BSIM4xtsd;
    double BSIM4xtssws;
    double BSIM4xtsswd;
    double BSIM4xtsswgs;
    double BSIM4xtsswgd;
    double BSIM4tnjts;
    double BSIM4tnjtssw;
    double BSIM4tnjtsswg;
    double BSIM4tnjtsd;
    double BSIM4tnjtsswd;
    double BSIM4tnjtsswgd;
    double BSIM4vtss;
    double BSIM4vtsd;
    double BSIM4vtssws;
    double BSIM4vtsswd;
    double BSIM4vtsswgs;
    double BSIM4vtsswgd;

    double BSIM4xrcrg1;
    double BSIM4xrcrg2;
    double BSIM4lambda;
    double BSIM4vtl;
    double BSIM4lc;
    double BSIM4xn;
    double BSIM4vfbsdoff;  /* S/D flatband offset voltage  */
    double BSIM4lintnoi;  /* lint offset for noise calculation  */
    double BSIM4tvfbsdoff;

    double BSIM4vfb;
    double BSIM4gbmin;
    double BSIM4rbdb;
    double BSIM4rbsb;
    double BSIM4rbpb;
    double BSIM4rbps;
    double BSIM4rbpd;

    double BSIM4rbps0;
    double BSIM4rbpsl;
    double BSIM4rbpsw;
    double BSIM4rbpsnf;

    double BSIM4rbpd0;
    double BSIM4rbpdl;
    double BSIM4rbpdw;
    double BSIM4rbpdnf;

    double BSIM4rbpbx0;
    double BSIM4rbpbxl;
    double BSIM4rbpbxw;
    double BSIM4rbpbxnf;
    double BSIM4rbpby0;
    double BSIM4rbpbyl;
    double BSIM4rbpbyw;
    double BSIM4rbpbynf;

    double BSIM4rbsbx0;
    double BSIM4rbsby0;
    double BSIM4rbdbx0;
    double BSIM4rbdby0;

    double BSIM4rbsdbxl;
    double BSIM4rbsdbxw;
    double BSIM4rbsdbxnf;
    double BSIM4rbsdbyl;
    double BSIM4rbsdbyw;
    double BSIM4rbsdbynf;

    double BSIM4tnoia;
    double BSIM4tnoib;
    double BSIM4tnoic;
    double BSIM4rnoia;
    double BSIM4rnoib;
    double BSIM4rnoic;
    double BSIM4ntnoi;

    /* CV model and Parasitics */
    double BSIM4cgsl;
    double BSIM4cgdl;
    double BSIM4ckappas;
    double BSIM4ckappad;
    double BSIM4cf;
    double BSIM4vfbcv;
    double BSIM4clc;
    double BSIM4cle;
    double BSIM4dwc;
    double BSIM4dlc;
    double BSIM4xw;
    double BSIM4xl;
    double BSIM4dlcig;
    double BSIM4dlcigd;
    double BSIM4dwj;
    double BSIM4noff;
    double BSIM4voffcv;
    double BSIM4acde;
    double BSIM4moin;
    double BSIM4tcj;
    double BSIM4tcjsw;
    double BSIM4tcjswg;
    double BSIM4tpb;
    double BSIM4tpbsw;
    double BSIM4tpbswg;
    double BSIM4dmcg;
    double BSIM4dmci;
    double BSIM4dmdg;
    double BSIM4dmcgt;
    double BSIM4xgw;
    double BSIM4xgl;
    double BSIM4rshg;
    double BSIM4ngcon;

    /* Length Dependence */
    double BSIM4lcdsc;
    double BSIM4lcdscb;
    double BSIM4lcdscd;
    double BSIM4lcit;
    double BSIM4lnfactor;
    double BSIM4lxj;
    double BSIM4lvsat;
    double BSIM4lat;
    double BSIM4la0;
    double BSIM4lags;
    double BSIM4la1;
    double BSIM4la2;
    double BSIM4lketa;
    double BSIM4lnsub;
    double BSIM4lndep;
    double BSIM4lnsd;
    double BSIM4lphin;
    double BSIM4lngate;
    double BSIM4lgamma1;
    double BSIM4lgamma2;
    double BSIM4lvbx;
    double BSIM4lvbm;
    double BSIM4lxt;
    double BSIM4lk1;
    double BSIM4lkt1;
    double BSIM4lkt1l;
    double BSIM4lkt2;
    double BSIM4lk2;
    double BSIM4lk3;
    double BSIM4lk3b;
    double BSIM4lw0;
    double BSIM4ldvtp0;
    double BSIM4ldvtp1;
    double BSIM4ldvtp2; /* New DIBL/Rout */
    double BSIM4ldvtp3;
    double BSIM4ldvtp4;
    double BSIM4ldvtp5;
    double BSIM4llpe0;
    double BSIM4llpeb;
    double BSIM4ldvt0;
    double BSIM4ldvt1;
    double BSIM4ldvt2;
    double BSIM4ldvt0w;
    double BSIM4ldvt1w;
    double BSIM4ldvt2w;
    double BSIM4ldrout;
    double BSIM4ldsub;
    double BSIM4lvth0;
    double BSIM4lua;
    double BSIM4lua1;
    double BSIM4lub;
    double BSIM4lub1;
    double BSIM4luc;
    double BSIM4luc1;
    double BSIM4lud;
    double BSIM4lud1;
    double BSIM4lup;
    double BSIM4llp;
    double BSIM4lu0;
    double BSIM4leu;
    double BSIM4lucs;
    double BSIM4lute;
    double BSIM4lucste;
    double BSIM4lvoff;
    double BSIM4ltvoff;
    double BSIM4ltnfactor;  /* v4.7 Temp dep of leakage current */
    double BSIM4lteta0;     /* v4.7 temp dep of leakage current */
    double BSIM4ltvoffcv;   /* v4.7 temp dep of leakage current */
    double BSIM4lminv;
    double BSIM4lminvcv;
    double BSIM4ldelta;
    double BSIM4lrdsw;
    double BSIM4lrsw;
    double BSIM4lrdw;
    double BSIM4lprwg;
    double BSIM4lprwb;
    double BSIM4lprt;
    double BSIM4leta0;
    double BSIM4letab;
    double BSIM4lpclm;
    double BSIM4lpdibl1;
    double BSIM4lpdibl2;
    double BSIM4lpdiblb;
    double BSIM4lfprout;
    double BSIM4lpdits;
    double BSIM4lpditsd;
    double BSIM4lpscbe1;
    double BSIM4lpscbe2;
    double BSIM4lpvag;
    double BSIM4lwr;
    double BSIM4ldwg;
    double BSIM4ldwb;
    double BSIM4lb0;
    double BSIM4lb1;
    double BSIM4lalpha0;
    double BSIM4lalpha1;
    double BSIM4lbeta0;
    double BSIM4lvfb;
    double BSIM4lagidl;
    double BSIM4lbgidl;
    double BSIM4lcgidl;
    double BSIM4legidl;
    double BSIM4lfgidl; /* v4.7 New GIDL/GISL */
    double BSIM4lkgidl; /* v4.7 New GIDL/GISL */
    double BSIM4lrgidl; /* v4.7 New GIDL/GISL */
    double BSIM4lagisl;
    double BSIM4lbgisl;
    double BSIM4lcgisl;
    double BSIM4legisl;
    double BSIM4lfgisl; /* v4.7 New GIDL/GISL */
    double BSIM4lkgisl; /* v4.7 New GIDL/GISL */
    double BSIM4lrgisl; /* v4.7 New GIDL/GISL */
    double BSIM4laigc;
    double BSIM4lbigc;
    double BSIM4lcigc;
    double BSIM4laigsd;
    double BSIM4lbigsd;
    double BSIM4lcigsd;
    double BSIM4laigs;
    double BSIM4lbigs;
    double BSIM4lcigs;
    double BSIM4laigd;
    double BSIM4lbigd;
    double BSIM4lcigd;
    double BSIM4laigbacc;
    double BSIM4lbigbacc;
    double BSIM4lcigbacc;
    double BSIM4laigbinv;
    double BSIM4lbigbinv;
    double BSIM4lcigbinv;
    double BSIM4lnigc;
    double BSIM4lnigbacc;
    double BSIM4lnigbinv;
    double BSIM4lntox;
    double BSIM4leigbinv;
    double BSIM4lpigcd;
    double BSIM4lpoxedge;
    double BSIM4lxrcrg1;
    double BSIM4lxrcrg2;
    double BSIM4llambda;
    double BSIM4lvtl;
    double BSIM4lxn;
    double BSIM4lvfbsdoff;
    double BSIM4ltvfbsdoff;

    /* CV model */
    double BSIM4lcgsl;
    double BSIM4lcgdl;
    double BSIM4lckappas;
    double BSIM4lckappad;
    double BSIM4lcf;
    double BSIM4lclc;
    double BSIM4lcle;
    double BSIM4lvfbcv;
    double BSIM4lnoff;
    double BSIM4lvoffcv;
    double BSIM4lacde;
    double BSIM4lmoin;

    /* Width Dependence */
    double BSIM4wcdsc;
    double BSIM4wcdscb;
    double BSIM4wcdscd;
    double BSIM4wcit;
    double BSIM4wnfactor;
    double BSIM4wxj;
    double BSIM4wvsat;
    double BSIM4wat;
    double BSIM4wa0;
    double BSIM4wags;
    double BSIM4wa1;
    double BSIM4wa2;
    double BSIM4wketa;
    double BSIM4wnsub;
    double BSIM4wndep;
    double BSIM4wnsd;
    double BSIM4wphin;
    double BSIM4wngate;
    double BSIM4wgamma1;
    double BSIM4wgamma2;
    double BSIM4wvbx;
    double BSIM4wvbm;
    double BSIM4wxt;
    double BSIM4wk1;
    double BSIM4wkt1;
    double BSIM4wkt1l;
    double BSIM4wkt2;
    double BSIM4wk2;
    double BSIM4wk3;
    double BSIM4wk3b;
    double BSIM4ww0;
    double BSIM4wdvtp0;
    double BSIM4wdvtp1;
    double BSIM4wdvtp2; /* New DIBL/Rout */
    double BSIM4wdvtp3;
    double BSIM4wdvtp4;
    double BSIM4wdvtp5;
    double BSIM4wlpe0;
    double BSIM4wlpeb;
    double BSIM4wdvt0;
    double BSIM4wdvt1;
    double BSIM4wdvt2;
    double BSIM4wdvt0w;
    double BSIM4wdvt1w;
    double BSIM4wdvt2w;
    double BSIM4wdrout;
    double BSIM4wdsub;
    double BSIM4wvth0;
    double BSIM4wua;
    double BSIM4wua1;
    double BSIM4wub;
    double BSIM4wub1;
    double BSIM4wuc;
    double BSIM4wuc1;
    double BSIM4wud;
    double BSIM4wud1;
    double BSIM4wup;
    double BSIM4wlp;
    double BSIM4wu0;
    double BSIM4weu;
    double BSIM4wucs;
    double BSIM4wute;
    double BSIM4wucste;
    double BSIM4wvoff;
    double BSIM4wtvoff;
    double BSIM4wtnfactor;  /* v4.7 Temp dep of leakage current */
    double BSIM4wteta0;     /* v4.7 temp dep of leakage current */
    double BSIM4wtvoffcv;   /* v4.7 temp dep of leakage current */
    double BSIM4wminv;
    double BSIM4wminvcv;
    double BSIM4wdelta;
    double BSIM4wrdsw;
    double BSIM4wrsw;
    double BSIM4wrdw;
    double BSIM4wprwg;
    double BSIM4wprwb;
    double BSIM4wprt;
    double BSIM4weta0;
    double BSIM4wetab;
    double BSIM4wpclm;
    double BSIM4wpdibl1;
    double BSIM4wpdibl2;
    double BSIM4wpdiblb;
    double BSIM4wfprout;
    double BSIM4wpdits;
    double BSIM4wpditsd;
    double BSIM4wpscbe1;
    double BSIM4wpscbe2;
    double BSIM4wpvag;
    double BSIM4wwr;
    double BSIM4wdwg;
    double BSIM4wdwb;
    double BSIM4wb0;
    double BSIM4wb1;
    double BSIM4walpha0;
    double BSIM4walpha1;
    double BSIM4wbeta0;
    double BSIM4wvfb;
    double BSIM4wagidl;
    double BSIM4wbgidl;
    double BSIM4wcgidl;
    double BSIM4wegidl;
    double BSIM4wfgidl; /* v4.7 New GIDL/GISL */
    double BSIM4wkgidl; /* v4.7 New GIDL/GISL */
    double BSIM4wrgidl; /* v4.7 New GIDL/GISL */
    double BSIM4wagisl;
    double BSIM4wbgisl;
    double BSIM4wcgisl;
    double BSIM4wegisl;
    double BSIM4wfgisl; /* v4.7 New GIDL/GISL */
    double BSIM4wkgisl; /* v4.7 New GIDL/GISL */
    double BSIM4wrgisl; /* v4.7 New GIDL/GISL */
    double BSIM4waigc;
    double BSIM4wbigc;
    double BSIM4wcigc;
    double BSIM4waigsd;
    double BSIM4wbigsd;
    double BSIM4wcigsd;
    double BSIM4waigs;
    double BSIM4wbigs;
    double BSIM4wcigs;
    double BSIM4waigd;
    double BSIM4wbigd;
    double BSIM4wcigd;
    double BSIM4waigbacc;
    double BSIM4wbigbacc;
    double BSIM4wcigbacc;
    double BSIM4waigbinv;
    double BSIM4wbigbinv;
    double BSIM4wcigbinv;
    double BSIM4wnigc;
    double BSIM4wnigbacc;
    double BSIM4wnigbinv;
    double BSIM4wntox;
    double BSIM4weigbinv;
    double BSIM4wpigcd;
    double BSIM4wpoxedge;
    double BSIM4wxrcrg1;
    double BSIM4wxrcrg2;
    double BSIM4wlambda;
    double BSIM4wvtl;
    double BSIM4wxn;
    double BSIM4wvfbsdoff;
    double BSIM4wtvfbsdoff;

    /* CV model */
    double BSIM4wcgsl;
    double BSIM4wcgdl;
    double BSIM4wckappas;
    double BSIM4wckappad;
    double BSIM4wcf;
    double BSIM4wclc;
    double BSIM4wcle;
    double BSIM4wvfbcv;
    double BSIM4wnoff;
    double BSIM4wvoffcv;
    double BSIM4wacde;
    double BSIM4wmoin;

    /* Cross-term Dependence */
    double BSIM4pcdsc;
    double BSIM4pcdscb;
    double BSIM4pcdscd;
    double BSIM4pcit;
    double BSIM4pnfactor;
    double BSIM4pxj;
    double BSIM4pvsat;
    double BSIM4pat;
    double BSIM4pa0;
    double BSIM4pags;
    double BSIM4pa1;
    double BSIM4pa2;
    double BSIM4pketa;
    double BSIM4pnsub;
    double BSIM4pndep;
    double BSIM4pnsd;
    double BSIM4pphin;
    double BSIM4pngate;
    double BSIM4pgamma1;
    double BSIM4pgamma2;
    double BSIM4pvbx;
    double BSIM4pvbm;
    double BSIM4pxt;
    double BSIM4pk1;
    double BSIM4pkt1;
    double BSIM4pkt1l;
    double BSIM4pkt2;
    double BSIM4pk2;
    double BSIM4pk3;
    double BSIM4pk3b;
    double BSIM4pw0;
    double BSIM4pdvtp0;
    double BSIM4pdvtp1;
    double BSIM4pdvtp2; /* New DIBL/Rout */
    double BSIM4pdvtp3;
    double BSIM4pdvtp4;
    double BSIM4pdvtp5;
    double BSIM4plpe0;
    double BSIM4plpeb;
    double BSIM4pdvt0;
    double BSIM4pdvt1;
    double BSIM4pdvt2;
    double BSIM4pdvt0w;
    double BSIM4pdvt1w;
    double BSIM4pdvt2w;
    double BSIM4pdrout;
    double BSIM4pdsub;
    double BSIM4pvth0;
    double BSIM4pua;
    double BSIM4pua1;
    double BSIM4pub;
    double BSIM4pub1;
    double BSIM4puc;
    double BSIM4puc1;
    double BSIM4pud;
    double BSIM4pud1;
    double BSIM4pup;
    double BSIM4plp;
    double BSIM4pu0;
    double BSIM4peu;
    double BSIM4pucs;
    double BSIM4pute;
    double BSIM4pucste;
    double BSIM4pvoff;
    double BSIM4ptvoff;
    double BSIM4ptnfactor;  /* v4.7 Temp dep of leakage current */
    double BSIM4pteta0;     /* v4.7 temp dep of leakage current */
    double BSIM4ptvoffcv;   /* v4.7 temp dep of leakage current */
    double BSIM4pminv;
    double BSIM4pminvcv;
    double BSIM4pdelta;
    double BSIM4prdsw;
    double BSIM4prsw;
    double BSIM4prdw;
    double BSIM4pprwg;
    double BSIM4pprwb;
    double BSIM4pprt;
    double BSIM4peta0;
    double BSIM4petab;
    double BSIM4ppclm;
    double BSIM4ppdibl1;
    double BSIM4ppdibl2;
    double BSIM4ppdiblb;
    double BSIM4pfprout;
    double BSIM4ppdits;
    double BSIM4ppditsd;
    double BSIM4ppscbe1;
    double BSIM4ppscbe2;
    double BSIM4ppvag;
    double BSIM4pwr;
    double BSIM4pdwg;
    double BSIM4pdwb;
    double BSIM4pb0;
    double BSIM4pb1;
    double BSIM4palpha0;
    double BSIM4palpha1;
    double BSIM4pbeta0;
    double BSIM4pvfb;
    double BSIM4pagidl;
    double BSIM4pbgidl;
    double BSIM4pcgidl;
    double BSIM4pegidl;
    double BSIM4pfgidl; /* v4.7 New GIDL/GISL */
    double BSIM4pkgidl; /* v4.7 New GIDL/GISL */
    double BSIM4prgidl; /* v4.7 New GIDL/GISL */
    double BSIM4pagisl;
    double BSIM4pbgisl;
    double BSIM4pcgisl;
    double BSIM4pegisl;
    double BSIM4pfgisl; /* v4.7 New GIDL/GISL */
    double BSIM4pkgisl; /* v4.7 New GIDL/GISL */
    double BSIM4prgisl; /* v4.7 New GIDL/GISL */
    double BSIM4paigc;
    double BSIM4pbigc;
    double BSIM4pcigc;
    double BSIM4paigsd;
    double BSIM4pbigsd;
    double BSIM4pcigsd;
    double BSIM4paigs;
    double BSIM4pbigs;
    double BSIM4pcigs;
    double BSIM4paigd;
    double BSIM4pbigd;
    double BSIM4pcigd;
    double BSIM4paigbacc;
    double BSIM4pbigbacc;
    double BSIM4pcigbacc;
    double BSIM4paigbinv;
    double BSIM4pbigbinv;
    double BSIM4pcigbinv;
    double BSIM4pnigc;
    double BSIM4pnigbacc;
    double BSIM4pnigbinv;
    double BSIM4pntox;
    double BSIM4peigbinv;
    double BSIM4ppigcd;
    double BSIM4ppoxedge;
    double BSIM4pxrcrg1;
    double BSIM4pxrcrg2;
    double BSIM4plambda;
    double BSIM4pvtl;
    double BSIM4pxn;
    double BSIM4pvfbsdoff;
    double BSIM4ptvfbsdoff;

    /* CV model */
    double BSIM4pcgsl;
    double BSIM4pcgdl;
    double BSIM4pckappas;
    double BSIM4pckappad;
    double BSIM4pcf;
    double BSIM4pclc;
    double BSIM4pcle;
    double BSIM4pvfbcv;
    double BSIM4pnoff;
    double BSIM4pvoffcv;
    double BSIM4pacde;
    double BSIM4pmoin;

    double BSIM4tnom;
    double BSIM4cgso;
    double BSIM4cgdo;
    double BSIM4cgbo;
    double BSIM4xpart;
    double BSIM4cFringOut;
    double BSIM4cFringMax;

    double BSIM4sheetResistance;
    double BSIM4SjctSatCurDensity;
    double BSIM4DjctSatCurDensity;
    double BSIM4SjctSidewallSatCurDensity;
    double BSIM4DjctSidewallSatCurDensity;
    double BSIM4SjctGateSidewallSatCurDensity;
    double BSIM4DjctGateSidewallSatCurDensity;
    double BSIM4SbulkJctPotential;
    double BSIM4DbulkJctPotential;
    double BSIM4SbulkJctBotGradingCoeff;
    double BSIM4DbulkJctBotGradingCoeff;
    double BSIM4SbulkJctSideGradingCoeff;
    double BSIM4DbulkJctSideGradingCoeff;
    double BSIM4SbulkJctGateSideGradingCoeff;
    double BSIM4DbulkJctGateSideGradingCoeff;
    double BSIM4SsidewallJctPotential;
    double BSIM4DsidewallJctPotential;
    double BSIM4SGatesidewallJctPotential;
    double BSIM4DGatesidewallJctPotential;
    double BSIM4SunitAreaJctCap;
    double BSIM4DunitAreaJctCap;
    double BSIM4SunitLengthSidewallJctCap;
    double BSIM4DunitLengthSidewallJctCap;
    double BSIM4SunitLengthGateSidewallJctCap;
    double BSIM4DunitLengthGateSidewallJctCap;
    double BSIM4SjctEmissionCoeff;
    double BSIM4DjctEmissionCoeff;
    double BSIM4SjctTempExponent;
    double BSIM4DjctTempExponent;
    double BSIM4njtsstemp;
    double BSIM4njtsswstemp;
    double BSIM4njtsswgstemp;
    double BSIM4njtsdtemp;
    double BSIM4njtsswdtemp;
    double BSIM4njtsswgdtemp;

    double BSIM4Lint;
    double BSIM4Ll;
    double BSIM4Llc;
    double BSIM4Lln;
    double BSIM4Lw;
    double BSIM4Lwc;
    double BSIM4Lwn;
    double BSIM4Lwl;
    double BSIM4Lwlc;
    double BSIM4Lmin;
    double BSIM4Lmax;

    double BSIM4Wint;
    double BSIM4Wl;
    double BSIM4Wlc;
    double BSIM4Wln;
    double BSIM4Ww;
    double BSIM4Wwc;
    double BSIM4Wwn;
    double BSIM4Wwl;
    double BSIM4Wwlc;
    double BSIM4Wmin;
    double BSIM4Wmax;

    /* added for stress effect */
    double BSIM4saref;
    double BSIM4sbref;
    double BSIM4wlod;
    double BSIM4ku0;
    double BSIM4kvsat;
    double BSIM4kvth0;
    double BSIM4tku0;
    double BSIM4llodku0;
    double BSIM4wlodku0;
    double BSIM4llodvth;
    double BSIM4wlodvth;
    double BSIM4lku0;
    double BSIM4wku0;
    double BSIM4pku0;
    double BSIM4lkvth0;
    double BSIM4wkvth0;
    double BSIM4pkvth0;
    double BSIM4stk2;
    double BSIM4lodk2;
    double BSIM4steta0;
    double BSIM4lodeta0;

    double BSIM4web;
    double BSIM4wec;
    double BSIM4kvth0we;
    double BSIM4k2we;
    double BSIM4ku0we;
    double BSIM4scref;
    double BSIM4wpemod;
    double BSIM4lkvth0we;
    double BSIM4lk2we;
    double BSIM4lku0we;
    double BSIM4wkvth0we;
    double BSIM4wk2we;
    double BSIM4wku0we;
    double BSIM4pkvth0we;
    double BSIM4pk2we;
    double BSIM4pku0we;

/* Pre-calculated constants
 * move to size-dependent param */

    double BSIM4Eg0;
    double BSIM4vtm;
    double BSIM4vtm0;
    double BSIM4coxe;
    double BSIM4coxp;
    double BSIM4cof1;
    double BSIM4cof2;
    double BSIM4cof3;
    double BSIM4cof4;
    double BSIM4vcrit;
    double BSIM4factor1;
    double BSIM4PhiBS;
    double BSIM4PhiBSWS;
    double BSIM4PhiBSWGS;
    double BSIM4SjctTempSatCurDensity;
    double BSIM4SjctSidewallTempSatCurDensity;
    double BSIM4SjctGateSidewallTempSatCurDensity;
    double BSIM4PhiBD;
    double BSIM4PhiBSWD;
    double BSIM4PhiBSWGD;
    double BSIM4DjctTempSatCurDensity;
    double BSIM4DjctSidewallTempSatCurDensity;
    double BSIM4DjctGateSidewallTempSatCurDensity;
    double BSIM4SunitAreaTempJctCap;
    double BSIM4DunitAreaTempJctCap;
    double BSIM4SunitLengthSidewallTempJctCap;
    double BSIM4DunitLengthSidewallTempJctCap;
    double BSIM4SunitLengthGateSidewallTempJctCap;
    double BSIM4DunitLengthGateSidewallTempJctCap;

    double BSIM4oxideTrapDensityA;
    double BSIM4oxideTrapDensityB;
    double BSIM4oxideTrapDensityC;
    double BSIM4em;
    double BSIM4ef;
    double BSIM4af;
    double BSIM4kf;
    double BSIM4gidlclamp;
    double BSIM4idovvdsc;
    // struct bsim4SizeDependParam *pSizeDependParamKnot;
    std::vector<std::shared_ptr<bsim4SizeDependParam>> vSizeDependParamKnot;

    double DMCGeff, DMCIeff, DMDGeff, Rtot;

    /* Flags */
    bool  BSIM4mobModGiven = false;
    bool  BSIM4binUnitGiven = false;
    bool  BSIM4cvchargeModGiven = false;
    bool  BSIM4capModGiven = false;
    bool  BSIM4dioModGiven = false;
    bool  BSIM4rdsModGiven = false;
    bool  BSIM4rbodyModGiven = false;
    bool  BSIM4rgateModGiven = false;
    bool  BSIM4perModGiven = false;
    bool  BSIM4geoModGiven = false;
    bool  BSIM4paramChkGiven = false;
    bool  BSIM4trnqsModGiven = false;
    bool  BSIM4acnqsModGiven = false;
    bool  BSIM4fnoiModGiven = false;
    bool  BSIM4tnoiModGiven = false;
    bool  BSIM4mtrlModGiven = false;
    bool  BSIM4mtrlCompatModGiven = false;
    bool  BSIM4gidlModGiven = false;    /* v4.7 New GIDL/GISL */
    bool  BSIM4igcModGiven = false;
    bool  BSIM4igbModGiven = false;
    bool  BSIM4tempModGiven = false;
    bool  BSIM4typeGiven   = false;
    bool  BSIM4toxrefGiven   = false;
    bool  BSIM4eotGiven   = false;
    bool  BSIM4vddeotGiven   = false;
    bool  BSIM4tempeotGiven  = false;
    bool  BSIM4leffeotGiven  = false;
    bool  BSIM4weffeotGiven  = false;
    bool  BSIM4adosGiven   = false;
    bool  BSIM4bdosGiven   = false;
    bool  BSIM4toxeGiven   = false;
    bool  BSIM4toxpGiven   = false;
    bool  BSIM4toxmGiven   = false;
    bool  BSIM4dtoxGiven   = false;
    bool  BSIM4epsroxGiven   = false;
    bool  BSIM4versionGiven   = false;
    bool  BSIM4cdscGiven   = false;
    bool  BSIM4cdscbGiven   = false;
    bool  BSIM4cdscdGiven   = false;
    bool  BSIM4citGiven   = false;
    bool  BSIM4nfactorGiven   = false;
    bool  BSIM4xjGiven   = false;
    bool  BSIM4vsatGiven   = false;
    bool  BSIM4atGiven   = false;
    bool  BSIM4a0Given   = false;
    bool  BSIM4agsGiven   = false;
    bool  BSIM4a1Given   = false;
    bool  BSIM4a2Given   = false;
    bool  BSIM4ketaGiven   = false;
    bool  BSIM4nsubGiven   = false;
    bool  BSIM4phigGiven   = false;
    bool  BSIM4epsrgateGiven   = false;
    bool  BSIM4easubGiven   = false;
    bool  BSIM4epsrsubGiven   = false;
    bool  BSIM4ni0subGiven   = false;
    bool  BSIM4bg0subGiven   = false;
    bool  BSIM4tbgasubGiven   = false;
    bool  BSIM4tbgbsubGiven   = false;
    bool  BSIM4ndepGiven   = false;
    bool  BSIM4nsdGiven    = false;
    bool  BSIM4phinGiven   = false;
    bool  BSIM4ngateGiven   = false;
    bool  BSIM4gamma1Given   = false;
    bool  BSIM4gamma2Given   = false;
    bool  BSIM4vbxGiven   = false;
    bool  BSIM4vbmGiven   = false;
    bool  BSIM4xtGiven   = false;
    bool  BSIM4k1Given   = false;
    bool  BSIM4kt1Given   = false;
    bool  BSIM4kt1lGiven   = false;
    bool  BSIM4kt2Given   = false;
    bool  BSIM4k2Given   = false;
    bool  BSIM4k3Given   = false;
    bool  BSIM4k3bGiven   = false;
    bool  BSIM4w0Given   = false;
    bool  BSIM4dvtp0Given = false;
    bool  BSIM4dvtp1Given = false;
    bool  BSIM4dvtp2Given = false;   /* New DIBL/Rout */
    bool  BSIM4dvtp3Given = false;
    bool  BSIM4dvtp4Given = false;
    bool  BSIM4dvtp5Given = false;
    bool  BSIM4lpe0Given   = false;
    bool  BSIM4lpebGiven   = false;
    bool  BSIM4dvt0Given   = false;
    bool  BSIM4dvt1Given   = false;
    bool  BSIM4dvt2Given   = false;
    bool  BSIM4dvt0wGiven   = false;
    bool  BSIM4dvt1wGiven   = false;
    bool  BSIM4dvt2wGiven   = false;
    bool  BSIM4droutGiven   = false;
    bool  BSIM4dsubGiven   = false;
    bool  BSIM4vth0Given   = false;
    bool  BSIM4euGiven   = false;
    bool  BSIM4ucsGiven  = false;
    bool  BSIM4uaGiven   = false;
    bool  BSIM4ua1Given   = false;
    bool  BSIM4ubGiven   = false;
    bool  BSIM4ub1Given   = false;
    bool  BSIM4ucGiven   = false;
    bool  BSIM4uc1Given   = false;
    bool  BSIM4udGiven     = false;
    bool  BSIM4ud1Given     = false;
    bool  BSIM4upGiven     = false;
    bool  BSIM4lpGiven     = false;
    bool  BSIM4u0Given   = false;
    bool  BSIM4uteGiven   = false;
    bool  BSIM4ucsteGiven = false;
    bool  BSIM4voffGiven   = false;
    bool  BSIM4tvoffGiven   = false;
    bool  BSIM4tnfactorGiven  = false;   /* v4.7 Temp dep of leakage current */
    bool  BSIM4teta0Given   = false;     /* v4.7 temp dep of leakage current */
    bool  BSIM4tvoffcvGiven   = false;   /* v4.7 temp dep of leakage current */
    bool  BSIM4vofflGiven  = false;
    bool  BSIM4voffcvlGiven  = false;
    bool  BSIM4minvGiven   = false;
    bool  BSIM4minvcvGiven   = false;
    bool  BSIM4rdswGiven   = false;
    bool  BSIM4rdswminGiven = false;
    bool  BSIM4rdwminGiven = false;
    bool  BSIM4rswminGiven = false;
    bool  BSIM4rswGiven   = false;
    bool  BSIM4rdwGiven   = false;
    bool  BSIM4prwgGiven   = false;
    bool  BSIM4prwbGiven   = false;
    bool  BSIM4prtGiven   = false;
    bool  BSIM4eta0Given   = false;
    bool  BSIM4etabGiven   = false;
    bool  BSIM4pclmGiven   = false;
    bool  BSIM4pdibl1Given   = false;
    bool  BSIM4pdibl2Given   = false;
    bool  BSIM4pdiblbGiven   = false;
    bool  BSIM4fproutGiven   = false;
    bool  BSIM4pditsGiven    = false;
    bool  BSIM4pditsdGiven    = false;
    bool  BSIM4pditslGiven    = false;
    bool  BSIM4pscbe1Given   = false;
    bool  BSIM4pscbe2Given   = false;
    bool  BSIM4pvagGiven   = false;
    bool  BSIM4deltaGiven  = false;
    bool  BSIM4wrGiven   = false;
    bool  BSIM4dwgGiven   = false;
    bool  BSIM4dwbGiven   = false;
    bool  BSIM4b0Given   = false;
    bool  BSIM4b1Given   = false;
    bool  BSIM4alpha0Given   = false;
    bool  BSIM4alpha1Given   = false;
    bool  BSIM4beta0Given   = false;
    bool  BSIM4agidlGiven   = false;
    bool  BSIM4bgidlGiven   = false;
    bool  BSIM4cgidlGiven   = false;
    bool  BSIM4egidlGiven   = false;
    bool  BSIM4fgidlGiven   = false;    /* v4.7 New GIDL/GISL */
    bool  BSIM4kgidlGiven   = false;    /* v4.7 New GIDL/GISL */
    bool  BSIM4rgidlGiven   = false;    /* v4.7 New GIDL/GISL */
    bool  BSIM4agislGiven   = false;
    bool  BSIM4bgislGiven   = false;
    bool  BSIM4cgislGiven   = false;
    bool  BSIM4egislGiven   = false;
    bool  BSIM4fgislGiven   = false;    /* v4.7 New GIDL/GISL */
    bool  BSIM4kgislGiven   = false;    /* v4.7 New GIDL/GISL */
    bool  BSIM4rgislGiven   = false;    /* v4.7 New GIDL/GISL */
    bool  BSIM4aigcGiven   = false;
    bool  BSIM4bigcGiven   = false;
    bool  BSIM4cigcGiven   = false;
    bool  BSIM4aigsdGiven   = false;
    bool  BSIM4bigsdGiven   = false;
    bool  BSIM4cigsdGiven   = false;
    bool  BSIM4aigsGiven   = false;
    bool  BSIM4bigsGiven   = false;
    bool  BSIM4cigsGiven   = false;
    bool  BSIM4aigdGiven   = false;
    bool  BSIM4bigdGiven   = false;
    bool  BSIM4cigdGiven   = false;
    bool  BSIM4aigbaccGiven   = false;
    bool  BSIM4bigbaccGiven   = false;
    bool  BSIM4cigbaccGiven   = false;
    bool  BSIM4aigbinvGiven   = false;
    bool  BSIM4bigbinvGiven   = false;
    bool  BSIM4cigbinvGiven   = false;
    bool  BSIM4nigcGiven   = false;
    bool  BSIM4nigbinvGiven   = false;
    bool  BSIM4nigbaccGiven   = false;
    bool  BSIM4ntoxGiven   = false;
    bool  BSIM4eigbinvGiven   = false;
    bool  BSIM4pigcdGiven   = false;
    bool  BSIM4poxedgeGiven   = false;
    bool  BSIM4ijthdfwdGiven  = false;
    bool  BSIM4ijthsfwdGiven  = false;
    bool  BSIM4ijthdrevGiven  = false;
    bool  BSIM4ijthsrevGiven  = false;
    bool  BSIM4xjbvdGiven   = false;
    bool  BSIM4xjbvsGiven   = false;
    bool  BSIM4bvdGiven   = false;
    bool  BSIM4bvsGiven   = false;

    bool  BSIM4jtssGiven   = false;
    bool  BSIM4jtsdGiven   = false;
    bool  BSIM4jtsswsGiven   = false;
    bool  BSIM4jtsswdGiven   = false;
    bool  BSIM4jtsswgsGiven   = false;
    bool  BSIM4jtsswgdGiven   = false;
    bool  BSIM4jtweffGiven    = false;
    bool  BSIM4njtsGiven   = false;
    bool  BSIM4njtsswGiven   = false;
    bool  BSIM4njtsswgGiven   = false;
    bool  BSIM4njtsdGiven   = false;
    bool  BSIM4njtsswdGiven   = false;
    bool  BSIM4njtsswgdGiven   = false;
    bool  BSIM4xtssGiven   = false;
    bool  BSIM4xtsdGiven   = false;
    bool  BSIM4xtsswsGiven   = false;
    bool  BSIM4xtsswdGiven   = false;
    bool  BSIM4xtsswgsGiven   = false;
    bool  BSIM4xtsswgdGiven   = false;
    bool  BSIM4tnjtsGiven   = false;
    bool  BSIM4tnjtsswGiven   = false;
    bool  BSIM4tnjtsswgGiven   = false;
    bool  BSIM4tnjtsdGiven   = false;
    bool  BSIM4tnjtsswdGiven   = false;
    bool  BSIM4tnjtsswgdGiven   = false;
    bool  BSIM4vtssGiven   = false;
    bool  BSIM4vtsdGiven   = false;
    bool  BSIM4vtsswsGiven   = false;
    bool  BSIM4vtsswdGiven   = false;
    bool  BSIM4vtsswgsGiven   = false;
    bool  BSIM4vtsswgdGiven   = false;

    bool  BSIM4vfbGiven   = false;
    bool  BSIM4gbminGiven = false;
    bool  BSIM4rbdbGiven = false;
    bool  BSIM4rbsbGiven = false;
    bool  BSIM4rbpsGiven = false;
    bool  BSIM4rbpdGiven = false;
    bool  BSIM4rbpbGiven = false;

    bool BSIM4rbps0Given = false;
    bool BSIM4rbpslGiven = false;
    bool BSIM4rbpswGiven = false;
    bool BSIM4rbpsnfGiven = false;

    bool BSIM4rbpd0Given = false;
    bool BSIM4rbpdlGiven = false;
    bool BSIM4rbpdwGiven = false;
    bool BSIM4rbpdnfGiven = false;

    bool BSIM4rbpbx0Given = false;
    bool BSIM4rbpbxlGiven = false;
    bool BSIM4rbpbxwGiven = false;
    bool BSIM4rbpbxnfGiven = false;
    bool BSIM4rbpby0Given = false;
    bool BSIM4rbpbylGiven = false;
    bool BSIM4rbpbywGiven = false;
    bool BSIM4rbpbynfGiven = false;

    bool BSIM4rbsbx0Given = false;
    bool BSIM4rbsby0Given = false;
    bool BSIM4rbdbx0Given = false;
    bool BSIM4rbdby0Given = false;

    bool BSIM4rbsdbxlGiven = false;
    bool BSIM4rbsdbxwGiven = false;
    bool BSIM4rbsdbxnfGiven = false;
    bool BSIM4rbsdbylGiven = false;
    bool BSIM4rbsdbywGiven = false;
    bool BSIM4rbsdbynfGiven = false;

    bool  BSIM4xrcrg1Given   = false;
    bool  BSIM4xrcrg2Given   = false;
    bool  BSIM4tnoiaGiven    = false;
    bool  BSIM4tnoibGiven    = false;
    bool  BSIM4tnoicGiven    = false;
    bool  BSIM4rnoiaGiven    = false;
    bool  BSIM4rnoibGiven    = false;
    bool  BSIM4rnoicGiven    = false;
    bool  BSIM4ntnoiGiven    = false;

    bool  BSIM4lambdaGiven    = false;
    bool  BSIM4vtlGiven    = false;
    bool  BSIM4lcGiven    = false;
    bool  BSIM4xnGiven    = false;
    bool  BSIM4vfbsdoffGiven    = false;
    bool  BSIM4lintnoiGiven    = false;
    bool  BSIM4tvfbsdoffGiven    = false;

    /* CV model and parasitics */
    bool  BSIM4cgslGiven   = false;
    bool  BSIM4cgdlGiven   = false;
    bool  BSIM4ckappasGiven   = false;
    bool  BSIM4ckappadGiven   = false;
    bool  BSIM4cfGiven   = false;
    bool  BSIM4vfbcvGiven   = false;
    bool  BSIM4clcGiven   = false;
    bool  BSIM4cleGiven   = false;
    bool  BSIM4dwcGiven   = false;
    bool  BSIM4dlcGiven   = false;
    bool  BSIM4xwGiven    = false;
    bool  BSIM4xlGiven    = false;
    bool  BSIM4dlcigGiven   = false;
    bool  BSIM4dlcigdGiven   = false;
    bool  BSIM4dwjGiven   = false;
    bool  BSIM4noffGiven  = false;
    bool  BSIM4voffcvGiven = false;
    bool  BSIM4acdeGiven  = false;
    bool  BSIM4moinGiven  = false;
    bool  BSIM4tcjGiven   = false;
    bool  BSIM4tcjswGiven = false;
    bool  BSIM4tcjswgGiven = false;
    bool  BSIM4tpbGiven    = false;
    bool  BSIM4tpbswGiven  = false;
    bool  BSIM4tpbswgGiven = false;
    bool  BSIM4dmcgGiven = false;
    bool  BSIM4dmciGiven = false;
    bool  BSIM4dmdgGiven = false;
    bool  BSIM4dmcgtGiven = false;
    bool  BSIM4xgwGiven = false;
    bool  BSIM4xglGiven = false;
    bool  BSIM4rshgGiven = false;
    bool  BSIM4ngconGiven = false;


    /* Length dependence */
    bool  BSIM4lcdscGiven   = false;
    bool  BSIM4lcdscbGiven   = false;
    bool  BSIM4lcdscdGiven   = false;
    bool  BSIM4lcitGiven   = false;
    bool  BSIM4lnfactorGiven   = false;
    bool  BSIM4lxjGiven   = false;
    bool  BSIM4lvsatGiven   = false;
    bool  BSIM4latGiven   = false;
    bool  BSIM4la0Given   = false;
    bool  BSIM4lagsGiven   = false;
    bool  BSIM4la1Given   = false;
    bool  BSIM4la2Given   = false;
    bool  BSIM4lketaGiven   = false;
    bool  BSIM4lnsubGiven   = false;
    bool  BSIM4lndepGiven   = false;
    bool  BSIM4lnsdGiven    = false;
    bool  BSIM4lphinGiven   = false;
    bool  BSIM4lngateGiven   = false;
    bool  BSIM4lgamma1Given   = false;
    bool  BSIM4lgamma2Given   = false;
    bool  BSIM4lvbxGiven   = false;
    bool  BSIM4lvbmGiven   = false;
    bool  BSIM4lxtGiven   = false;
    bool  BSIM4lk1Given   = false;
    bool  BSIM4lkt1Given   = false;
    bool  BSIM4lkt1lGiven   = false;
    bool  BSIM4lkt2Given   = false;
    bool  BSIM4lk2Given   = false;
    bool  BSIM4lk3Given   = false;
    bool  BSIM4lk3bGiven   = false;
    bool  BSIM4lw0Given   = false;
    bool  BSIM4ldvtp0Given = false;
    bool  BSIM4ldvtp1Given = false;
    bool  BSIM4ldvtp2Given = false;  /* New DIBL/Rout */
    bool  BSIM4ldvtp3Given = false;
    bool  BSIM4ldvtp4Given = false;
    bool  BSIM4ldvtp5Given = false;
    bool  BSIM4llpe0Given   = false;
    bool  BSIM4llpebGiven   = false;
    bool  BSIM4ldvt0Given   = false;
    bool  BSIM4ldvt1Given   = false;
    bool  BSIM4ldvt2Given   = false;
    bool  BSIM4ldvt0wGiven   = false;
    bool  BSIM4ldvt1wGiven   = false;
    bool  BSIM4ldvt2wGiven   = false;
    bool  BSIM4ldroutGiven   = false;
    bool  BSIM4ldsubGiven   = false;
    bool  BSIM4lvth0Given   = false;
    bool  BSIM4luaGiven   = false;
    bool  BSIM4lua1Given   = false;
    bool  BSIM4lubGiven   = false;
    bool  BSIM4lub1Given   = false;
    bool  BSIM4lucGiven   = false;
    bool  BSIM4luc1Given   = false;
    bool  BSIM4ludGiven     = false;
    bool  BSIM4lud1Given     = false;
    bool  BSIM4lupGiven     = false;
    bool  BSIM4llpGiven     = false;
    bool  BSIM4lu0Given   = false;
    bool  BSIM4leuGiven   = false;
    bool  BSIM4lucsGiven   = false;
    bool  BSIM4luteGiven   = false;
    bool  BSIM4lucsteGiven  = false;
    bool  BSIM4lvoffGiven   = false;
    bool  BSIM4ltvoffGiven   = false;
    bool  BSIM4ltnfactorGiven  = false;  /* v4.7 Temp dep of leakage current */
    bool  BSIM4lteta0Given   = false;    /* v4.7 temp dep of leakage current */
    bool  BSIM4ltvoffcvGiven   = false;  /* v4.7 temp dep of leakage current */
    bool  BSIM4lminvGiven   = false;
    bool  BSIM4lminvcvGiven   = false;
    bool  BSIM4lrdswGiven   = false;
    bool  BSIM4lrswGiven   = false;
    bool  BSIM4lrdwGiven   = false;
    bool  BSIM4lprwgGiven   = false;
    bool  BSIM4lprwbGiven   = false;
    bool  BSIM4lprtGiven   = false;
    bool  BSIM4leta0Given   = false;
    bool  BSIM4letabGiven   = false;
    bool  BSIM4lpclmGiven   = false;
    bool  BSIM4lpdibl1Given   = false;
    bool  BSIM4lpdibl2Given   = false;
    bool  BSIM4lpdiblbGiven   = false;
    bool  BSIM4lfproutGiven   = false;
    bool  BSIM4lpditsGiven    = false;
    bool  BSIM4lpditsdGiven    = false;
    bool  BSIM4lpscbe1Given   = false;
    bool  BSIM4lpscbe2Given   = false;
    bool  BSIM4lpvagGiven   = false;
    bool  BSIM4ldeltaGiven  = false;
    bool  BSIM4lwrGiven   = false;
    bool  BSIM4ldwgGiven   = false;
    bool  BSIM4ldwbGiven   = false;
    bool  BSIM4lb0Given   = false;
    bool  BSIM4lb1Given   = false;
    bool  BSIM4lalpha0Given   = false;
    bool  BSIM4lalpha1Given   = false;
    bool  BSIM4lbeta0Given   = false;
    bool  BSIM4lvfbGiven   = false;
    bool  BSIM4lagidlGiven   = false;
    bool  BSIM4lbgidlGiven   = false;
    bool  BSIM4lcgidlGiven   = false;
    bool  BSIM4legidlGiven   = false;
    bool  BSIM4lfgidlGiven   = false;    /* v4.7 New GIDL/GISL */
    bool  BSIM4lkgidlGiven   = false;    /* v4.7 New GIDL/GISL */
    bool  BSIM4lrgidlGiven   = false;    /* v4.7 New GIDL/GISL */
    bool  BSIM4lagislGiven   = false;
    bool  BSIM4lbgislGiven   = false;
    bool  BSIM4lcgislGiven   = false;
    bool  BSIM4legislGiven   = false;
    bool  BSIM4lfgislGiven   = false;    /* v4.7 New GIDL/GISL */
    bool  BSIM4lkgislGiven   = false;    /* v4.7 New GIDL/GISL */
    bool  BSIM4lrgislGiven   = false;    /* v4.7 New GIDL/GISL */
    bool  BSIM4laigcGiven   = false;
    bool  BSIM4lbigcGiven   = false;
    bool  BSIM4lcigcGiven   = false;
    bool  BSIM4laigsdGiven   = false;
    bool  BSIM4lbigsdGiven   = false;
    bool  BSIM4lcigsdGiven   = false;
    bool  BSIM4laigsGiven   = false;
    bool  BSIM4lbigsGiven   = false;
    bool  BSIM4lcigsGiven   = false;
    bool  BSIM4laigdGiven   = false;
    bool  BSIM4lbigdGiven   = false;
    bool  BSIM4lcigdGiven   = false;
    bool  BSIM4laigbaccGiven   = false;
    bool  BSIM4lbigbaccGiven   = false;
    bool  BSIM4lcigbaccGiven   = false;
    bool  BSIM4laigbinvGiven   = false;
    bool  BSIM4lbigbinvGiven   = false;
    bool  BSIM4lcigbinvGiven   = false;
    bool  BSIM4lnigcGiven   = false;
    bool  BSIM4lnigbinvGiven   = false;
    bool  BSIM4lnigbaccGiven   = false;
    bool  BSIM4lntoxGiven   = false;
    bool  BSIM4leigbinvGiven   = false;
    bool  BSIM4lpigcdGiven   = false;
    bool  BSIM4lpoxedgeGiven   = false;
    bool  BSIM4lxrcrg1Given   = false;
    bool  BSIM4lxrcrg2Given   = false;
    bool  BSIM4llambdaGiven    = false;
    bool  BSIM4lvtlGiven    = false;
    bool  BSIM4lxnGiven    = false;
    bool  BSIM4lvfbsdoffGiven    = false;
    bool  BSIM4ltvfbsdoffGiven    = false;

    /* CV model */
    bool  BSIM4lcgslGiven   = false;
    bool  BSIM4lcgdlGiven   = false;
    bool  BSIM4lckappasGiven   = false;
    bool  BSIM4lckappadGiven   = false;
    bool  BSIM4lcfGiven   = false;
    bool  BSIM4lclcGiven   = false;
    bool  BSIM4lcleGiven   = false;
    bool  BSIM4lvfbcvGiven   = false;
    bool  BSIM4lnoffGiven   = false;
    bool  BSIM4lvoffcvGiven = false;
    bool  BSIM4lacdeGiven   = false;
    bool  BSIM4lmoinGiven   = false;

    /* Width dependence */
    bool  BSIM4wcdscGiven   = false;
    bool  BSIM4wcdscbGiven   = false;
    bool  BSIM4wcdscdGiven   = false;
    bool  BSIM4wcitGiven   = false;
    bool  BSIM4wnfactorGiven   = false;
    bool  BSIM4wxjGiven   = false;
    bool  BSIM4wvsatGiven   = false;
    bool  BSIM4watGiven   = false;
    bool  BSIM4wa0Given   = false;
    bool  BSIM4wagsGiven   = false;
    bool  BSIM4wa1Given   = false;
    bool  BSIM4wa2Given   = false;
    bool  BSIM4wketaGiven   = false;
    bool  BSIM4wnsubGiven   = false;
    bool  BSIM4wndepGiven   = false;
    bool  BSIM4wnsdGiven    = false;
    bool  BSIM4wphinGiven   = false;
    bool  BSIM4wngateGiven   = false;
    bool  BSIM4wgamma1Given   = false;
    bool  BSIM4wgamma2Given   = false;
    bool  BSIM4wvbxGiven   = false;
    bool  BSIM4wvbmGiven   = false;
    bool  BSIM4wxtGiven   = false;
    bool  BSIM4wk1Given   = false;
    bool  BSIM4wkt1Given   = false;
    bool  BSIM4wkt1lGiven   = false;
    bool  BSIM4wkt2Given   = false;
    bool  BSIM4wk2Given   = false;
    bool  BSIM4wk3Given   = false;
    bool  BSIM4wk3bGiven   = false;
    bool  BSIM4ww0Given   = false;
    bool  BSIM4wdvtp0Given = false;
    bool  BSIM4wdvtp1Given = false;
    bool  BSIM4wdvtp2Given = false;  /* New DIBL/Rout */
    bool  BSIM4wdvtp3Given = false;
    bool  BSIM4wdvtp4Given = false;
    bool  BSIM4wdvtp5Given = false;
    bool  BSIM4wlpe0Given   = false;
    bool  BSIM4wlpebGiven   = false;
    bool  BSIM4wdvt0Given   = false;
    bool  BSIM4wdvt1Given   = false;
    bool  BSIM4wdvt2Given   = false;
    bool  BSIM4wdvt0wGiven   = false;
    bool  BSIM4wdvt1wGiven   = false;
    bool  BSIM4wdvt2wGiven   = false;
    bool  BSIM4wdroutGiven   = false;
    bool  BSIM4wdsubGiven   = false;
    bool  BSIM4wvth0Given   = false;
    bool  BSIM4wuaGiven   = false;
    bool  BSIM4wua1Given   = false;
    bool  BSIM4wubGiven   = false;
    bool  BSIM4wub1Given   = false;
    bool  BSIM4wucGiven   = false;
    bool  BSIM4wuc1Given   = false;
    bool  BSIM4wudGiven     = false;
    bool  BSIM4wud1Given     = false;
    bool  BSIM4wupGiven     = false;
    bool  BSIM4wlpGiven     = false;
    bool  BSIM4wu0Given   = false;
    bool  BSIM4weuGiven   = false;
    bool  BSIM4wucsGiven  = false;
    bool  BSIM4wuteGiven   = false;
    bool  BSIM4wucsteGiven  = false;
    bool  BSIM4wvoffGiven   = false;
    bool  BSIM4wtvoffGiven   = false;
    bool  BSIM4wtnfactorGiven  = false;  /* v4.7 Temp dep of leakage current */
    bool  BSIM4wteta0Given   = false;    /* v4.7 temp dep of leakage current */
    bool  BSIM4wtvoffcvGiven   = false;  /* v4.7 temp dep of leakage current */
    bool  BSIM4wminvGiven   = false;
    bool  BSIM4wminvcvGiven   = false;
    bool  BSIM4wrdswGiven   = false;
    bool  BSIM4wrswGiven   = false;
    bool  BSIM4wrdwGiven   = false;
    bool  BSIM4wprwgGiven   = false;
    bool  BSIM4wprwbGiven   = false;
    bool  BSIM4wprtGiven   = false;
    bool  BSIM4weta0Given   = false;
    bool  BSIM4wetabGiven   = false;
    bool  BSIM4wpclmGiven   = false;
    bool  BSIM4wpdibl1Given   = false;
    bool  BSIM4wpdibl2Given   = false;
    bool  BSIM4wpdiblbGiven   = false;
    bool  BSIM4wfproutGiven   = false;
    bool  BSIM4wpditsGiven    = false;
    bool  BSIM4wpditsdGiven    = false;
    bool  BSIM4wpscbe1Given   = false;
    bool  BSIM4wpscbe2Given   = false;
    bool  BSIM4wpvagGiven   = false;
    bool  BSIM4wdeltaGiven  = false;
    bool  BSIM4wwrGiven   = false;
    bool  BSIM4wdwgGiven   = false;
    bool  BSIM4wdwbGiven   = false;
    bool  BSIM4wb0Given   = false;
    bool  BSIM4wb1Given   = false;
    bool  BSIM4walpha0Given   = false;
    bool  BSIM4walpha1Given   = false;
    bool  BSIM4wbeta0Given   = false;
    bool  BSIM4wvfbGiven   = false;
    bool  BSIM4wagidlGiven   = false;
    bool  BSIM4wbgidlGiven   = false;
    bool  BSIM4wcgidlGiven   = false;
    bool  BSIM4wegidlGiven   = false;
    bool  BSIM4wfgidlGiven   = false;    /* v4.7 New GIDL/GISL */
    bool  BSIM4wkgidlGiven   = false;    /* v4.7 New GIDL/GISL */
    bool  BSIM4wrgidlGiven   = false;    /* v4.7 New GIDL/GISL */
    bool  BSIM4wagislGiven   = false;
    bool  BSIM4wbgislGiven   = false;
    bool  BSIM4wcgislGiven   = false;
    bool  BSIM4wegislGiven   = false;
    bool  BSIM4wfgislGiven   = false;    /* v4.7 New GIDL/GISL */
    bool  BSIM4wkgislGiven   = false;    /* v4.7 New GIDL/GISL */
    bool  BSIM4wrgislGiven   = false;    /* v4.7 New GIDL/GISL */
    bool  BSIM4waigcGiven   = false;
    bool  BSIM4wbigcGiven   = false;
    bool  BSIM4wcigcGiven   = false;
    bool  BSIM4waigsdGiven   = false;
    bool  BSIM4wbigsdGiven   = false;
    bool  BSIM4wcigsdGiven   = false;
    bool  BSIM4waigsGiven   = false;
    bool  BSIM4wbigsGiven   = false;
    bool  BSIM4wcigsGiven   = false;
    bool  BSIM4waigdGiven   = false;
    bool  BSIM4wbigdGiven   = false;
    bool  BSIM4wcigdGiven   = false;
    bool  BSIM4waigbaccGiven   = false;
    bool  BSIM4wbigbaccGiven   = false;
    bool  BSIM4wcigbaccGiven   = false;
    bool  BSIM4waigbinvGiven   = false;
    bool  BSIM4wbigbinvGiven   = false;
    bool  BSIM4wcigbinvGiven   = false;
    bool  BSIM4wnigcGiven   = false;
    bool  BSIM4wnigbinvGiven   = false;
    bool  BSIM4wnigbaccGiven   = false;
    bool  BSIM4wntoxGiven   = false;
    bool  BSIM4weigbinvGiven   = false;
    bool  BSIM4wpigcdGiven   = false;
    bool  BSIM4wpoxedgeGiven   = false;
    bool  BSIM4wxrcrg1Given   = false;
    bool  BSIM4wxrcrg2Given   = false;
    bool  BSIM4wlambdaGiven    = false;
    bool  BSIM4wvtlGiven    = false;
    bool  BSIM4wxnGiven    = false;
    bool  BSIM4wvfbsdoffGiven    = false;
    bool  BSIM4wtvfbsdoffGiven    = false;

    /* CV model */
    bool  BSIM4wcgslGiven   = false;
    bool  BSIM4wcgdlGiven   = false;
    bool  BSIM4wckappasGiven   = false;
    bool  BSIM4wckappadGiven   = false;
    bool  BSIM4wcfGiven   = false;
    bool  BSIM4wclcGiven   = false;
    bool  BSIM4wcleGiven   = false;
    bool  BSIM4wvfbcvGiven   = false;
    bool  BSIM4wnoffGiven   = false;
    bool  BSIM4wvoffcvGiven = false;
    bool  BSIM4wacdeGiven   = false;
    bool  BSIM4wmoinGiven   = false;

    /* Cross-term dependence */
    bool  BSIM4pcdscGiven   = false;
    bool  BSIM4pcdscbGiven   = false;
    bool  BSIM4pcdscdGiven   = false;
    bool  BSIM4pcitGiven   = false;
    bool  BSIM4pnfactorGiven   = false;
    bool  BSIM4pxjGiven   = false;
    bool  BSIM4pvsatGiven   = false;
    bool  BSIM4patGiven   = false;
    bool  BSIM4pa0Given   = false;
    bool  BSIM4pagsGiven   = false;
    bool  BSIM4pa1Given   = false;
    bool  BSIM4pa2Given   = false;
    bool  BSIM4pketaGiven   = false;
    bool  BSIM4pnsubGiven   = false;
    bool  BSIM4pndepGiven   = false;
    bool  BSIM4pnsdGiven    = false;
    bool  BSIM4pphinGiven   = false;
    bool  BSIM4pngateGiven   = false;
    bool  BSIM4pgamma1Given   = false;
    bool  BSIM4pgamma2Given   = false;
    bool  BSIM4pvbxGiven   = false;
    bool  BSIM4pvbmGiven   = false;
    bool  BSIM4pxtGiven   = false;
    bool  BSIM4pk1Given   = false;
    bool  BSIM4pkt1Given   = false;
    bool  BSIM4pkt1lGiven   = false;
    bool  BSIM4pkt2Given   = false;
    bool  BSIM4pk2Given   = false;
    bool  BSIM4pk3Given   = false;
    bool  BSIM4pk3bGiven   = false;
    bool  BSIM4pw0Given   = false;
    bool  BSIM4pdvtp0Given = false;
    bool  BSIM4pdvtp1Given = false;
    bool  BSIM4pdvtp2Given = false;  /* New DIBL/Rout */
    bool  BSIM4pdvtp3Given = false;
    bool  BSIM4pdvtp4Given = false;
    bool  BSIM4pdvtp5Given = false;
    bool  BSIM4plpe0Given   = false;
    bool  BSIM4plpebGiven   = false;
    bool  BSIM4pdvt0Given   = false;
    bool  BSIM4pdvt1Given   = false;
    bool  BSIM4pdvt2Given   = false;
    bool  BSIM4pdvt0wGiven   = false;
    bool  BSIM4pdvt1wGiven   = false;
    bool  BSIM4pdvt2wGiven   = false;
    bool  BSIM4pdroutGiven   = false;
    bool  BSIM4pdsubGiven   = false;
    bool  BSIM4pvth0Given   = false;
    bool  BSIM4puaGiven   = false;
    bool  BSIM4pua1Given   = false;
    bool  BSIM4pubGiven   = false;
    bool  BSIM4pub1Given   = false;
    bool  BSIM4pucGiven   = false;
    bool  BSIM4puc1Given   = false;
    bool  BSIM4pudGiven     = false;
    bool  BSIM4pud1Given     = false;
    bool  BSIM4pupGiven     = false;
    bool  BSIM4plpGiven     = false;
    bool  BSIM4pu0Given   = false;
    bool  BSIM4peuGiven   = false;
    bool  BSIM4pucsGiven   = false;
    bool  BSIM4puteGiven   = false;
    bool  BSIM4pucsteGiven  = false;
    bool  BSIM4pvoffGiven   = false;
    bool  BSIM4ptvoffGiven   = false;
    bool  BSIM4ptnfactorGiven  = false;  /* v4.7 Temp dep of leakage current */
    bool  BSIM4pteta0Given   = false;    /* v4.7 temp dep of leakage current */
    bool  BSIM4ptvoffcvGiven   = false;  /* v4.7 temp dep of leakage current */
    bool  BSIM4pminvGiven   = false;
    bool  BSIM4pminvcvGiven   = false;
    bool  BSIM4prdswGiven   = false;
    bool  BSIM4prswGiven   = false;
    bool  BSIM4prdwGiven   = false;
    bool  BSIM4pprwgGiven   = false;
    bool  BSIM4pprwbGiven   = false;
    bool  BSIM4pprtGiven   = false;
    bool  BSIM4peta0Given   = false;
    bool  BSIM4petabGiven   = false;
    bool  BSIM4ppclmGiven   = false;
    bool  BSIM4ppdibl1Given   = false;
    bool  BSIM4ppdibl2Given   = false;
    bool  BSIM4ppdiblbGiven   = false;
    bool  BSIM4pfproutGiven   = false;
    bool  BSIM4ppditsGiven    = false;
    bool  BSIM4ppditsdGiven    = false;
    bool  BSIM4ppscbe1Given   = false;
    bool  BSIM4ppscbe2Given   = false;
    bool  BSIM4ppvagGiven   = false;
    bool  BSIM4pdeltaGiven  = false;
    bool  BSIM4pwrGiven   = false;
    bool  BSIM4pdwgGiven   = false;
    bool  BSIM4pdwbGiven   = false;
    bool  BSIM4pb0Given   = false;
    bool  BSIM4pb1Given   = false;
    bool  BSIM4palpha0Given   = false;
    bool  BSIM4palpha1Given   = false;
    bool  BSIM4pbeta0Given   = false;
    bool  BSIM4pvfbGiven   = false;
    bool  BSIM4pagidlGiven   = false;
    bool  BSIM4pbgidlGiven   = false;
    bool  BSIM4pcgidlGiven   = false;
    bool  BSIM4pegidlGiven   = false;
    bool  BSIM4pfgidlGiven   = false;    /* v4.7 New GIDL/GISL */
    bool  BSIM4pkgidlGiven   = false;    /* v4.7 New GIDL/GISL */
    bool  BSIM4prgidlGiven   = false;    /* v4.7 New GIDL/GISL */
    bool  BSIM4pagislGiven   = false;
    bool  BSIM4pbgislGiven   = false;
    bool  BSIM4pcgislGiven   = false;
    bool  BSIM4pegislGiven   = false;
    bool  BSIM4pfgislGiven   = false;    /* v4.7 New GIDL/GISL */
    bool  BSIM4pkgislGiven   = false;    /* v4.7 New GIDL/GISL */
    bool  BSIM4prgislGiven   = false;    /* v4.7 New GIDL/GISL */
    bool  BSIM4paigcGiven   = false;
    bool  BSIM4pbigcGiven   = false;
    bool  BSIM4pcigcGiven   = false;
    bool  BSIM4paigsdGiven   = false;
    bool  BSIM4pbigsdGiven   = false;
    bool  BSIM4pcigsdGiven   = false;
    bool  BSIM4paigsGiven   = false;
    bool  BSIM4pbigsGiven   = false;
    bool  BSIM4pcigsGiven   = false;
    bool  BSIM4paigdGiven   = false;
    bool  BSIM4pbigdGiven   = false;
    bool  BSIM4pcigdGiven   = false;
    bool  BSIM4paigbaccGiven   = false;
    bool  BSIM4pbigbaccGiven   = false;
    bool  BSIM4pcigbaccGiven   = false;
    bool  BSIM4paigbinvGiven   = false;
    bool  BSIM4pbigbinvGiven   = false;
    bool  BSIM4pcigbinvGiven   = false;
    bool  BSIM4pnigcGiven   = false;
    bool  BSIM4pnigbinvGiven   = false;
    bool  BSIM4pnigbaccGiven   = false;
    bool  BSIM4pntoxGiven   = false;
    bool  BSIM4peigbinvGiven   = false;
    bool  BSIM4ppigcdGiven   = false;
    bool  BSIM4ppoxedgeGiven   = false;
    bool  BSIM4pxrcrg1Given   = false;
    bool  BSIM4pxrcrg2Given   = false;
    bool  BSIM4plambdaGiven    = false;
    bool  BSIM4pvtlGiven    = false;
    bool  BSIM4pxnGiven    = false;
    bool  BSIM4pvfbsdoffGiven    = false;
    bool  BSIM4ptvfbsdoffGiven    = false;

    /* CV model */
    bool  BSIM4pcgslGiven   = false;
    bool  BSIM4pcgdlGiven   = false;
    bool  BSIM4pckappasGiven   = false;
    bool  BSIM4pckappadGiven   = false;
    bool  BSIM4pcfGiven   = false;
    bool  BSIM4pclcGiven   = false;
    bool  BSIM4pcleGiven   = false;
    bool  BSIM4pvfbcvGiven   = false;
    bool  BSIM4pnoffGiven   = false;
    bool  BSIM4pvoffcvGiven = false;
    bool  BSIM4pacdeGiven   = false;
    bool  BSIM4pmoinGiven   = false;

    bool  BSIM4useFringeGiven   = false;

    bool  BSIM4tnomGiven   = false;
    bool  BSIM4cgsoGiven   = false;
    bool  BSIM4cgdoGiven   = false;
    bool  BSIM4cgboGiven   = false;
    bool  BSIM4xpartGiven   = false;
    bool  BSIM4sheetResistanceGiven   = false;

    bool  BSIM4SjctSatCurDensityGiven   = false;
    bool  BSIM4SjctSidewallSatCurDensityGiven   = false;
    bool  BSIM4SjctGateSidewallSatCurDensityGiven   = false;
    bool  BSIM4SbulkJctPotentialGiven   = false;
    bool  BSIM4SbulkJctBotGradingCoeffGiven   = false;
    bool  BSIM4SsidewallJctPotentialGiven   = false;
    bool  BSIM4SGatesidewallJctPotentialGiven   = false;
    bool  BSIM4SbulkJctSideGradingCoeffGiven   = false;
    bool  BSIM4SunitAreaJctCapGiven   = false;
    bool  BSIM4SunitLengthSidewallJctCapGiven   = false;
    bool  BSIM4SbulkJctGateSideGradingCoeffGiven   = false;
    bool  BSIM4SunitLengthGateSidewallJctCapGiven   = false;
    bool  BSIM4SjctEmissionCoeffGiven = false;
    bool  BSIM4SjctTempExponentGiven    = false;

    bool  BSIM4DjctSatCurDensityGiven   = false;
    bool  BSIM4DjctSidewallSatCurDensityGiven   = false;
    bool  BSIM4DjctGateSidewallSatCurDensityGiven   = false;
    bool  BSIM4DbulkJctPotentialGiven   = false;
    bool  BSIM4DbulkJctBotGradingCoeffGiven   = false;
    bool  BSIM4DsidewallJctPotentialGiven   = false;
    bool  BSIM4DGatesidewallJctPotentialGiven   = false;
    bool  BSIM4DbulkJctSideGradingCoeffGiven   = false;
    bool  BSIM4DunitAreaJctCapGiven   = false;
    bool  BSIM4DunitLengthSidewallJctCapGiven   = false;
    bool  BSIM4DbulkJctGateSideGradingCoeffGiven   = false;
    bool  BSIM4DunitLengthGateSidewallJctCapGiven   = false;
    bool  BSIM4DjctEmissionCoeffGiven = false;
    bool  BSIM4DjctTempExponentGiven = false;

    bool  BSIM4oxideTrapDensityAGiven  = false;
    bool  BSIM4oxideTrapDensityBGiven  = false;
    bool  BSIM4oxideTrapDensityCGiven  = false;
    bool  BSIM4emGiven  = false;
    bool  BSIM4efGiven  = false;
    bool  BSIM4afGiven  = false;
    bool  BSIM4kfGiven  = false;

    bool  BSIM4LintGiven   = false;
    bool  BSIM4LlGiven   = false;
    bool  BSIM4LlcGiven   = false;
    bool  BSIM4LlnGiven   = false;
    bool  BSIM4LwGiven   = false;
    bool  BSIM4LwcGiven   = false;
    bool  BSIM4LwnGiven   = false;
    bool  BSIM4LwlGiven   = false;
    bool  BSIM4LwlcGiven   = false;
    bool  BSIM4LminGiven   = false;
    bool  BSIM4LmaxGiven   = false;

    bool  BSIM4WintGiven   = false;
    bool  BSIM4WlGiven   = false;
    bool  BSIM4WlcGiven   = false;
    bool  BSIM4WlnGiven   = false;
    bool  BSIM4WwGiven   = false;
    bool  BSIM4WwcGiven   = false;
    bool  BSIM4WwnGiven   = false;
    bool  BSIM4WwlGiven   = false;
    bool  BSIM4WwlcGiven   = false;
    bool  BSIM4WminGiven   = false;
    bool  BSIM4WmaxGiven   = false;

    /* added for stress effect */
    bool  BSIM4sarefGiven   = false;
    bool  BSIM4sbrefGiven   = false;
    bool  BSIM4wlodGiven  = false;
    bool  BSIM4ku0Given   = false;
    bool  BSIM4kvsatGiven  = false;
    bool  BSIM4kvth0Given  = false;
    bool  BSIM4tku0Given   = false;
    bool  BSIM4llodku0Given   = false;
    bool  BSIM4wlodku0Given   = false;
    bool  BSIM4llodvthGiven   = false;
    bool  BSIM4wlodvthGiven   = false;
    bool  BSIM4lku0Given   = false;
    bool  BSIM4wku0Given   = false;
    bool  BSIM4pku0Given   = false;
    bool  BSIM4lkvth0Given   = false;
    bool  BSIM4wkvth0Given   = false;
    bool  BSIM4pkvth0Given   = false;
    bool  BSIM4stk2Given   = false;
    bool  BSIM4lodk2Given  = false;
    bool  BSIM4steta0Given = false;
    bool  BSIM4lodeta0Given = false;

    bool  BSIM4webGiven   = false;
    bool  BSIM4wecGiven   = false;
    bool  BSIM4kvth0weGiven   = false;
    bool  BSIM4k2weGiven   = false;
    bool  BSIM4ku0weGiven   = false;
    bool  BSIM4screfGiven   = false;
    bool  BSIM4wpemodGiven   = false;
    bool  BSIM4lkvth0weGiven   = false;
    bool  BSIM4lk2weGiven   = false;
    bool  BSIM4lku0weGiven   = false;
    bool  BSIM4wkvth0weGiven   = false;
    bool  BSIM4wk2weGiven   = false;
    bool  BSIM4wku0weGiven   = false;
    bool  BSIM4pkvth0weGiven   = false;
    bool  BSIM4pk2weGiven   = false;
    bool  BSIM4pku0weGiven   = false;
    bool  BSIM4gidlclampGiven   = false;
    bool  BSIM4idovvdscGiven   = false;

    // If nomodcheck is true, skip the model check
    bool modcheck = false;
};

struct BSIM4V82{
    // BSIM4V82(std::string name): BSIM4name(name) {};
    // struct sBSIM4model *BSIM4modPtr;
    std::shared_ptr<BSIM4model> BSIM4modPtr = nullptr;

    std::string BSIM4name;
    // int BSIM4states;     /* index into state table for this device */
    std::array<double, 29> BSIM4states0;
    std::array<double, 29> BSIM4states1;
    int BSIM4dNode;
    int BSIM4gNodeExt;
    int BSIM4sNode;
    int BSIM4bNode;
    int BSIM4dNodePrime;
    int BSIM4gNodePrime;
    int BSIM4gNodeMid;
    int BSIM4sNodePrime;
    int BSIM4bNodePrime;
    int BSIM4dbNode;
    int BSIM4sbNode;
    int BSIM4qNode;

    double BSIM4ueff;
    double BSIM4thetavth;
    double BSIM4von;
    double BSIM4vdsat;
    double BSIM4cgdo;
    double BSIM4qgdo;
    double BSIM4cgso;
    double BSIM4qgso;
    double BSIM4grbsb;
    double BSIM4grbdb;
    double BSIM4grbpb;
    double BSIM4grbps;
    double BSIM4grbpd;

    double BSIM4vjsmFwd;
    double BSIM4vjsmRev;
    double BSIM4vjdmFwd;
    double BSIM4vjdmRev;
    double BSIM4XExpBVS;
    double BSIM4XExpBVD;
    double BSIM4SslpFwd;
    double BSIM4SslpRev;
    double BSIM4DslpFwd;
    double BSIM4DslpRev;
    double BSIM4IVjsmFwd;
    double BSIM4IVjsmRev;
    double BSIM4IVjdmFwd;
    double BSIM4IVjdmRev;

    double BSIM4grgeltd;
    double BSIM4Pseff;
    double BSIM4Pdeff;
    double BSIM4Aseff;
    double BSIM4Adeff;

    double BSIM4l;
    double BSIM4w;
    double BSIM4drainArea;
    double BSIM4sourceArea;
    double BSIM4drainSquares;
    double BSIM4sourceSquares;
    double BSIM4drainPerimeter;
    double BSIM4sourcePerimeter;
    double BSIM4sourceConductance;
    double BSIM4drainConductance;

    /* stress effect instance param */
    double BSIM4sa;
    double BSIM4sb;
    double BSIM4sd;
    double BSIM4sca;
    double BSIM4scb;
    double BSIM4scc;
    double BSIM4sc;

    double BSIM4rbdb;
    double BSIM4rbsb;
    double BSIM4rbpb;
    double BSIM4rbps;
    double BSIM4rbpd;

    double BSIM4delvto;
    double BSIM4xgw;
    double BSIM4ngcon;

    /* added here to account stress effect instance dependence */

    double BSIM4u0temp;
    double BSIM4vsattemp;
    double BSIM4vth0;
    double BSIM4vfb;
    double BSIM4vfbzb;
    double BSIM4vtfbphi1;
    double BSIM4vtfbphi2;
    double BSIM4k2;
    double BSIM4vbsc;
    double BSIM4k2ox;
    double BSIM4eta0;
    double BSIM4toxp;
    double BSIM4coxp;

    double BSIM4icVDS;
    double BSIM4icVGS;
    double BSIM4icVBS;
    double BSIM4nf;
    int BSIM4off;

    int BSIM4mode;
    int BSIM4trnqsMod;
    int BSIM4acnqsMod;
    int BSIM4rbodyMod;
    int BSIM4rgateMod;
    int BSIM4geoMod;
    int BSIM4rgeoMod;
    int BSIM4min;


    /* OP point */
    double BSIM4Vgsteff;
    double BSIM4vgs_eff;
    double BSIM4vgd_eff;
    double BSIM4dvgs_eff_dvg;
    double BSIM4dvgd_eff_dvg;
    double BSIM4Vdseff;
    double BSIM4nstar;
    double BSIM4Abulk;
    double BSIM4EsatL;
    double BSIM4AbovVgst2Vtm;
    double BSIM4qinv;
    double BSIM4cd;
    double BSIM4cbs;
    double BSIM4cbd;
    double BSIM4csub;
    double BSIM4Igidl;
    double BSIM4Igisl;
    double BSIM4gm;
    double BSIM4gds;
    double BSIM4gmbs;
    double BSIM4gbd;
    double BSIM4gbs;
    double BSIM4noiGd0;   /* tnoiMod=2 (v4.7) */
    double BSIM4Coxeff;

    double BSIM4gbbs;
    double BSIM4gbgs;
    double BSIM4gbds;
    double BSIM4ggidld;
    double BSIM4ggidlg;
    double BSIM4ggidls;
    double BSIM4ggidlb;
    double BSIM4ggisld;
    double BSIM4ggislg;
    double BSIM4ggisls;
    double BSIM4ggislb;

    double BSIM4Igcs;
    double BSIM4gIgcsg;
    double BSIM4gIgcsd;
    double BSIM4gIgcss;
    double BSIM4gIgcsb;
    double BSIM4Igcd;
    double BSIM4gIgcdg;
    double BSIM4gIgcdd;
    double BSIM4gIgcds;
    double BSIM4gIgcdb;

    double BSIM4Igs;
    double BSIM4gIgsg;
    double BSIM4gIgss;
    double BSIM4Igd;
    double BSIM4gIgdg;
    double BSIM4gIgdd;

    double BSIM4Igb;
    double BSIM4gIgbg;
    double BSIM4gIgbd;
    double BSIM4gIgbs;
    double BSIM4gIgbb;

    double BSIM4grdsw;
    double BSIM4IdovVds;
    double BSIM4gcrg;
    double BSIM4gcrgd;
    double BSIM4gcrgg;
    double BSIM4gcrgs;
    double BSIM4gcrgb;

    double BSIM4gstot;
    double BSIM4gstotd;
    double BSIM4gstotg;
    double BSIM4gstots;
    double BSIM4gstotb;

    double BSIM4gdtot;
    double BSIM4gdtotd;
    double BSIM4gdtotg;
    double BSIM4gdtots;
    double BSIM4gdtotb;

    double BSIM4cggb;
    double BSIM4cgdb;
    double BSIM4cgsb;
    double BSIM4cbgb;
    double BSIM4cbdb;
    double BSIM4cbsb;
    double BSIM4cdgb;
    double BSIM4cddb;
    double BSIM4cdsb;
    double BSIM4csgb;
    double BSIM4csdb;
    double BSIM4cssb;
    double BSIM4cgbb;
    double BSIM4cdbb;
    double BSIM4csbb;
    double BSIM4cbbb;
    double BSIM4capbd;
    double BSIM4capbs;

    double BSIM4cqgb;
    double BSIM4cqdb;
    double BSIM4cqsb;
    double BSIM4cqbb;

    double BSIM4qgate;
    double BSIM4qbulk;
    double BSIM4qdrn;
    double BSIM4qsrc;
    double BSIM4qdef;

    double BSIM4qchqs;
    double BSIM4taunet;
    double BSIM4gtau;
    double BSIM4gtg;
    double BSIM4gtd;
    double BSIM4gts;
    double BSIM4gtb;
    double BSIM4SjctTempRevSatCur;
    double BSIM4DjctTempRevSatCur;
    double BSIM4SswTempRevSatCur;
    double BSIM4DswTempRevSatCur;
    double BSIM4SswgTempRevSatCur;
    double BSIM4DswgTempRevSatCur;

    // struct bsim4SizeDependParam  *pParam;
    std::shared_ptr<bsim4SizeDependParam> pParam = nullptr;

    bool BSIM4lGiven = false;
    bool BSIM4wGiven = false;
    bool BSIM4nfGiven = false;
    bool BSIM4minGiven = false;
    bool BSIM4drainAreaGiven = false;
    bool BSIM4sourceAreaGiven    = false;
    bool BSIM4drainSquaresGiven  = false;
    bool BSIM4sourceSquaresGiven = false;
    bool BSIM4drainPerimeterGiven    = false;
    bool BSIM4sourcePerimeterGiven   = false;
    bool BSIM4saGiven = false;
    bool BSIM4sbGiven = false;
    bool BSIM4sdGiven = false;
    bool BSIM4scaGiven = false;
    bool BSIM4scbGiven = false;
    bool BSIM4sccGiven = false;
    bool BSIM4scGiven = false;
    bool BSIM4rbdbGiven   = false;
    bool BSIM4rbsbGiven   = false;
    bool BSIM4rbpbGiven   = false;
    bool BSIM4rbpdGiven   = false;
    bool BSIM4rbpsGiven   = false;
    bool BSIM4delvtoGiven   = false;
    bool BSIM4xgwGiven   = false;
    bool BSIM4ngconGiven   = false;
    bool BSIM4icVDSGiven = false;
    bool BSIM4icVGSGiven = false;
    bool BSIM4icVBSGiven = false;
    bool BSIM4trnqsModGiven = false;
    bool BSIM4acnqsModGiven = false;
    bool BSIM4rbodyModGiven = false;
    bool BSIM4rgateModGiven = false;
    bool BSIM4geoModGiven = false;
    bool BSIM4rgeoModGiven = false;  
};

// BSIM4states index
    constexpr int BSIM4vbd = 0;
    constexpr int BSIM4vbs = 1;
    constexpr int BSIM4vgs = 2;
    constexpr int BSIM4vds = 3;
    constexpr int BSIM4vdbs = 4;
    constexpr int BSIM4vdbd = 5;
    constexpr int BSIM4vsbs = 6;
    constexpr int BSIM4vges = 7;
    constexpr int BSIM4vgms = 8;
    constexpr int BSIM4vses = 9;
    constexpr int BSIM4vdes = 10;

    constexpr int BSIM4qb = 11;
    constexpr int BSIM4cqb = 12;
    constexpr int BSIM4qg = 13;
    constexpr int BSIM4cqg = 14;
    constexpr int BSIM4qd = 15;
    constexpr int BSIM4cqd = 16;
    constexpr int BSIM4qgmid = 17;
    constexpr int BSIM4cqgmid = 18;

    constexpr int BSIM4qbs  = 19;
    constexpr int BSIM4cqbs  = 20;
    constexpr int BSIM4qbd  = 21;
    constexpr int BSIM4cqbd  = 22;

    constexpr int BSIM4qcheq = 23;
    constexpr int BSIM4cqcheq = 24;
    constexpr int BSIM4qcdump = 25;
    constexpr int BSIM4cqcdump = 26;
    constexpr int BSIM4qdef = 27;
    constexpr int BSIM4qs = 28;

    constexpr int BSIM4numStates = 29;
// end of BSIM4states index

/* Instance parameters */
constexpr int BSIM4_NMOS =  1;
constexpr int BSIM4_PMOS = -1;

constexpr int BSIM4_W                   = 1;
constexpr int BSIM4_L                   = 2;
constexpr int BSIM4_AS                  = 3;
constexpr int BSIM4_AD                  = 4;
constexpr int BSIM4_PS                  = 5;
constexpr int BSIM4_PD                  = 6;
constexpr int BSIM4_NRS                 = 7;
constexpr int BSIM4_NRD                 = 8;
constexpr int BSIM4_OFF                 = 9;
constexpr int BSIM4_IC                  = 10;
constexpr int BSIM4_IC_VDS              = 11;
constexpr int BSIM4_IC_VGS              = 12;
constexpr int BSIM4_IC_VBS              = 13;
constexpr int BSIM4_TRNQSMOD            = 14;
constexpr int BSIM4_RBODYMOD            = 15;
constexpr int BSIM4_RGATEMOD            = 16;
constexpr int BSIM4_GEOMOD              = 17;
constexpr int BSIM4_RGEOMOD             = 18;
constexpr int BSIM4_NF                  = 19;
constexpr int BSIM4_MIN                 = 20;
constexpr int BSIM4_ACNQSMOD            = 22;
constexpr int BSIM4_RBDB                = 23;
constexpr int BSIM4_RBSB                = 24;
constexpr int BSIM4_RBPB                = 25;
constexpr int BSIM4_RBPS                = 26;
constexpr int BSIM4_RBPD                = 27;
constexpr int BSIM4_SA                  = 28;
constexpr int BSIM4_SB                  = 29;
constexpr int BSIM4_SD                  = 30;
constexpr int BSIM4_DELVTO              = 31;
constexpr int BSIM4_XGW                 = 32;
constexpr int BSIM4_NGCON               = 33;
constexpr int BSIM4_SCA                 = 34;
constexpr int BSIM4_SCB                 = 35;
constexpr int BSIM4_SCC                 = 36;
constexpr int BSIM4_SC                  = 37;

/* Global parameters */
constexpr int BSIM4_MOD_TEMPEOT         = 66;
constexpr int BSIM4_MOD_LEFFEOT         = 67;
constexpr int BSIM4_MOD_WEFFEOT         = 68;
constexpr int BSIM4_MOD_UCSTE           = 69;
constexpr int BSIM4_MOD_LUCSTE          = 70;
constexpr int BSIM4_MOD_WUCSTE          = 71;
constexpr int BSIM4_MOD_PUCSTE          = 72;
constexpr int BSIM4_MOD_UCS             = 73;
constexpr int BSIM4_MOD_LUCS            = 74;
constexpr int BSIM4_MOD_WUCS            = 75;
constexpr int BSIM4_MOD_PUCS            = 76;
constexpr int BSIM4_MOD_CVCHARGEMOD     = 77;
constexpr int BSIM4_MOD_ADOS            = 78;
constexpr int BSIM4_MOD_BDOS            = 79;
constexpr int BSIM4_MOD_TEMPMOD         = 80;
constexpr int BSIM4_MOD_MTRLMOD         = 81;
constexpr int BSIM4_MOD_IGCMOD          = 82;
constexpr int BSIM4_MOD_IGBMOD          = 83;
constexpr int BSIM4_MOD_ACNQSMOD        = 84;
constexpr int BSIM4_MOD_FNOIMOD         = 85;
constexpr int BSIM4_MOD_RDSMOD          = 86;
constexpr int BSIM4_MOD_DIOMOD          = 87;
constexpr int BSIM4_MOD_PERMOD          = 88;
constexpr int BSIM4_MOD_GEOMOD          = 89;
constexpr int BSIM4_MOD_RGATEMOD        = 90;
constexpr int BSIM4_MOD_RBODYMOD        = 91;
constexpr int BSIM4_MOD_CAPMOD          = 92;
constexpr int BSIM4_MOD_TRNQSMOD        = 93;
constexpr int BSIM4_MOD_MOBMOD          = 94;
constexpr int BSIM4_MOD_TNOIMOD         = 95;
constexpr int BSIM4_MOD_EOT             = 96;
constexpr int BSIM4_MOD_VDDEOT          = 97;
constexpr int BSIM4_MOD_TOXE            = 98;
constexpr int BSIM4_MOD_CDSC            = 99;
constexpr int BSIM4_MOD_CDSCB           = 100;
constexpr int BSIM4_MOD_CIT             = 101;
constexpr int BSIM4_MOD_NFACTOR         = 102;
constexpr int BSIM4_MOD_XJ              = 103;
constexpr int BSIM4_MOD_VSAT            = 104;
constexpr int BSIM4_MOD_AT              = 105;
constexpr int BSIM4_MOD_A0              = 106;
constexpr int BSIM4_MOD_A1              = 107;
constexpr int BSIM4_MOD_A2              = 108;
constexpr int BSIM4_MOD_KETA            = 109;
constexpr int BSIM4_MOD_NSUB            = 110;
constexpr int BSIM4_MOD_PHIG            = 111;
constexpr int BSIM4_MOD_EPSRGATE        = 112;
constexpr int BSIM4_MOD_EASUB           = 113;
constexpr int BSIM4_MOD_EPSRSUB         = 114;
constexpr int BSIM4_MOD_NI0SUB          = 115;
constexpr int BSIM4_MOD_BG0SUB          = 116;
constexpr int BSIM4_MOD_TBGASUB         = 117;
constexpr int BSIM4_MOD_TBGBSUB         = 118;
constexpr int BSIM4_MOD_NDEP            = 119;
constexpr int BSIM4_MOD_NGATE           = 120;
constexpr int BSIM4_MOD_GAMMA1          = 121;
constexpr int BSIM4_MOD_GAMMA2          = 122;
constexpr int BSIM4_MOD_VBX             = 123;
constexpr int BSIM4_MOD_BINUNIT         = 124;
constexpr int BSIM4_MOD_VBM             = 125;
constexpr int BSIM4_MOD_XT              = 126;
constexpr int BSIM4_MOD_K1              = 129;
constexpr int BSIM4_MOD_KT1             = 130;
constexpr int BSIM4_MOD_KT1L            = 131;
constexpr int BSIM4_MOD_K2              = 132;
constexpr int BSIM4_MOD_KT2             = 133;
constexpr int BSIM4_MOD_K3              = 134;
constexpr int BSIM4_MOD_K3B             = 135;
constexpr int BSIM4_MOD_W0              = 136;
constexpr int BSIM4_MOD_LPE0            = 137;
constexpr int BSIM4_MOD_DVT0            = 138;
constexpr int BSIM4_MOD_DVT1            = 139;
constexpr int BSIM4_MOD_DVT2            = 140;
constexpr int BSIM4_MOD_DVT0W           = 141;
constexpr int BSIM4_MOD_DVT1W           = 142;
constexpr int BSIM4_MOD_DVT2W           = 143;
constexpr int BSIM4_MOD_DROUT           = 144;
constexpr int BSIM4_MOD_DSUB            = 145;
constexpr int BSIM4_MOD_VTH0            = 146;
constexpr int BSIM4_MOD_UA              = 147;
constexpr int BSIM4_MOD_UA1             = 148;
constexpr int BSIM4_MOD_UB              = 149;
constexpr int BSIM4_MOD_UB1             = 150;
constexpr int BSIM4_MOD_UC              = 151;
constexpr int BSIM4_MOD_UC1             = 152;
constexpr int BSIM4_MOD_U0              = 153;
constexpr int BSIM4_MOD_UTE             = 154;
constexpr int BSIM4_MOD_VOFF            = 155;
constexpr int BSIM4_MOD_DELTA           = 156;
constexpr int BSIM4_MOD_RDSW            = 157;
constexpr int BSIM4_MOD_PRT             = 158;
constexpr int BSIM4_MOD_LDD             = 159;
constexpr int BSIM4_MOD_ETA             = 160;
constexpr int BSIM4_MOD_ETA0            = 161;
constexpr int BSIM4_MOD_ETAB            = 162;
constexpr int BSIM4_MOD_PCLM            = 163;
constexpr int BSIM4_MOD_PDIBL1          = 164;
constexpr int BSIM4_MOD_PDIBL2          = 165;
constexpr int BSIM4_MOD_PSCBE1          = 166;
constexpr int BSIM4_MOD_PSCBE2          = 167;
constexpr int BSIM4_MOD_PVAG            = 168;
constexpr int BSIM4_MOD_WR              = 169;
constexpr int BSIM4_MOD_DWG             = 170;
constexpr int BSIM4_MOD_DWB             = 171;
constexpr int BSIM4_MOD_B0              = 172;
constexpr int BSIM4_MOD_B1              = 173;
constexpr int BSIM4_MOD_ALPHA0          = 174;
constexpr int BSIM4_MOD_BETA0           = 175;
constexpr int BSIM4_MOD_PDIBLB          = 178;
constexpr int BSIM4_MOD_PRWG            = 179;
constexpr int BSIM4_MOD_PRWB            = 180;
constexpr int BSIM4_MOD_CDSCD           = 181;
constexpr int BSIM4_MOD_AGS             = 182;
constexpr int BSIM4_MOD_FRINGE          = 184;
constexpr int BSIM4_MOD_CGSL            = 186;
constexpr int BSIM4_MOD_CGDL            = 187;
constexpr int BSIM4_MOD_CKAPPAS         = 188;
constexpr int BSIM4_MOD_CF              = 189;
constexpr int BSIM4_MOD_CLC             = 190;
constexpr int BSIM4_MOD_CLE             = 191;
constexpr int BSIM4_MOD_PARAMCHK        = 192;
constexpr int BSIM4_MOD_VERSION         = 193;
constexpr int BSIM4_MOD_VFBCV           = 194;
constexpr int BSIM4_MOD_ACDE            = 195;
constexpr int BSIM4_MOD_MOIN            = 196;
constexpr int BSIM4_MOD_NOFF            = 197;
constexpr int BSIM4_MOD_IJTHDFWD        = 198;
constexpr int BSIM4_MOD_ALPHA1          = 199;
constexpr int BSIM4_MOD_VFB             = 200;
constexpr int BSIM4_MOD_TOXM            = 201;
constexpr int BSIM4_MOD_TCJ             = 202;
constexpr int BSIM4_MOD_TCJSW           = 203;
constexpr int BSIM4_MOD_TCJSWG          = 204;
constexpr int BSIM4_MOD_TPB             = 205;
constexpr int BSIM4_MOD_TPBSW           = 206;
constexpr int BSIM4_MOD_TPBSWG          = 207;
constexpr int BSIM4_MOD_VOFFCV          = 208;
constexpr int BSIM4_MOD_GBMIN           = 209;
constexpr int BSIM4_MOD_RBDB            = 210;
constexpr int BSIM4_MOD_RBSB            = 211;
constexpr int BSIM4_MOD_RBPB            = 212;
constexpr int BSIM4_MOD_RBPS            = 213;
constexpr int BSIM4_MOD_RBPD            = 214;
constexpr int BSIM4_MOD_DMCG            = 215;
constexpr int BSIM4_MOD_DMCI            = 216;
constexpr int BSIM4_MOD_DMDG            = 217;
constexpr int BSIM4_MOD_XGW             = 218;
constexpr int BSIM4_MOD_XGL             = 219;
constexpr int BSIM4_MOD_RSHG            = 220;
constexpr int BSIM4_MOD_NGCON           = 221;
constexpr int BSIM4_MOD_AGIDL           = 222;
constexpr int BSIM4_MOD_BGIDL           = 223;
constexpr int BSIM4_MOD_EGIDL           = 224;
constexpr int BSIM4_MOD_IJTHSFWD        = 225;
constexpr int BSIM4_MOD_XJBVD           = 226;
constexpr int BSIM4_MOD_XJBVS           = 227;
constexpr int BSIM4_MOD_BVD             = 228;
constexpr int BSIM4_MOD_BVS             = 229;
constexpr int BSIM4_MOD_TOXP            = 230;
constexpr int BSIM4_MOD_DTOX            = 231;
constexpr int BSIM4_MOD_XRCRG1          = 232;
constexpr int BSIM4_MOD_XRCRG2          = 233;
constexpr int BSIM4_MOD_EU              = 234;
constexpr int BSIM4_MOD_IJTHSREV        = 235;
constexpr int BSIM4_MOD_IJTHDREV        = 236;
constexpr int BSIM4_MOD_MINV            = 237;
constexpr int BSIM4_MOD_VOFFL           = 238;
constexpr int BSIM4_MOD_PDITS           = 239;
constexpr int BSIM4_MOD_PDITSD          = 240;
constexpr int BSIM4_MOD_PDITSL          = 241;
constexpr int BSIM4_MOD_TNOIA           = 242;
constexpr int BSIM4_MOD_TNOIB           = 243;
constexpr int BSIM4_MOD_NTNOI           = 244;
constexpr int BSIM4_MOD_FPROUT          = 245;
constexpr int BSIM4_MOD_LPEB            = 246;
constexpr int BSIM4_MOD_DVTP0           = 247;
constexpr int BSIM4_MOD_DVTP1           = 248;
constexpr int BSIM4_MOD_CGIDL           = 249;
constexpr int BSIM4_MOD_PHIN            = 250;
constexpr int BSIM4_MOD_RDSWMIN         = 251;
constexpr int BSIM4_MOD_RSW             = 252;
constexpr int BSIM4_MOD_RDW             = 253;
constexpr int BSIM4_MOD_RDWMIN          = 254;
constexpr int BSIM4_MOD_RSWMIN          = 255;
constexpr int BSIM4_MOD_NSD             = 256;
constexpr int BSIM4_MOD_CKAPPAD         = 257;
constexpr int BSIM4_MOD_DMCGT           = 258;
constexpr int BSIM4_MOD_AIGC            = 259;
constexpr int BSIM4_MOD_BIGC            = 260;
constexpr int BSIM4_MOD_CIGC            = 261;
constexpr int BSIM4_MOD_AIGBACC         = 262;
constexpr int BSIM4_MOD_BIGBACC         = 263;
constexpr int BSIM4_MOD_CIGBACC         = 264;
constexpr int BSIM4_MOD_AIGBINV         = 265;
constexpr int BSIM4_MOD_BIGBINV         = 266;
constexpr int BSIM4_MOD_CIGBINV         = 267;
constexpr int BSIM4_MOD_NIGC            = 268;
constexpr int BSIM4_MOD_NIGBACC         = 269;
constexpr int BSIM4_MOD_NIGBINV         = 270;
constexpr int BSIM4_MOD_NTOX            = 271;
constexpr int BSIM4_MOD_TOXREF          = 272;
constexpr int BSIM4_MOD_EIGBINV         = 273;
constexpr int BSIM4_MOD_PIGCD           = 274;
constexpr int BSIM4_MOD_POXEDGE         = 275;
constexpr int BSIM4_MOD_EPSROX          = 276;
constexpr int BSIM4_MOD_AIGSD           = 277;
constexpr int BSIM4_MOD_BIGSD           = 278;
constexpr int BSIM4_MOD_CIGSD           = 279;
constexpr int BSIM4_MOD_JSWGS           = 280;
constexpr int BSIM4_MOD_JSWGD           = 281;
constexpr int BSIM4_MOD_LAMBDA          = 282;
constexpr int BSIM4_MOD_VTL             = 283;
constexpr int BSIM4_MOD_LC              = 284;
constexpr int BSIM4_MOD_XN              = 285;
constexpr int BSIM4_MOD_RNOIA           = 286;
constexpr int BSIM4_MOD_RNOIB           = 287;
constexpr int BSIM4_MOD_VFBSDOFF        = 288;
constexpr int BSIM4_MOD_LINTNOI         = 289;
constexpr int BSIM4_MOD_UD              = 290;
constexpr int BSIM4_MOD_UD1             = 291;
constexpr int BSIM4_MOD_UP              = 292;
constexpr int BSIM4_MOD_LP              = 293;
constexpr int BSIM4_MOD_TVOFF           = 294;
constexpr int BSIM4_MOD_TVFBSDOFF       = 295;
constexpr int BSIM4_MOD_MINVCV          = 296;
constexpr int BSIM4_MOD_VOFFCVL         = 297;
constexpr int BSIM4_MOD_MTRLCOMPATMOD   = 380;

/* Length dependence */
constexpr int BSIM4_MOD_LCDSC            = 301;
constexpr int BSIM4_MOD_LCDSCB           = 302;
constexpr int BSIM4_MOD_LCIT             = 303;
constexpr int BSIM4_MOD_LNFACTOR         = 304;
constexpr int BSIM4_MOD_LXJ              = 305;
constexpr int BSIM4_MOD_LVSAT            = 306;
constexpr int BSIM4_MOD_LAT              = 307;
constexpr int BSIM4_MOD_LA0              = 308;
constexpr int BSIM4_MOD_LA1              = 309;
constexpr int BSIM4_MOD_LA2              = 310;
constexpr int BSIM4_MOD_LKETA            = 311;
constexpr int BSIM4_MOD_LNSUB            = 312;
constexpr int BSIM4_MOD_LNDEP            = 313;
constexpr int BSIM4_MOD_LNGATE           = 315;
constexpr int BSIM4_MOD_LGAMMA1          = 316;
constexpr int BSIM4_MOD_LGAMMA2          = 317;
constexpr int BSIM4_MOD_LVBX             = 318;
constexpr int BSIM4_MOD_LVBM             = 320;
constexpr int BSIM4_MOD_LXT              = 322;
constexpr int BSIM4_MOD_LK1              = 325;
constexpr int BSIM4_MOD_LKT1             = 326;
constexpr int BSIM4_MOD_LKT1L            = 327;
constexpr int BSIM4_MOD_LK2              = 328;
constexpr int BSIM4_MOD_LKT2             = 329;
constexpr int BSIM4_MOD_LK3              = 330;
constexpr int BSIM4_MOD_LK3B             = 331;
constexpr int BSIM4_MOD_LW0              = 332;
constexpr int BSIM4_MOD_LLPE0            = 333;
constexpr int BSIM4_MOD_LDVT0            = 334;
constexpr int BSIM4_MOD_LDVT1            = 335;
constexpr int BSIM4_MOD_LDVT2            = 336;
constexpr int BSIM4_MOD_LDVT0W           = 337;
constexpr int BSIM4_MOD_LDVT1W           = 338;
constexpr int BSIM4_MOD_LDVT2W           = 339;
constexpr int BSIM4_MOD_LDROUT           = 340;
constexpr int BSIM4_MOD_LDSUB            = 341;
constexpr int BSIM4_MOD_LVTH0            = 342;
constexpr int BSIM4_MOD_LUA              = 343;
constexpr int BSIM4_MOD_LUA1             = 344;
constexpr int BSIM4_MOD_LUB              = 345;
constexpr int BSIM4_MOD_LUB1             = 346;
constexpr int BSIM4_MOD_LUC              = 347;
constexpr int BSIM4_MOD_LUC1             = 348;
constexpr int BSIM4_MOD_LU0              = 349;
constexpr int BSIM4_MOD_LUTE             = 350;
constexpr int BSIM4_MOD_LVOFF            = 351;
constexpr int BSIM4_MOD_LDELTA           = 352;
constexpr int BSIM4_MOD_LRDSW            = 353;
constexpr int BSIM4_MOD_LPRT             = 354;
constexpr int BSIM4_MOD_LLDD             = 355;
constexpr int BSIM4_MOD_LETA             = 356;
constexpr int BSIM4_MOD_LETA0            = 357;
constexpr int BSIM4_MOD_LETAB            = 358;
constexpr int BSIM4_MOD_LPCLM            = 359;
constexpr int BSIM4_MOD_LPDIBL1          = 360;
constexpr int BSIM4_MOD_LPDIBL2          = 361;
constexpr int BSIM4_MOD_LPSCBE1          = 362;
constexpr int BSIM4_MOD_LPSCBE2          = 363;
constexpr int BSIM4_MOD_LPVAG            = 364;
constexpr int BSIM4_MOD_LWR              = 365;
constexpr int BSIM4_MOD_LDWG             = 366;
constexpr int BSIM4_MOD_LDWB             = 367;
constexpr int BSIM4_MOD_LB0              = 368;
constexpr int BSIM4_MOD_LB1              = 369;
constexpr int BSIM4_MOD_LALPHA0          = 370;
constexpr int BSIM4_MOD_LBETA0           = 371;
constexpr int BSIM4_MOD_LPDIBLB          = 374;
constexpr int BSIM4_MOD_LPRWG            = 375;
constexpr int BSIM4_MOD_LPRWB            = 376;
constexpr int BSIM4_MOD_LCDSCD           = 377;
constexpr int BSIM4_MOD_LAGS             = 378;

constexpr int BSIM4_MOD_LFRINGE          = 381;
constexpr int BSIM4_MOD_LCGSL            = 383;
constexpr int BSIM4_MOD_LCGDL            = 384;
constexpr int BSIM4_MOD_LCKAPPAS         = 385;
constexpr int BSIM4_MOD_LCF              = 386;
constexpr int BSIM4_MOD_LCLC             = 387;
constexpr int BSIM4_MOD_LCLE             = 388;
constexpr int BSIM4_MOD_LVFBCV           = 389;
constexpr int BSIM4_MOD_LACDE            = 390;
constexpr int BSIM4_MOD_LMOIN            = 391;
constexpr int BSIM4_MOD_LNOFF            = 392;
constexpr int BSIM4_MOD_LALPHA1          = 394;
constexpr int BSIM4_MOD_LVFB             = 395;
constexpr int BSIM4_MOD_LVOFFCV          = 396;
constexpr int BSIM4_MOD_LAGIDL           = 397;
constexpr int BSIM4_MOD_LBGIDL           = 398;
constexpr int BSIM4_MOD_LEGIDL           = 399;
constexpr int BSIM4_MOD_LXRCRG1          = 400;
constexpr int BSIM4_MOD_LXRCRG2          = 401;
constexpr int BSIM4_MOD_LEU              = 402;
constexpr int BSIM4_MOD_LMINV            = 403;
constexpr int BSIM4_MOD_LPDITS           = 404;
constexpr int BSIM4_MOD_LPDITSD          = 405;
constexpr int BSIM4_MOD_LFPROUT          = 406;
constexpr int BSIM4_MOD_LLPEB            = 407;
constexpr int BSIM4_MOD_LDVTP0           = 408;
constexpr int BSIM4_MOD_LDVTP1           = 409;
constexpr int BSIM4_MOD_LCGIDL           = 410;
constexpr int BSIM4_MOD_LPHIN            = 411;
constexpr int BSIM4_MOD_LRSW             = 412;
constexpr int BSIM4_MOD_LRDW             = 413;
constexpr int BSIM4_MOD_LNSD             = 414;
constexpr int BSIM4_MOD_LCKAPPAD         = 415;
constexpr int BSIM4_MOD_LAIGC            = 416;
constexpr int BSIM4_MOD_LBIGC            = 417;
constexpr int BSIM4_MOD_LCIGC            = 418;
constexpr int BSIM4_MOD_LAIGBACC         = 419;
constexpr int BSIM4_MOD_LBIGBACC         = 420;
constexpr int BSIM4_MOD_LCIGBACC         = 421;
constexpr int BSIM4_MOD_LAIGBINV         = 422;
constexpr int BSIM4_MOD_LBIGBINV         = 423;
constexpr int BSIM4_MOD_LCIGBINV         = 424;
constexpr int BSIM4_MOD_LNIGC            = 425;
constexpr int BSIM4_MOD_LNIGBACC         = 426;
constexpr int BSIM4_MOD_LNIGBINV         = 427;
constexpr int BSIM4_MOD_LNTOX            = 428;
constexpr int BSIM4_MOD_LEIGBINV         = 429;
constexpr int BSIM4_MOD_LPIGCD           = 430;
constexpr int BSIM4_MOD_LPOXEDGE         = 431;
constexpr int BSIM4_MOD_LAIGSD           = 432;
constexpr int BSIM4_MOD_LBIGSD           = 433;
constexpr int BSIM4_MOD_LCIGSD           = 434;

constexpr int BSIM4_MOD_LLAMBDA          = 435;
constexpr int BSIM4_MOD_LVTL             = 436;
constexpr int BSIM4_MOD_LXN              = 437;
constexpr int BSIM4_MOD_LVFBSDOFF        = 438;
constexpr int BSIM4_MOD_LUD              = 439;
constexpr int BSIM4_MOD_LUD1             = 440;
constexpr int BSIM4_MOD_LUP              = 441;
constexpr int BSIM4_MOD_LLP              = 442;
constexpr int BSIM4_MOD_LMINVCV          = 443;

constexpr int BSIM4_MOD_FGIDL            = 444;          /* v4.7 New GIDL/GISL*/
constexpr int BSIM4_MOD_KGIDL            = 445;         /* v4.7 New GIDL/GISL*/
constexpr int BSIM4_MOD_RGIDL            = 446;          /* v4.7 New GIDL/GISL*/
constexpr int BSIM4_MOD_FGISL            = 447;          /* v4.7 New GIDL/GISL*/
constexpr int BSIM4_MOD_KGISL            = 448;          /* v4.7 New GIDL/GISL*/
constexpr int BSIM4_MOD_RGISL            = 449;          /* v4.7 New GIDL/GISL*/
constexpr int BSIM4_MOD_LFGIDL           = 450;          /* v4.7 New GIDL/GISL*/
constexpr int BSIM4_MOD_LKGIDL           = 451;          /* v4.7 New GIDL/GISL*/
constexpr int BSIM4_MOD_LRGIDL           = 452;          /* v4.7 New GIDL/GISL*/
constexpr int BSIM4_MOD_LFGISL           = 453;          /* v4.7 New GIDL/GISL*/
constexpr int BSIM4_MOD_LKGISL           = 454;          /* v4.7 New GIDL/GISL*/
constexpr int BSIM4_MOD_LRGISL           = 455;          /* v4.7 New GIDL/GISL*/
constexpr int BSIM4_MOD_WFGIDL           = 456;          /* v4.7 New GIDL/GISL*/
constexpr int BSIM4_MOD_WKGIDL           = 457;          /* v4.7 New GIDL/GISL*/
constexpr int BSIM4_MOD_WRGIDL           = 458;         /* v4.7 New GIDL/GISL*/
constexpr int BSIM4_MOD_WFGISL           = 459;          /* v4.7 New GIDL/GISL*/
constexpr int BSIM4_MOD_WKGISL           = 460;          /* v4.7 New GIDL/GISL*/
constexpr int BSIM4_MOD_WRGISL           = 461;          /* v4.7 New GIDL/GISL*/
constexpr int BSIM4_MOD_PFGIDL           = 462;          /* v4.7 New GIDL/GISL*/
constexpr int BSIM4_MOD_PKGIDL           = 463;          /* v4.7 New GIDL/GISL*/
constexpr int BSIM4_MOD_PRGIDL           = 464;          /* v4.7 New GIDL/GISL*/
constexpr int BSIM4_MOD_PFGISL           = 465;          /* v4.7 New GIDL/GISL*/
constexpr int BSIM4_MOD_PKGISL           = 466;          /* v4.7 New GIDL/GISL*/
constexpr int BSIM4_MOD_PRGISL           = 467;          /* v4.7 New GIDL/GISL*/
constexpr int BSIM4_MOD_GIDLMOD          = 379;          /* v4.7 New GIDL/GISL*/
constexpr int BSIM4_MOD_DVTP2           = 468;           /* v4.7 NEW DIBL/Rout*/
constexpr int BSIM4_MOD_DVTP3           = 469;           /* v4.7 NEW DIBL/Rout*/
constexpr int BSIM4_MOD_DVTP4           = 470;           /* v4.7 NEW DIBL/Rout*/
constexpr int BSIM4_MOD_DVTP5           = 471;           /* v4.7 NEW DIBL/Rout*/
constexpr int BSIM4_MOD_LDVTP2          = 472;           /* v4.7 NEW DIBL/Rout*/
constexpr int BSIM4_MOD_LDVTP3          = 473;           /* v4.7 NEW DIBL/Rout*/
constexpr int BSIM4_MOD_LDVTP4          = 474;           /* v4.7 NEW DIBL/Rout*/
constexpr int BSIM4_MOD_LDVTP5          = 475;           /* v4.7 NEW DIBL/Rout*/
constexpr int BSIM4_MOD_WDVTP2          = 476;           /* v4.7 NEW DIBL/Rout*/
constexpr int BSIM4_MOD_WDVTP3          = 477;           /* v4.7 NEW DIBL/Rout*/
constexpr int BSIM4_MOD_WDVTP4          = 478;           /* v4.7 NEW DIBL/Rout*/
constexpr int BSIM4_MOD_WDVTP5          = 479;           /* v4.7 NEW DIBL/Rout*/
constexpr int BSIM4_MOD_PDVTP2          = 480;           /* v4.7 NEW DIBL/Rout*/
constexpr int BSIM4_MOD_PDVTP3          = 298;           /* v4.7 NEW DIBL/Rout*/
constexpr int BSIM4_MOD_PDVTP4          = 299;           /* v4.7 NEW DIBL/Rout*/
constexpr int BSIM4_MOD_PDVTP5          = 300;           /* v4.7 NEW DIBL/Rout*/

/* Width dependence */
constexpr int BSIM4_MOD_WCDSC            = 481;
constexpr int BSIM4_MOD_WCDSCB           = 482;
constexpr int BSIM4_MOD_WCIT             = 483;
constexpr int BSIM4_MOD_WNFACTOR         = 484;
constexpr int BSIM4_MOD_WXJ              = 485;
constexpr int BSIM4_MOD_WVSAT            = 486;
constexpr int BSIM4_MOD_WAT              = 487;
constexpr int BSIM4_MOD_WA0              = 488;
constexpr int BSIM4_MOD_WA1              = 489;
constexpr int BSIM4_MOD_WA2              = 490;
constexpr int BSIM4_MOD_WKETA            = 491;
constexpr int BSIM4_MOD_WNSUB            = 492;
constexpr int BSIM4_MOD_WNDEP            = 493;
constexpr int BSIM4_MOD_WNGATE           = 495;
constexpr int BSIM4_MOD_WGAMMA1          = 496;
constexpr int BSIM4_MOD_WGAMMA2          = 497;
constexpr int BSIM4_MOD_WVBX             = 498;
constexpr int BSIM4_MOD_WVBM             = 500;
constexpr int BSIM4_MOD_WXT              = 502;
constexpr int BSIM4_MOD_WK1              = 505;
constexpr int BSIM4_MOD_WKT1             = 506;
constexpr int BSIM4_MOD_WKT1L            = 507;
constexpr int BSIM4_MOD_WK2              = 508;
constexpr int BSIM4_MOD_WKT2             = 509;
constexpr int BSIM4_MOD_WK3              = 510;
constexpr int BSIM4_MOD_WK3B             = 511;
constexpr int BSIM4_MOD_WW0              = 512;
constexpr int BSIM4_MOD_WLPE0            = 513;
constexpr int BSIM4_MOD_WDVT0            = 514;
constexpr int BSIM4_MOD_WDVT1            = 515;
constexpr int BSIM4_MOD_WDVT2            = 516;
constexpr int BSIM4_MOD_WDVT0W           = 517;
constexpr int BSIM4_MOD_WDVT1W           = 518;
constexpr int BSIM4_MOD_WDVT2W           = 519;
constexpr int BSIM4_MOD_WDROUT           = 520;
constexpr int BSIM4_MOD_WDSUB            = 521;
constexpr int BSIM4_MOD_WVTH0            = 522;
constexpr int BSIM4_MOD_WUA              = 523;
constexpr int BSIM4_MOD_WUA1             = 524;
constexpr int BSIM4_MOD_WUB              = 525;
constexpr int BSIM4_MOD_WUB1             = 526;
constexpr int BSIM4_MOD_WUC              = 527;
constexpr int BSIM4_MOD_WUC1             = 528;
constexpr int BSIM4_MOD_WU0              = 529;
constexpr int BSIM4_MOD_WUTE             = 530;
constexpr int BSIM4_MOD_WVOFF            = 531;
constexpr int BSIM4_MOD_WDELTA           = 532;
constexpr int BSIM4_MOD_WRDSW            = 533;
constexpr int BSIM4_MOD_WPRT             = 534;
constexpr int BSIM4_MOD_WLDD             = 535;
constexpr int BSIM4_MOD_WETA             = 536;
constexpr int BSIM4_MOD_WETA0            = 537;
constexpr int BSIM4_MOD_WETAB            = 538;
constexpr int BSIM4_MOD_WPCLM            = 539;
constexpr int BSIM4_MOD_WPDIBL1          = 540;
constexpr int BSIM4_MOD_WPDIBL2          = 541;
constexpr int BSIM4_MOD_WPSCBE1          = 542;
constexpr int BSIM4_MOD_WPSCBE2          = 543;
constexpr int BSIM4_MOD_WPVAG            = 544;
constexpr int BSIM4_MOD_WWR              = 545;
constexpr int BSIM4_MOD_WDWG             = 546;
constexpr int BSIM4_MOD_WDWB             = 547;
constexpr int BSIM4_MOD_WB0              = 548;
constexpr int BSIM4_MOD_WB1              = 549;
constexpr int BSIM4_MOD_WALPHA0          = 550;
constexpr int BSIM4_MOD_WBETA0           = 551;
constexpr int BSIM4_MOD_WPDIBLB          = 554;
constexpr int BSIM4_MOD_WPRWG            = 555;
constexpr int BSIM4_MOD_WPRWB            = 556;
constexpr int BSIM4_MOD_WCDSCD           = 557;
constexpr int BSIM4_MOD_WAGS             = 558;

constexpr int BSIM4_MOD_WFRINGE          = 561;
constexpr int BSIM4_MOD_WCGSL            = 563;
constexpr int BSIM4_MOD_WCGDL            = 564;
constexpr int BSIM4_MOD_WCKAPPAS         = 565;
constexpr int BSIM4_MOD_WCF              = 566;
constexpr int BSIM4_MOD_WCLC             = 567;
constexpr int BSIM4_MOD_WCLE             = 568;
constexpr int BSIM4_MOD_WVFBCV           = 569;
constexpr int BSIM4_MOD_WACDE            = 570;
constexpr int BSIM4_MOD_WMOIN            = 571;
constexpr int BSIM4_MOD_WNOFF            = 572;
constexpr int BSIM4_MOD_WALPHA1          = 574;
constexpr int BSIM4_MOD_WVFB             = 575;
constexpr int BSIM4_MOD_WVOFFCV          = 576;
constexpr int BSIM4_MOD_WAGIDL           = 577;
constexpr int BSIM4_MOD_WBGIDL           = 578;
constexpr int BSIM4_MOD_WEGIDL           = 579;
constexpr int BSIM4_MOD_WXRCRG1          = 580;
constexpr int BSIM4_MOD_WXRCRG2          = 581;
constexpr int BSIM4_MOD_WEU              = 582;
constexpr int BSIM4_MOD_WMINV            = 583;
constexpr int BSIM4_MOD_WPDITS           = 584;
constexpr int BSIM4_MOD_WPDITSD          = 585;
constexpr int BSIM4_MOD_WFPROUT          = 586;
constexpr int BSIM4_MOD_WLPEB            = 587;
constexpr int BSIM4_MOD_WDVTP0           = 588;
constexpr int BSIM4_MOD_WDVTP1           = 589;
constexpr int BSIM4_MOD_WCGIDL           = 590;
constexpr int BSIM4_MOD_WPHIN            = 591;
constexpr int BSIM4_MOD_WRSW             = 592;
constexpr int BSIM4_MOD_WRDW             = 593;
constexpr int BSIM4_MOD_WNSD             = 594;
constexpr int BSIM4_MOD_WCKAPPAD         = 595;
constexpr int BSIM4_MOD_WAIGC            = 596;
constexpr int BSIM4_MOD_WBIGC            = 597;
constexpr int BSIM4_MOD_WCIGC            = 598;
constexpr int BSIM4_MOD_WAIGBACC         = 599;
constexpr int BSIM4_MOD_WBIGBACC         = 600;
constexpr int BSIM4_MOD_WCIGBACC         = 601;
constexpr int BSIM4_MOD_WAIGBINV         = 602;
constexpr int BSIM4_MOD_WBIGBINV         = 603;
constexpr int BSIM4_MOD_WCIGBINV         = 604;
constexpr int BSIM4_MOD_WNIGC            = 605;
constexpr int BSIM4_MOD_WNIGBACC         = 606;
constexpr int BSIM4_MOD_WNIGBINV         = 607;
constexpr int BSIM4_MOD_WNTOX            = 608;
constexpr int BSIM4_MOD_WEIGBINV         = 609;
constexpr int BSIM4_MOD_WPIGCD           = 610;
constexpr int BSIM4_MOD_WPOXEDGE         = 611;
constexpr int BSIM4_MOD_WAIGSD           = 612;
constexpr int BSIM4_MOD_WBIGSD           = 613;
constexpr int BSIM4_MOD_WCIGSD           = 614;
constexpr int BSIM4_MOD_WLAMBDA          = 615;
constexpr int BSIM4_MOD_WVTL             = 616;
constexpr int BSIM4_MOD_WXN              = 617;
constexpr int BSIM4_MOD_WVFBSDOFF        = 618;
constexpr int BSIM4_MOD_WUD              = 619;
constexpr int BSIM4_MOD_WUD1             = 620;
constexpr int BSIM4_MOD_WUP              = 621;
constexpr int BSIM4_MOD_WLP              = 622;
constexpr int BSIM4_MOD_WMINVCV          = 623;

/* Cross-term dependence */
constexpr int BSIM4_MOD_PCDSC            = 661;
constexpr int BSIM4_MOD_PCDSCB           = 662;
constexpr int BSIM4_MOD_PCIT             = 663;
constexpr int BSIM4_MOD_PNFACTOR         = 664;
constexpr int BSIM4_MOD_PXJ              = 665;
constexpr int BSIM4_MOD_PVSAT            = 666;
constexpr int BSIM4_MOD_PAT              = 667;
constexpr int BSIM4_MOD_PA0              = 668;
constexpr int BSIM4_MOD_PA1              = 669;
constexpr int BSIM4_MOD_PA2              = 670;
constexpr int BSIM4_MOD_PKETA            = 671;
constexpr int BSIM4_MOD_PNSUB            = 672;
constexpr int BSIM4_MOD_PNDEP            = 673;
constexpr int BSIM4_MOD_PNGATE           = 675;
constexpr int BSIM4_MOD_PGAMMA1          = 676;
constexpr int BSIM4_MOD_PGAMMA2          = 677;
constexpr int BSIM4_MOD_PVBX             = 678;

constexpr int BSIM4_MOD_PVBM             = 680;

constexpr int BSIM4_MOD_PXT              = 682;
constexpr int BSIM4_MOD_PK1              = 685;
constexpr int BSIM4_MOD_PKT1             = 686;
constexpr int BSIM4_MOD_PKT1L            = 687;
constexpr int BSIM4_MOD_PK2              = 688;
constexpr int BSIM4_MOD_PKT2             = 689;
constexpr int BSIM4_MOD_PK3              = 690;
constexpr int BSIM4_MOD_PK3B             = 691;
constexpr int BSIM4_MOD_PW0              = 692;
constexpr int BSIM4_MOD_PLPE0            = 693;

constexpr int BSIM4_MOD_PDVT0            = 694;
constexpr int BSIM4_MOD_PDVT1            = 695;
constexpr int BSIM4_MOD_PDVT2            = 696;

constexpr int BSIM4_MOD_PDVT0W           = 697;
constexpr int BSIM4_MOD_PDVT1W           = 698;
constexpr int BSIM4_MOD_PDVT2W           = 699;

constexpr int BSIM4_MOD_PDROUT           = 700;
constexpr int BSIM4_MOD_PDSUB            = 701;
constexpr int BSIM4_MOD_PVTH0            = 702;
constexpr int BSIM4_MOD_PUA              = 703;
constexpr int BSIM4_MOD_PUA1             = 704;
constexpr int BSIM4_MOD_PUB              = 705;
constexpr int BSIM4_MOD_PUB1             = 706;
constexpr int BSIM4_MOD_PUC              = 707;
constexpr int BSIM4_MOD_PUC1             = 708;
constexpr int BSIM4_MOD_PU0              = 709;
constexpr int BSIM4_MOD_PUTE             = 710;
constexpr int BSIM4_MOD_PVOFF            = 711;
constexpr int BSIM4_MOD_PDELTA           = 712;
constexpr int BSIM4_MOD_PRDSW            = 713;
constexpr int BSIM4_MOD_PPRT             = 714;
constexpr int BSIM4_MOD_PLDD             = 715;
constexpr int BSIM4_MOD_PETA             = 716;
constexpr int BSIM4_MOD_PETA0            = 717;
constexpr int BSIM4_MOD_PETAB            = 718;
constexpr int BSIM4_MOD_PPCLM            = 719;
constexpr int BSIM4_MOD_PPDIBL1          = 720;
constexpr int BSIM4_MOD_PPDIBL2          = 721;
constexpr int BSIM4_MOD_PPSCBE1          = 722;
constexpr int BSIM4_MOD_PPSCBE2          = 723;
constexpr int BSIM4_MOD_PPVAG            = 724;
constexpr int BSIM4_MOD_PWR              = 725;
constexpr int BSIM4_MOD_PDWG             = 726;
constexpr int BSIM4_MOD_PDWB             = 727;
constexpr int BSIM4_MOD_PB0              = 728;
constexpr int BSIM4_MOD_PB1              = 729;
constexpr int BSIM4_MOD_PALPHA0          = 730;
constexpr int BSIM4_MOD_PBETA0           = 731;
constexpr int BSIM4_MOD_PPDIBLB          = 734;

constexpr int BSIM4_MOD_PPRWG            = 735;
constexpr int BSIM4_MOD_PPRWB            = 736;

constexpr int BSIM4_MOD_PCDSCD           = 737;
constexpr int BSIM4_MOD_PAGS             = 738;

constexpr int BSIM4_MOD_PFRINGE          = 741;
constexpr int BSIM4_MOD_PCGSL            = 743;
constexpr int BSIM4_MOD_PCGDL            = 744;
constexpr int BSIM4_MOD_PCKAPPAS         = 745;
constexpr int BSIM4_MOD_PCF              = 746;
constexpr int BSIM4_MOD_PCLC             = 747;
constexpr int BSIM4_MOD_PCLE             = 748;
constexpr int BSIM4_MOD_PVFBCV           = 749;
constexpr int BSIM4_MOD_PACDE            = 750;
constexpr int BSIM4_MOD_PMOIN            = 751;
constexpr int BSIM4_MOD_PNOFF            = 752;
constexpr int BSIM4_MOD_PALPHA1          = 754;
constexpr int BSIM4_MOD_PVFB             = 755;
constexpr int BSIM4_MOD_PVOFFCV          = 756;
constexpr int BSIM4_MOD_PAGIDL           = 757;
constexpr int BSIM4_MOD_PBGIDL           = 758;
constexpr int BSIM4_MOD_PEGIDL           = 759;
constexpr int BSIM4_MOD_PXRCRG1          = 760;
constexpr int BSIM4_MOD_PXRCRG2          = 761;
constexpr int BSIM4_MOD_PEU              = 762;
constexpr int BSIM4_MOD_PMINV            = 763;
constexpr int BSIM4_MOD_PPDITS           = 764;
constexpr int BSIM4_MOD_PPDITSD          = 765;
constexpr int BSIM4_MOD_PFPROUT          = 766;
constexpr int BSIM4_MOD_PLPEB            = 767;
constexpr int BSIM4_MOD_PDVTP0           = 768;
constexpr int BSIM4_MOD_PDVTP1           = 769;
constexpr int BSIM4_MOD_PCGIDL           = 770;
constexpr int BSIM4_MOD_PPHIN            = 771;
constexpr int BSIM4_MOD_PRSW             = 772;
constexpr int BSIM4_MOD_PRDW             = 773;
constexpr int BSIM4_MOD_PNSD             = 774;
constexpr int BSIM4_MOD_PCKAPPAD         = 775;
constexpr int BSIM4_MOD_PAIGC            = 776;
constexpr int BSIM4_MOD_PBIGC            = 777;
constexpr int BSIM4_MOD_PCIGC            = 778;
constexpr int BSIM4_MOD_PAIGBACC         = 779;
constexpr int BSIM4_MOD_PBIGBACC         = 780;
constexpr int BSIM4_MOD_PCIGBACC         = 781;
constexpr int BSIM4_MOD_PAIGBINV         = 782;
constexpr int BSIM4_MOD_PBIGBINV         = 783;
constexpr int BSIM4_MOD_PCIGBINV         = 784;
constexpr int BSIM4_MOD_PNIGC            = 785;
constexpr int BSIM4_MOD_PNIGBACC         = 786;
constexpr int BSIM4_MOD_PNIGBINV         = 787;
constexpr int BSIM4_MOD_PNTOX            = 788;
constexpr int BSIM4_MOD_PEIGBINV         = 789;
constexpr int BSIM4_MOD_PPIGCD           = 790;
constexpr int BSIM4_MOD_PPOXEDGE         = 791;
constexpr int BSIM4_MOD_PAIGSD           = 792;
constexpr int BSIM4_MOD_PBIGSD           = 793;
constexpr int BSIM4_MOD_PCIGSD           = 794;

constexpr int BSIM4_MOD_SAREF            = 795;
constexpr int BSIM4_MOD_SBREF            = 796;
constexpr int BSIM4_MOD_KU0              = 797;
constexpr int BSIM4_MOD_KVSAT            = 798;
constexpr int BSIM4_MOD_TKU0             = 799;
constexpr int BSIM4_MOD_LLODKU0          = 800;
constexpr int BSIM4_MOD_WLODKU0          = 801;
constexpr int BSIM4_MOD_LLODVTH          = 802;
constexpr int BSIM4_MOD_WLODVTH          = 803;
constexpr int BSIM4_MOD_LKU0             = 804;
constexpr int BSIM4_MOD_WKU0             = 805;
constexpr int BSIM4_MOD_PKU0             = 806;
constexpr int BSIM4_MOD_KVTH0            = 807;
constexpr int BSIM4_MOD_LKVTH0           = 808;
constexpr int BSIM4_MOD_WKVTH0           = 809;
constexpr int BSIM4_MOD_PKVTH0           = 810;
constexpr int BSIM4_MOD_WLOD         = 811;
constexpr int BSIM4_MOD_STK2         = 812;
constexpr int BSIM4_MOD_LODK2        = 813;
constexpr int BSIM4_MOD_STETA0       = 814;
constexpr int BSIM4_MOD_LODETA0      = 815;

constexpr int BSIM4_MOD_WEB              = 816;
constexpr int BSIM4_MOD_WEC              = 817;
constexpr int BSIM4_MOD_KVTH0WE          = 818;
constexpr int BSIM4_MOD_K2WE             = 819;
constexpr int BSIM4_MOD_KU0WE            = 820;
constexpr int BSIM4_MOD_SCREF            = 821;
constexpr int BSIM4_MOD_WPEMOD           = 822;
constexpr int BSIM4_MOD_PMINVCV          = 823;

constexpr int BSIM4_MOD_PLAMBDA          = 825;
constexpr int BSIM4_MOD_PVTL             = 826;
constexpr int BSIM4_MOD_PXN              = 827;
constexpr int BSIM4_MOD_PVFBSDOFF        = 828;

constexpr int BSIM4_MOD_TNOM             = 831;
constexpr int BSIM4_MOD_CGSO             = 832;
constexpr int BSIM4_MOD_CGDO             = 833;
constexpr int BSIM4_MOD_CGBO             = 834;
constexpr int BSIM4_MOD_XPART            = 835;
constexpr int BSIM4_MOD_RSH              = 836;
constexpr int BSIM4_MOD_JSS              = 837;
constexpr int BSIM4_MOD_PBS              = 838;
constexpr int BSIM4_MOD_MJS              = 839;
constexpr int BSIM4_MOD_PBSWS            = 840;
constexpr int BSIM4_MOD_MJSWS            = 841;
constexpr int BSIM4_MOD_CJS              = 842;
constexpr int BSIM4_MOD_CJSWS            = 843;
constexpr int BSIM4_MOD_NMOS             = 844;
constexpr int BSIM4_MOD_PMOS             = 845;
constexpr int BSIM4_MOD_NOIA             = 846;
constexpr int BSIM4_MOD_NOIB             = 847;
constexpr int BSIM4_MOD_NOIC             = 848;
constexpr int BSIM4_MOD_LINT             = 849;
constexpr int BSIM4_MOD_LL               = 850;
constexpr int BSIM4_MOD_LLN              = 851;
constexpr int BSIM4_MOD_LW               = 852;
constexpr int BSIM4_MOD_LWN              = 853;
constexpr int BSIM4_MOD_LWL              = 854;
constexpr int BSIM4_MOD_LMIN             = 855;
constexpr int BSIM4_MOD_LMAX             = 856;
constexpr int BSIM4_MOD_WINT             = 857;
constexpr int BSIM4_MOD_WL               = 858;
constexpr int BSIM4_MOD_WLN              = 859;
constexpr int BSIM4_MOD_WW               = 860;
constexpr int BSIM4_MOD_WWN              = 861;
constexpr int BSIM4_MOD_WWL              = 862;
constexpr int BSIM4_MOD_WMIN             = 863;
constexpr int BSIM4_MOD_WMAX             = 864;
constexpr int BSIM4_MOD_DWC              = 865;
constexpr int BSIM4_MOD_DLC              = 866;
constexpr int BSIM4_MOD_XL               = 867;
constexpr int BSIM4_MOD_XW               = 868;
constexpr int BSIM4_MOD_EM               = 869;
constexpr int BSIM4_MOD_EF               = 870;
constexpr int BSIM4_MOD_AF               = 871;
constexpr int BSIM4_MOD_KF               = 872;
constexpr int BSIM4_MOD_NJS              = 873;
constexpr int BSIM4_MOD_XTIS             = 874;
constexpr int BSIM4_MOD_PBSWGS           = 875;
constexpr int BSIM4_MOD_MJSWGS           = 876;
constexpr int BSIM4_MOD_CJSWGS           = 877;
constexpr int BSIM4_MOD_JSWS             = 878;
constexpr int BSIM4_MOD_LLC              = 879;
constexpr int BSIM4_MOD_LWC              = 880;
constexpr int BSIM4_MOD_LWLC             = 881;
constexpr int BSIM4_MOD_WLC              = 882;
constexpr int BSIM4_MOD_WWC              = 883;
constexpr int BSIM4_MOD_WWLC             = 884;
constexpr int BSIM4_MOD_DWJ              = 885;
constexpr int BSIM4_MOD_JSD              = 886;
constexpr int BSIM4_MOD_PBD              = 887;
constexpr int BSIM4_MOD_MJD              = 888;
constexpr int BSIM4_MOD_PBSWD            = 889;
constexpr int BSIM4_MOD_MJSWD            = 890;
constexpr int BSIM4_MOD_CJD              = 891;
constexpr int BSIM4_MOD_CJSWD            = 892;
constexpr int BSIM4_MOD_NJD              = 893;
constexpr int BSIM4_MOD_XTID             = 894;
constexpr int BSIM4_MOD_PBSWGD           = 895;
constexpr int BSIM4_MOD_MJSWGD           = 896;
constexpr int BSIM4_MOD_CJSWGD           = 897;
constexpr int BSIM4_MOD_JSWD             = 898;
constexpr int BSIM4_MOD_DLCIG            = 899;

/* trap-assisted tunneling */

constexpr int BSIM4_MOD_JTSS             = 900;
constexpr int BSIM4_MOD_JTSD         = 901;
constexpr int BSIM4_MOD_JTSSWS       = 902;
constexpr int BSIM4_MOD_JTSSWD       = 903;
constexpr int BSIM4_MOD_JTSSWGS      = 904;
constexpr int BSIM4_MOD_JTSSWGD      = 905;
constexpr int BSIM4_MOD_NJTS         = 906;
constexpr int BSIM4_MOD_NJTSSW       = 907;
constexpr int BSIM4_MOD_NJTSSWG      = 908;
constexpr int BSIM4_MOD_XTSS         = 909;
constexpr int BSIM4_MOD_XTSD         = 910;
constexpr int BSIM4_MOD_XTSSWS       = 911;
constexpr int BSIM4_MOD_XTSSWD       = 912;
constexpr int BSIM4_MOD_XTSSWGS      = 913;
constexpr int BSIM4_MOD_XTSSWGD      = 914;
constexpr int BSIM4_MOD_TNJTS        = 915;
constexpr int BSIM4_MOD_TNJTSSW      = 916;
constexpr int BSIM4_MOD_TNJTSSWG     = 917;
constexpr int BSIM4_MOD_VTSS             = 918;
constexpr int BSIM4_MOD_VTSD         = 919;
constexpr int BSIM4_MOD_VTSSWS       = 920;
constexpr int BSIM4_MOD_VTSSWD       = 921;
constexpr int BSIM4_MOD_VTSSWGS      = 922;
constexpr int BSIM4_MOD_VTSSWGD      = 923;
constexpr int BSIM4_MOD_PUD              = 924;
constexpr int BSIM4_MOD_PUD1             = 925;
constexpr int BSIM4_MOD_PUP              = 926;
constexpr int BSIM4_MOD_PLP              = 927;
constexpr int BSIM4_MOD_JTWEFF           = 928;

/* device questions */
constexpr int BSIM4_DNODE                = 945;
constexpr int BSIM4_GNODEEXT             = 946;
constexpr int BSIM4_SNODE                = 947;
constexpr int BSIM4_BNODE                = 948;
constexpr int BSIM4_DNODEPRIME           = 949;
constexpr int BSIM4_GNODEPRIME           = 950;
constexpr int BSIM4_GNODEMIDE            = 951;
constexpr int BSIM4_GNODEMID             = 952;
constexpr int BSIM4_SNODEPRIME           = 953;
constexpr int BSIM4_BNODEPRIME           = 954;
constexpr int BSIM4_DBNODE               = 955;
constexpr int BSIM4_SBNODE               = 956;
constexpr int BSIM4_VBD                  = 957;
constexpr int BSIM4_VBS                  = 958;
constexpr int BSIM4_VGS                  = 959;
constexpr int BSIM4_VDS                  = 960;
constexpr int BSIM4_CD                   = 961;
constexpr int BSIM4_CBS                  = 962;
constexpr int BSIM4_CBD                  = 963;
constexpr int BSIM4_GM                   = 964;
constexpr int BSIM4_GDS                  = 965;
constexpr int BSIM4_GMBS                 = 966;
constexpr int BSIM4_GBD                  = 967;
constexpr int BSIM4_GBS                  = 968;
constexpr int BSIM4_QB                   = 969;
constexpr int BSIM4_CQB                  = 970;
constexpr int BSIM4_QG                   = 971;
constexpr int BSIM4_CQG                  = 972;
constexpr int BSIM4_QD                   = 973;
constexpr int BSIM4_CQD                  = 974;
constexpr int BSIM4_CGGB                 = 975;
constexpr int BSIM4_CGDB                 = 976;
constexpr int BSIM4_CGSB                 = 977;
constexpr int BSIM4_CBGB                 = 978;
constexpr int BSIM4_CAPBD                = 979;
constexpr int BSIM4_CQBD                 = 980;
constexpr int BSIM4_CAPBS                = 981;
constexpr int BSIM4_CQBS                 = 982;
constexpr int BSIM4_CDGB                 = 983;
constexpr int BSIM4_CDDB                 = 984;
constexpr int BSIM4_CDSB                 = 985;
constexpr int BSIM4_VON                  = 986;
constexpr int BSIM4_VDSAT                = 987;
constexpr int BSIM4_QBS                  = 988;
constexpr int BSIM4_QBD                  = 989;
constexpr int BSIM4_SOURCECONDUCT        = 990;
constexpr int BSIM4_DRAINCONDUCT         = 991;
constexpr int BSIM4_CBDB                 = 992;
constexpr int BSIM4_CBSB                 = 993;
constexpr int BSIM4_CSUB         = 994;
constexpr int BSIM4_QINV         = 995;
constexpr int BSIM4_IGIDL        = 996;
constexpr int BSIM4_CSGB                 = 997;
constexpr int BSIM4_CSDB                 = 998;
constexpr int BSIM4_CSSB                 = 999;
constexpr int BSIM4_CGBB                 = 1000;
constexpr int BSIM4_CDBB                 = 1001;
constexpr int BSIM4_CSBB                 = 1002;
constexpr int BSIM4_CBBB                 = 1003;
constexpr int BSIM4_QS                   = 1004;
constexpr int BSIM4_IGISL        = 1005;
constexpr int BSIM4_IGS          = 1006;
constexpr int BSIM4_IGD          = 1007;
constexpr int BSIM4_IGB          = 1008;
constexpr int BSIM4_IGCS         = 1009;
constexpr int BSIM4_IGCD         = 1010;
constexpr int BSIM4_QDEF         = 1011;
constexpr int BSIM4_DELVT0           = 1012;
constexpr int BSIM4_GCRG                 = 1013;
constexpr int BSIM4_GTAU                 = 1014;

constexpr int BSIM4_MOD_LTVOFF           = 1051;
constexpr int BSIM4_MOD_LTVFBSDOFF       = 1052;
constexpr int BSIM4_MOD_WTVOFF           = 1053;
constexpr int BSIM4_MOD_WTVFBSDOFF       = 1054;
constexpr int BSIM4_MOD_PTVOFF           = 1055;
constexpr int BSIM4_MOD_PTVFBSDOFF       = 1056;

constexpr int BSIM4_MOD_LKVTH0WE          = 1061;
constexpr int BSIM4_MOD_LK2WE             = 1062;
constexpr int BSIM4_MOD_LKU0WE        = 1063;
constexpr int BSIM4_MOD_WKVTH0WE          = 1064;
constexpr int BSIM4_MOD_WK2WE             = 1065;
constexpr int BSIM4_MOD_WKU0WE        = 1066;
constexpr int BSIM4_MOD_PKVTH0WE          = 1067;
constexpr int BSIM4_MOD_PK2WE             = 1068;
constexpr int BSIM4_MOD_PKU0WE        = 1069;

constexpr int BSIM4_MOD_RBPS0               = 1101;
constexpr int BSIM4_MOD_RBPSL               = 1102;
constexpr int BSIM4_MOD_RBPSW               = 1103;
constexpr int BSIM4_MOD_RBPSNF              = 1104;
constexpr int BSIM4_MOD_RBPD0               = 1105;
constexpr int BSIM4_MOD_RBPDL               = 1106;
constexpr int BSIM4_MOD_RBPDW               = 1107;
constexpr int BSIM4_MOD_RBPDNF              = 1108;

constexpr int BSIM4_MOD_RBPBX0              = 1109;
constexpr int BSIM4_MOD_RBPBXL              = 1110;
constexpr int BSIM4_MOD_RBPBXW              = 1111;
constexpr int BSIM4_MOD_RBPBXNF             = 1112;
constexpr int BSIM4_MOD_RBPBY0              = 1113;
constexpr int BSIM4_MOD_RBPBYL              = 1114;
constexpr int BSIM4_MOD_RBPBYW              = 1115;
constexpr int BSIM4_MOD_RBPBYNF             = 1116;

constexpr int BSIM4_MOD_RBSBX0              = 1117;
constexpr int BSIM4_MOD_RBSBY0              = 1118;
constexpr int BSIM4_MOD_RBDBX0              = 1119;
constexpr int BSIM4_MOD_RBDBY0              = 1120;

constexpr int BSIM4_MOD_RBSDBXL             = 1121;
constexpr int BSIM4_MOD_RBSDBXW             = 1122;
constexpr int BSIM4_MOD_RBSDBXNF            = 1123;
constexpr int BSIM4_MOD_RBSDBYL             = 1124;
constexpr int BSIM4_MOD_RBSDBYW             = 1125;
constexpr int BSIM4_MOD_RBSDBYNF            = 1126;

constexpr int BSIM4_MOD_AGISL               = 1200;
constexpr int BSIM4_MOD_BGISL               = 1201;
constexpr int BSIM4_MOD_EGISL               = 1202;
constexpr int BSIM4_MOD_CGISL               = 1203;
constexpr int BSIM4_MOD_LAGISL              = 1204;
constexpr int BSIM4_MOD_LBGISL              = 1205;
constexpr int BSIM4_MOD_LEGISL              = 1206;
constexpr int BSIM4_MOD_LCGISL              = 1207;
constexpr int BSIM4_MOD_WAGISL              = 1208;
constexpr int BSIM4_MOD_WBGISL              = 1209;
constexpr int BSIM4_MOD_WEGISL              = 1210;
constexpr int BSIM4_MOD_WCGISL              = 1211;
constexpr int BSIM4_MOD_PAGISL              = 1212;
constexpr int BSIM4_MOD_PBGISL              = 1213;
constexpr int BSIM4_MOD_PEGISL              = 1214;
constexpr int BSIM4_MOD_PCGISL              = 1215;

constexpr int BSIM4_MOD_AIGS                = 1220;
constexpr int BSIM4_MOD_BIGS                = 1221;
constexpr int BSIM4_MOD_CIGS                = 1222;
constexpr int BSIM4_MOD_LAIGS               = 1223;
constexpr int BSIM4_MOD_LBIGS               = 1224;
constexpr int BSIM4_MOD_LCIGS               = 1225;
constexpr int BSIM4_MOD_WAIGS               = 1226;
constexpr int BSIM4_MOD_WBIGS               = 1227;
constexpr int BSIM4_MOD_WCIGS               = 1228;
constexpr int BSIM4_MOD_PAIGS               = 1229;
constexpr int BSIM4_MOD_PBIGS               = 1230;
constexpr int BSIM4_MOD_PCIGS               = 1231;
constexpr int BSIM4_MOD_AIGD                = 1232;
constexpr int BSIM4_MOD_BIGD                = 1233;
constexpr int BSIM4_MOD_CIGD                = 1234;
constexpr int BSIM4_MOD_LAIGD               = 1235;
constexpr int BSIM4_MOD_LBIGD               = 1236;
constexpr int BSIM4_MOD_LCIGD               = 1237;
constexpr int BSIM4_MOD_WAIGD               = 1238;
constexpr int BSIM4_MOD_WBIGD               = 1239;
constexpr int BSIM4_MOD_WCIGD               = 1240;
constexpr int BSIM4_MOD_PAIGD               = 1241;
constexpr int BSIM4_MOD_PBIGD               = 1242;
constexpr int BSIM4_MOD_PCIGD               = 1243;
constexpr int BSIM4_MOD_DLCIGD              = 1244;

constexpr int BSIM4_MOD_NJTSD               = 1250;
constexpr int BSIM4_MOD_NJTSSWD             = 1251;
constexpr int BSIM4_MOD_NJTSSWGD            = 1252;
constexpr int BSIM4_MOD_TNJTSD              = 1253;
constexpr int BSIM4_MOD_TNJTSSWD            = 1254;
constexpr int BSIM4_MOD_TNJTSSWGD           = 1255;

/* v4.7 temp dep of leakage current  */

constexpr int BSIM4_MOD_TNFACTOR            = 1256;
constexpr int BSIM4_MOD_TETA0               = 1257;
constexpr int BSIM4_MOD_TVOFFCV             = 1258;
constexpr int BSIM4_MOD_LTNFACTOR           = 1260;
constexpr int BSIM4_MOD_LTETA0              = 1261;
constexpr int BSIM4_MOD_LTVOFFCV            = 1262;
constexpr int BSIM4_MOD_WTNFACTOR           = 1264;
constexpr int BSIM4_MOD_WTETA0              = 1265;
constexpr int BSIM4_MOD_WTVOFFCV            = 1266;
constexpr int BSIM4_MOD_PTNFACTOR           = 1268;
constexpr int BSIM4_MOD_PTETA0              = 1269;
constexpr int BSIM4_MOD_PTVOFFCV            = 1270;

/* tnoiMod=2 (v4.7) */
constexpr int BSIM4_MOD_TNOIC               = 1272;
constexpr int BSIM4_MOD_RNOIC               = 1273;
/* smoothing for gidl clamp (C.K.Dabhi) */
constexpr int BSIM4_MOD_GIDLCLAMP           = 1274;
/* Tuning for noise parameter BSIM4IdovVds (C.K.Dabhi) - request cadence */
constexpr int BSIM4_MOD_IDOVVDSC            = 1275;
} // namespace BSIM4
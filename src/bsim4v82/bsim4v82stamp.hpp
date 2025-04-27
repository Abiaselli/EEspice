/*
   b4ld.c - BSIM4v4.8.2
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
#include "bsim4v82.hpp"
#include "sim_variables.hpp"
#include "SPICEcompatible.hpp"
#include "devsup.hpp"
#include "bsim4v82const.hpp"

#include <cmath>
#include <armadillo>

namespace bsim4{

inline void dexp(const double A, double& B, double& C) noexcept {
    if (A > EXP_THRESHOLD) {
        B = MAX_EXP * (1.0 + (A) - EXP_THRESHOLD);
        C = MAX_EXP;
    } else if (A < -EXP_THRESHOLD) {
        B = MIN_EXP;
        C = 0;
    } else {
        B = std::exp(A);
        C = B;
    }
}

/* function to compute poly depletion effect */
int BSIM4polyDepletion(
    double  phi,
    double  ngate,
    double  epsgate,
    double  coxe,
    double  Vgs,
    double *Vgs_eff,
    double *dVgs_eff_dVg)
{
    double T1, T2, T3, T4, T5, T6, T7, T8;

    /* Poly Gate Si Depletion Effect */
    if ((ngate > 1.0e18) &&
        (ngate < 1.0e25) && (Vgs > phi) && (epsgate!=0)
       ){
        T1 = 1.0e6 * CHARGE * epsgate * ngate / (coxe * coxe);
        T8 = Vgs - phi;
        T4 = std::sqrt(1.0 + 2.0 * T8 / T1);
        T2 = 2.0 * T8 / (T4 + 1.0);
        T3 = 0.5 * T2 * T2 / T1; /* T3 = Vpoly */
        T7 = 1.12 - T3 - 0.05;
        T6 = std::sqrt(T7 * T7 + 0.224);
        T5 = 1.12 - 0.5 * (T7 + T6);
        *Vgs_eff = Vgs - T5;
        *dVgs_eff_dVg = 1.0 - (0.5 - 0.5 / T4) * (1.0 + T7 / T6);
    }
    else {
        *Vgs_eff = Vgs;
        *dVgs_eff_dVg = 1.0;
    }
    return(0);
}


int
BSIM4load(const BSIM4model &model, BSIM4V82 &instance, const SPICECompatible &spice, const arma::vec &presolution,
    const double CKTtemp, const double CKTgmin, arma::mat &LHS, arma::vec &RHS)
{
    double ceqgstot, dgstot_dvd, dgstot_dvg, dgstot_dvs, dgstot_dvb;
    double ceqgdtot, dgdtot_dvd, dgdtot_dvg, dgdtot_dvs, dgdtot_dvb;
    double gstot, gstotd, gstotg, gstots, gstotb, gspr, Rs, Rd;
    double gdtot, gdtotd, gdtotg, gdtots, gdtotb, gdpr;
    double vgs_eff, vgd_eff, dvgs_eff_dvg, dvgd_eff_dvg;
    double dRs_dvg, dRd_dvg, dRs_dvb, dRd_dvb;
    double dT0_dvg, dT1_dvb, dT3_dvg, dT3_dvb;
    double vses, vdes, vdedo, delvses, delvded, delvdes;
    double Isestot, cseshat, Idedtot, cdedhat;
    // #ifndef NEWCONV
    double tol0, tol1, tol2, tol3, tol4, tol5, tol6;
    
    double geltd, gcrg, gcrgg, gcrgd, gcrgs, gcrgb, ceqgcrg;
    double vges, vgms, vgedo, vgmdo, vged, vgmd, delvged, delvgmd;
    double delvges, delvgms, vgmb;
    double gcgmgmb=0.0, gcgmdb=0.0, gcgmsb=0.0, gcdgmb, gcsgmb;
    double gcgmbb=0.0, gcbgmb, qgmb, qgmid=0.0, ceqqgmid;
    
    double vbd, vbs, vds, vgb, vgd, vgs, vgdo;
    //  #ifndef PREDICTOR
    // double xfact;

    double vdbs, vdbd, vsbs, vsbdo, vsbd;
    double delvdbs, delvdbd, delvsbs;
    double delvbd_jct, delvbs_jct, vbs_jct, vbd_jct;
    
    double SourceSatCurrent, DrainSatCurrent;
    double ag0, qgb, von, cbhat, VgstNVt, ExpVgst;
    double ceqqb, ceqqd, ceqqg, ceqqjd=0.0, ceqqjs=0.0, ceq, geq;
    double cdrain, cdhat, ceqdrn, ceqbd, ceqbs, ceqjd, ceqjs, gjbd, gjbs;
    double czbd, czbdsw, czbdswg, czbs, czbssw, czbsswg, evbd, evbs, arg, sarg;
    double delvbd, delvbs, delvds, delvgd, delvgs;
    double Vfbeff, dVfbeff_dVg, dVfbeff_dVb, V3, V4;
    double gcbdb, gcbgb, gcbsb, gcddb, gcdgb, gcdsb, gcgdb, gcggb, gcgsb, gcsdb;
    double gcgbb, gcdbb, gcsbb, gcbbb;
    double gcdbdb, gcsbsb;
    double gcsgb, gcssb, MJD, MJSWD, MJSWGD, MJS, MJSWS, MJSWGS;
    double qgate=0.0, qbulk=0.0, qdrn=0.0, qsrc, cqgate, cqbody, cqdrn;
    double Vdb, Vds, Vgs, Vbs, Gmbs, FwdSum, RevSum;
    double Igidl, Ggidld, Ggidlg, Ggidlb;
    double Voxacc=0.0, dVoxacc_dVg=0.0, dVoxacc_dVb=0.0;
    double Voxdepinv=0.0, dVoxdepinv_dVg=0.0, dVoxdepinv_dVd=0.0, dVoxdepinv_dVb=0.0;
    double VxNVt=0.0, ExpVxNVt, Vaux=0.0, dVaux_dVg=0.0, dVaux_dVd=0.0, dVaux_dVb=0.0;
    double Igc, dIgc_dVg, dIgc_dVd, dIgc_dVb;
    double Igcs, dIgcs_dVg, dIgcs_dVd, dIgcs_dVb;
    double Igcd, dIgcd_dVg, dIgcd_dVd, dIgcd_dVb;
    double Igs, dIgs_dVg, dIgs_dVs, Igd, dIgd_dVg, dIgd_dVd;
    double Igbacc, dIgbacc_dVg, dIgbacc_dVb;
    double Igbinv, dIgbinv_dVg, dIgbinv_dVd, dIgbinv_dVb;
    double Pigcd, dPigcd_dVg, dPigcd_dVd, dPigcd_dVb;
    double Istoteq, gIstotg, gIstotd, gIstots, gIstotb;
    double Idtoteq, gIdtotg, gIdtotd, gIdtots, gIdtotb;
    double Ibtoteq, gIbtotg, gIbtotd, gIbtots, gIbtotb;
    double Igtoteq, gIgtotg, gIgtotd, gIgtots, gIgtotb;
    double Igstot, cgshat, Igdtot, cgdhat, Igbtot, cgbhat;
    double Vgs_eff, Vfb=0.0, Vth_NarrowW;
    /* double Vgd_eff, dVgd_eff_dVg;          v4.7.0 */
    double Phis, dPhis_dVb, sqrtPhis, dsqrtPhis_dVb, Vth, dVth_dVb, dVth_dVd;
    double Vgst, dVgst_dVg, dVgst_dVb, dVgs_eff_dVg, Nvtms, Nvtmd;
    double Vtm, Vtm0;
    double n, dn_dVb, dn_dVd, voffcv, noff, dnoff_dVd, dnoff_dVb;
    double V0, CoxWLcen, QovCox, LINK;
    double DeltaPhi, dDeltaPhi_dVg, VgDP, dVgDP_dVg;
    double Cox, Tox, Tcen, dTcen_dVg, dTcen_dVd, dTcen_dVb;
    double Ccen, Coxeff, dCoxeff_dVd, dCoxeff_dVg, dCoxeff_dVb;
    double Denomi, dDenomi_dVg, dDenomi_dVd, dDenomi_dVb;
    double ueff, dueff_dVg, dueff_dVd, dueff_dVb; 
    double Esat, Vdsat;
    double EsatL, dEsatL_dVg, dEsatL_dVd, dEsatL_dVb;
    double dVdsat_dVg, dVdsat_dVb, dVdsat_dVd, Vasat, dAlphaz_dVg, dAlphaz_dVb; 
    double dVasat_dVg, dVasat_dVb, dVasat_dVd, Va, dVa_dVd, dVa_dVg, dVa_dVb; 
    double Vbseff, dVbseff_dVb, VbseffCV, dVbseffCV_dVb; 
    double VgsteffVth, dT11_dVg;
    double Arg1, One_Third_CoxWL, Two_Third_CoxWL, Alphaz, CoxWL; 
    double T0=0.0, dT0_dVg, dT0_dVd, dT0_dVb;
    double T1, dT1_dVg, dT1_dVd, dT1_dVb;
    double T2, dT2_dVg, dT2_dVd, dT2_dVb;
    double T3, dT3_dVg, dT3_dVd, dT3_dVb;
    double T4, dT4_dVd, dT4_dVb;
    double T5, dT5_dVg, dT5_dVd, dT5_dVb;
    double T6, dT6_dVg, dT6_dVd, dT6_dVb;
    double T7, dT7_dVg, dT7_dVd, dT7_dVb;
    double T8, dT8_dVg, dT8_dVd, dT8_dVb;
    double T9, dT9_dVg, dT9_dVd, dT9_dVb;
    double T10, dT10_dVg, dT10_dVb, dT10_dVd; 
    double T11, T12, T13, T14;
    double tmp, Abulk, dAbulk_dVb, Abulk0, dAbulk0_dVb;
    double Cclm, dCclm_dVg, dCclm_dVd, dCclm_dVb;
    double FP, dFP_dVg, PvagTerm, dPvagTerm_dVg, dPvagTerm_dVd, dPvagTerm_dVb;
    double VADITS, dVADITS_dVg, dVADITS_dVd;
    double Lpe_Vb, dDITS_Sft_dVb, dDITS_Sft_dVd;
    double DITS_Sft2, dDITS_Sft2_dVd;        /* v4.7 New DITS */
    double VACLM, dVACLM_dVg, dVACLM_dVd, dVACLM_dVb;
    double VADIBL, dVADIBL_dVg, dVADIBL_dVd, dVADIBL_dVb;
    double Xdep, dXdep_dVb, lt1, dlt1_dVb, ltw, dltw_dVb, Delt_vth, dDelt_vth_dVb;
    double Theta0, dTheta0_dVb;
    double TempRatio, tmp1, tmp2, tmp3, tmp4;
    double DIBL_Sft, dDIBL_Sft_dVd, Lambda, dLambda_dVg;
    double Idtot, Ibtot, a1, ScalingFactor;
    
    double Vgsteff, dVgsteff_dVg, dVgsteff_dVd, dVgsteff_dVb; 
    double Vdseff, dVdseff_dVg, dVdseff_dVd, dVdseff_dVb; 
    double VdseffCV, dVdseffCV_dVg, dVdseffCV_dVd, dVdseffCV_dVb; 
    double diffVds, dAbulk_dVg;
    double beta, dbeta_dVg, dbeta_dVd, dbeta_dVb;
    double gche, dgche_dVg, dgche_dVd, dgche_dVb;
    double fgche1, dfgche1_dVg, dfgche1_dVd, dfgche1_dVb;
    double fgche2, dfgche2_dVg, dfgche2_dVd, dfgche2_dVb;
    double Idl, dIdl_dVg, dIdl_dVd, dIdl_dVb;
    double Idsa, dIdsa_dVg, dIdsa_dVd, dIdsa_dVb;
    double Ids, Gm, Gds, Gmb, devbs_dvb, devbd_dvb;
    double Isub, Gbd, Gbg, Gbb;
    double VASCBE, dVASCBE_dVg, dVASCBE_dVd, dVASCBE_dVb;
    double CoxeffWovL;
    double Rds, dRds_dVg, dRds_dVb, WVCox, WVCoxRds;
    double Vgst2Vtm, VdsatCV;
    double Leff, Weff, dWeff_dVg, dWeff_dVb;
    double AbulkCV, dAbulkCV_dVb;
    double qcheq, qdef, gqdef=0.0, cqdef=0.0, cqcheq=0.0;
    double gcqdb=0.0, gcqsb=0.0, gcqgb=0.0, gcqbb=0.0;
    double dxpart, sxpart, ggtg, ggtd, ggts, ggtb;
    double ddxpart_dVd, ddxpart_dVg, ddxpart_dVb, ddxpart_dVs;
    double dsxpart_dVd, dsxpart_dVg, dsxpart_dVb, dsxpart_dVs;
    double gbspsp, gbbdp, gbbsp, gbspg, gbspb, gbspdp; 
    double gbdpdp, gbdpg, gbdpb, gbdpsp; 
    double qgdo, qgso, cgdo, cgso;
    double Cgg, Cgd, Cgb, Cdg, Cdd, Cds;
    double Csg, Csd, Css, Csb, Cbg, Cbd, Cbb;
    double Cgg1, Cgd1, Cgb1, Cbg1, Cbb1, Cbd1, Qac0, Qsub0;
    double dQac0_dVg, dQac0_dVb, dQsub0_dVg, dQsub0_dVd, dQsub0_dVb;
    double ggidld, ggidlg, ggidlb, ggislg, ggislb, ggisls;
    double Igisl, Ggislg, Ggislb, Ggisls;
    double Nvtmrss, Nvtmrssws, Nvtmrsswgs;
    double Nvtmrsd, Nvtmrsswd, Nvtmrsswgd;
    
    double vs, Fsevl, dvs_dVg, dvs_dVd, dvs_dVb, dFsevl_dVg, dFsevl_dVd, dFsevl_dVb;
    double vgdx, vgsx, epssub, toxe, epsrox;
    struct bsim4SizeDependParam pParam;
    int ByPass, ChargeComputationNeeded, error, Check, Check1, Check2;

    // Node indexes in mna matrix
    auto get_NodeIndex = [](int node) -> int {
        return   node - 1;
    };
    const int dNode = get_NodeIndex(instance.BSIM4dNode);
    const int gNodeExt = get_NodeIndex(instance.BSIM4gNodeExt);
    const int sNode = get_NodeIndex(instance.BSIM4sNode);
    const int bNode = get_NodeIndex(instance.BSIM4bNode);
    const int dNodePrime = get_NodeIndex(instance.BSIM4dNodePrime);
    const int gNodePrime = get_NodeIndex(instance.BSIM4gNodePrime);
    const int gNodeMid = get_NodeIndex(instance.BSIM4gNodeMid);
    const int sNodePrime = get_NodeIndex(instance.BSIM4sNode);
    const int bNodePrime = get_NodeIndex(instance.BSIM4bNodePrime);
    const int dbNode = get_NodeIndex(instance.BSIM4dbNode);
    const int sbNode = get_NodeIndex(instance.BSIM4sbNode);
    const int qNode = get_NodeIndex(instance.BSIM4qNode);

    // Node voltages
    double volBSIM4dNode, volBSIM4dNodePrime;
    double volBSIM4gNodeExt, volBSIM4gNodePrime, volBSIM4gNodeMid;
    double volBSIM4sNode, volBSIM4sNodePrime;
    double volBSIM4bNode, volBSIM4bNodePrime;
    double volBSIM4dbNode, volBSIM4sbNode, volBSIM4qNode;

    // Imports all SPICEmode members into the current scope.
    using enum SPICECompatible::SPICEmode;
    auto CKTmode = spice.getMode();

    ScalingFactor = 1.0e-9;
    ChargeComputationNeeded =
                    ((CKTmode & (MODEAC | MODETRAN | MODEINITSMSIG)) ||
                    ((CKTmode & MODETRANOP) && (CKTmode & MODEUIC)))
                    ? 1 : 0;

    Check = Check1 = Check2 = 1;
    ByPass = 0;
    pParam = *instance.pParam;

    // Loading node voltages

    if ((CKTmode & MODEINITSMSIG))
    {   
        vds = instance.BSIM4states0[BSIM4vds];
        vgs = instance.BSIM4states0[BSIM4vgs];
        vbs = instance.BSIM4states0[BSIM4vbs];
        vges = instance.BSIM4states0[BSIM4vges];
        vgms = instance.BSIM4states0[BSIM4vgms];
        vdbs = instance.BSIM4states0[BSIM4vdbs];
        vsbs = instance.BSIM4states0[BSIM4vsbs];
        vses = instance.BSIM4states0[BSIM4vses];
        vdes = instance.BSIM4states0[BSIM4vdes];
        qdef = instance.BSIM4states0[BSIM4qdef];
    }
    else if ((CKTmode & MODEINITTRAN))
    {   
        vds = instance.BSIM4states1[BSIM4vds];
        vgs = instance.BSIM4states1[BSIM4vgs];
        vbs = instance.BSIM4states1[BSIM4vbs];
        vges = instance.BSIM4states1[BSIM4vges];
        vgms = instance.BSIM4states1[BSIM4vgms];
        vdbs = instance.BSIM4states1[BSIM4vdbs];
        vsbs = instance.BSIM4states1[BSIM4vsbs];
        vses = instance.BSIM4states1[BSIM4vses];
        vdes = instance.BSIM4states1[BSIM4vdes];
        qdef = instance.BSIM4states1[BSIM4qdef];
    }
    else if ((CKTmode & MODEINITJCT) && !instance.BSIM4off)
    {   
        vds = model.BSIM4type * instance.BSIM4icVDS;
        vgs = vges = vgms = model.BSIM4type * instance.BSIM4icVGS;
        vbs = vdbs = vsbs = model.BSIM4type * instance.BSIM4icVBS;
        if (vds > 0.0)
        {   vdes = vds + 0.01;
            vses = -0.01;
        }
        else if (vds < 0.0)
        {   vdes = vds - 0.01;
            vses = 0.01;
        }
        else
        {   vdes = vses = 0.0; }

        qdef = 0.0;

        if ((vds == 0.0) && (vgs == 0.0) && (vbs == 0.0) &&
            ((CKTmode & (MODETRAN | MODEAC|MODEDCOP |
            MODEDCTRANCURVE)) || (!(CKTmode & MODEUIC))))
        {   vds = 0.1;
            vdes = 0.11;
            vses = -0.01;
            vgs = vges = vgms = model.BSIM4type * instance.BSIM4vth0 + 0.1;
            vbs = vdbs = vsbs = 0.0;
        }
    }
    else if ((CKTmode & (MODEINITJCT | MODEINITFIX)) &&
            (instance.BSIM4off))
    {   vds = vgs = vbs = vges = vgms = 0.0;
        vdbs = vsbs = vdes = vses = qdef = 0.0;
    }
    else
    {   // Calculate node voltages
        auto get_NodeVoltage = [&presolution](int node) -> double {
            return (node > 0) ? presolution(node - 1) : 0.0;
        };
        
        volBSIM4dNode       = get_NodeVoltage(instance.BSIM4dNode);
        volBSIM4dNodePrime  = get_NodeVoltage(instance.BSIM4dNodePrime);
        volBSIM4gNodeExt    = get_NodeVoltage(instance.BSIM4gNodeExt);
        volBSIM4gNodePrime  = get_NodeVoltage(instance.BSIM4gNodePrime);
        volBSIM4gNodeMid    = get_NodeVoltage(instance.BSIM4gNodeMid);
        volBSIM4sNode       = get_NodeVoltage(instance.BSIM4sNode);
        volBSIM4sNodePrime  = get_NodeVoltage(instance.BSIM4sNodePrime);
        volBSIM4bNode       = get_NodeVoltage(instance.BSIM4bNode);
        volBSIM4bNodePrime  = get_NodeVoltage(instance.BSIM4bNodePrime);
        volBSIM4dbNode      = get_NodeVoltage(instance.BSIM4dbNode);
        volBSIM4sbNode      = get_NodeVoltage(instance.BSIM4sbNode);
        volBSIM4qNode       = get_NodeVoltage(instance.BSIM4qNode);


        vds = model.BSIM4type
            * (volBSIM4dNodePrime
            - volBSIM4sNodePrime);
        vgs = model.BSIM4type
            * (volBSIM4gNodePrime
            - volBSIM4sNodePrime);
        vbs = model.BSIM4type
            * (volBSIM4bNodePrime
            - volBSIM4sNodePrime);
        vges = model.BSIM4type
            * (volBSIM4gNodeExt
            - volBSIM4sNodePrime);
        vgms = model.BSIM4type
            * (volBSIM4gNodeMid
            - volBSIM4sNodePrime);
        vdbs = model.BSIM4type
            * (volBSIM4dbNode
            - volBSIM4sNodePrime);
        vsbs = model.BSIM4type
            * (volBSIM4sbNode
            - volBSIM4sNodePrime);
        vses = model.BSIM4type
            * (volBSIM4sNode
            - volBSIM4sNodePrime);
        vdes = model.BSIM4type
            * (volBSIM4dNode
            - volBSIM4sNodePrime);
        qdef = model.BSIM4type
            * (volBSIM4qNode);


        vgdo = instance.BSIM4states0[BSIM4vgs]
                - instance.BSIM4states0[BSIM4vds];
        vgedo = instance.BSIM4states0[BSIM4vges]
                - instance.BSIM4states0[BSIM4vds];
        vgmdo = instance.BSIM4states0[BSIM4vgms]
                - instance.BSIM4states0[BSIM4vds];
                
        vbd = vbs - vds;
        vdbd = vdbs - vds;
        vgd = vgs - vds;
        vged = vges - vds;
        vgmd = vgms - vds;

        delvbd = vbd - instance.BSIM4states0[BSIM4vbd];
        delvdbd = vdbd - instance.BSIM4states0[BSIM4vdbd];
        delvgd = vgd - vgdo;
        delvged = vged - vgedo;
        delvgmd = vgmd - vgmdo;

        delvds = vds - instance.BSIM4states0[BSIM4vds];
        delvgs = vgs - instance.BSIM4states0[BSIM4vgs];
        delvges = vges - instance.BSIM4states0[BSIM4vges];
        delvgms = vgms - instance.BSIM4states0[BSIM4vgms];
        delvbs = vbs - instance.BSIM4states0[BSIM4vbs];
        delvdbs = vdbs - instance.BSIM4states0[BSIM4vdbs];
        delvsbs = vsbs - instance.BSIM4states0[BSIM4vsbs];

        delvses = vses - instance.BSIM4states0[BSIM4vses];
        vdedo = instance.BSIM4states0[BSIM4vdes]
                - instance.BSIM4states0[BSIM4vds];
        delvdes = vdes - instance.BSIM4states0[BSIM4vdes];
        delvded = vdes - vds - vdedo;

        delvbd_jct = (!instance.BSIM4rbodyMod) ? delvbd : delvdbd;
        delvbs_jct = (!instance.BSIM4rbodyMod) ? delvbs : delvsbs;
        if (instance.BSIM4mode >= 0)
        {   Idtot = instance.BSIM4cd + instance.BSIM4csub - instance.BSIM4cbd
                + instance.BSIM4Igidl;
            cdhat = Idtot - instance.BSIM4gbd * delvbd_jct
                    + (instance.BSIM4gmbs + instance.BSIM4gbbs + instance.BSIM4ggidlb) * delvbs
                    + (instance.BSIM4gm + instance.BSIM4gbgs + instance.BSIM4ggidlg) * delvgs
                    + (instance.BSIM4gds + instance.BSIM4gbds + instance.BSIM4ggidld) * delvds;
            Ibtot = instance.BSIM4cbs + instance.BSIM4cbd
                - instance.BSIM4Igidl - instance.BSIM4Igisl - instance.BSIM4csub;
            cbhat = Ibtot + instance.BSIM4gbd * delvbd_jct
                    + instance.BSIM4gbs * delvbs_jct - (instance.BSIM4gbbs + instance.BSIM4ggidlb)
                    * delvbs - (instance.BSIM4gbgs + instance.BSIM4ggidlg) * delvgs
                    - (instance.BSIM4gbds + instance.BSIM4ggidld - instance.BSIM4ggisls) * delvds
                    - instance.BSIM4ggislg * delvgd - instance.BSIM4ggislb * delvbd;

           Igstot = instance.BSIM4Igs + instance.BSIM4Igcs;
           cgshat = Igstot + (instance.BSIM4gIgsg + instance.BSIM4gIgcsg) * delvgs
                    + instance.BSIM4gIgcsd * delvds + instance.BSIM4gIgcsb * delvbs;

           Igdtot = instance.BSIM4Igd + instance.BSIM4Igcd;
           cgdhat = Igdtot + instance.BSIM4gIgdg * delvgd + instance.BSIM4gIgcdg * delvgs
                          + instance.BSIM4gIgcdd * delvds + instance.BSIM4gIgcdb * delvbs;

           Igbtot = instance.BSIM4Igb;
           cgbhat = instance.BSIM4Igb + instance.BSIM4gIgbg * delvgs + instance.BSIM4gIgbd
              * delvds + instance.BSIM4gIgbb * delvbs;
        }
        else
        {   Idtot = instance.BSIM4cd + instance.BSIM4cbd - instance.BSIM4Igidl; /* bugfix */
            cdhat = Idtot + instance.BSIM4gbd * delvbd_jct + instance.BSIM4gmbs
                    * delvbd + instance.BSIM4gm * delvgd
                    - (instance.BSIM4gds + instance.BSIM4ggidls) * delvds
                    - instance.BSIM4ggidlg * delvgs - instance.BSIM4ggidlb * delvbs;
            Ibtot = instance.BSIM4cbs + instance.BSIM4cbd
                    - instance.BSIM4Igidl - instance.BSIM4Igisl - instance.BSIM4csub;
            cbhat = Ibtot + instance.BSIM4gbs * delvbs_jct + instance.BSIM4gbd
                    * delvbd_jct - (instance.BSIM4gbbs + instance.BSIM4ggislb) * delvbd
                    - (instance.BSIM4gbgs + instance.BSIM4ggislg) * delvgd
                    + (instance.BSIM4gbds + instance.BSIM4ggisld - instance.BSIM4ggidls) * delvds
                    - instance.BSIM4ggidlg * delvgs - instance.BSIM4ggidlb * delvbs;

            Igstot = instance.BSIM4Igs + instance.BSIM4Igcd;
            cgshat = Igstot + instance.BSIM4gIgsg * delvgs + instance.BSIM4gIgcdg * delvgd
                    - instance.BSIM4gIgcdd * delvds + instance.BSIM4gIgcdb * delvbd;

            Igdtot = instance.BSIM4Igd + instance.BSIM4Igcs;
            cgdhat = Igdtot + (instance.BSIM4gIgdg + instance.BSIM4gIgcsg) * delvgd
                    - instance.BSIM4gIgcsd * delvds + instance.BSIM4gIgcsb * delvbd;

            Igbtot = instance.BSIM4Igb;
            cgbhat = instance.BSIM4Igb + instance.BSIM4gIgbg * delvgd
                    - instance.BSIM4gIgbd * delvds + instance.BSIM4gIgbb * delvbd;
        }

        Isestot = instance.BSIM4gstot * instance.BSIM4states0[BSIM4vses];
        cseshat = Isestot + instance.BSIM4gstot * delvses
                + instance.BSIM4gstotd * delvds + instance.BSIM4gstotg * delvgs
                + instance.BSIM4gstotb * delvbs;

        Idedtot = instance.BSIM4gdtot * vdedo;
        cdedhat = Idedtot + instance.BSIM4gdtot * delvded
                + instance.BSIM4gdtotd * delvds + instance.BSIM4gdtotg * delvgs
                + instance.BSIM4gdtotb * delvbs;

        von = instance.BSIM4von;
        if (instance.BSIM4states0[BSIM4vds] >= 0.0)
        {   vgs = DEVfetlim(vgs, instance.BSIM4states0[BSIM4vgs], von);
            
            vds = vgs - vgd;
            vds = DEVlimvds(vds, instance.BSIM4states0[BSIM4vds]);
            vgd = vgs - vds;
            if (instance.BSIM4rgateMod == 3)
            {   vges = DEVfetlim(vges, instance.BSIM4states0[BSIM4vges], von);
                vgms = DEVfetlim(vgms, instance.BSIM4states0[BSIM4vgms], von);
                vged = vges - vds;
                vgmd = vgms - vds;
            }
            else if ((instance.BSIM4rgateMod == 1) || (instance.BSIM4rgateMod == 2))
            {   vges = DEVfetlim(vges, instance.BSIM4states0[BSIM4vges], von);
                vged = vges - vds;
            }

            if (model.BSIM4rdsMod)
            {   vdes = DEVlimvds(vdes, instance.BSIM4states0[BSIM4vdes]);
                vses = -DEVlimvds(-vses, -(instance.BSIM4states0[BSIM4vses]));
            }
        }
        else
        {    vgd = DEVfetlim(vgd, vgdo, von);
            vds = vgs - vgd;
            vds = -DEVlimvds(-vds, -(instance.BSIM4states0[BSIM4vds]));
            vgs = vgd + vds;

            if (instance.BSIM4rgateMod == 3)
            {   vged = DEVfetlim(vged, vgedo, von);
                vges = vged + vds;
                vgmd = DEVfetlim(vgmd, vgmdo, von);
                vgms = vgmd + vds;
            }
            if ((instance.BSIM4rgateMod == 1) || (instance.BSIM4rgateMod == 2))
            {   vged = DEVfetlim(vged, vgedo, von);
                vges = vged + vds;
            }
            if (model.BSIM4rdsMod)
            {   vdes = -DEVlimvds(-vdes, -(instance.BSIM4states0[BSIM4vdes]));
                vses = DEVlimvds(vses, instance.BSIM4states0[BSIM4vses]);
            }
        }

        if (vds >= 0.0)
        {   vbs = DEVpnjlim(vbs, instance.BSIM4states0[BSIM4vbs],
                                CONSTvt0, model.BSIM4vcrit, &Check);
            vbd = vbs - vds;
            if (instance.BSIM4rbodyMod)
            {   vdbs = DEVpnjlim(vdbs, instance.BSIM4states0[BSIM4vdbs],
                                CONSTvt0, model.BSIM4vcrit, &Check1);
                vdbd = vdbs - vds;
                vsbs = DEVpnjlim(vsbs, instance.BSIM4states0[BSIM4vsbs],
                                CONSTvt0, model.BSIM4vcrit, &Check2);
                if ((Check1 == 0) && (Check2 == 0))
                    Check = 0;
                else
                    Check = 1;
            }
        }
        else
        {   vbd = DEVpnjlim(vbd, instance.BSIM4states0[BSIM4vbd],
                            CONSTvt0, model.BSIM4vcrit, &Check);
            vbs = vbd + vds;
            if (instance.BSIM4rbodyMod)
            {   vdbd = DEVpnjlim(vdbd, instance.BSIM4states0[BSIM4vdbd],
                                CONSTvt0, model.BSIM4vcrit, &Check1);
                vdbs = vdbd + vds;
                vsbdo = instance.BSIM4states0[BSIM4vsbs]
                        - instance.BSIM4states0[BSIM4vds];
                vsbd = vsbs - vds;
                vsbd = DEVpnjlim(vsbd, vsbdo, CONSTvt0, model.BSIM4vcrit, &Check2);
                vsbs = vsbd + vds;
                if ((Check1 == 0) && (Check2 == 0))
                    Check = 0;
                else
                    Check = 1;
            }
        }
    }

    /* Calculate DC currents and their derivatives */
    vbd = vbs - vds;
    vgd = vgs - vds;
    vgb = vgs - vbs;
    vged = vges - vds;
    vgmd = vgms - vds;
    vgmb = vgms - vbs;
    vdbd = vdbs - vds;

    vbs_jct = (!instance.BSIM4rbodyMod) ? vbs : vsbs;
    vbd_jct = (!instance.BSIM4rbodyMod) ? vbd : vdbd;

    /* Source/drain junction diode DC model begins */
    Nvtms = model.BSIM4vtm * model.BSIM4SjctEmissionCoeff;
/*    if ((instance.BSIM4Aseff <= 0.0) && (instance.BSIM4Pseff <= 0.0))
      {   SourceSatCurrent = 1.0e-14;
      } v4.7 */
    if ((instance.BSIM4Aseff <= 0.0) && (instance.BSIM4Pseff <= 0.0))
    {   SourceSatCurrent = 0.0;
    }
    else
    {   SourceSatCurrent = instance.BSIM4Aseff * model.BSIM4SjctTempSatCurDensity
                        + instance.BSIM4Pseff * model.BSIM4SjctSidewallTempSatCurDensity
                        + pParam.BSIM4weffCJ * instance.BSIM4nf
                        * model.BSIM4SjctGateSidewallTempSatCurDensity;
    }

    if (SourceSatCurrent <= 0.0)
    {   instance.BSIM4gbs = CKTgmin;
        instance.BSIM4cbs = instance.BSIM4gbs * vbs_jct;
    }
    else
    {   switch(model.BSIM4dioMod)
        {   case 0:
                evbs = std::exp(vbs_jct / Nvtms);
                T1 = model.BSIM4xjbvs * std::exp(-(model.BSIM4bvs + vbs_jct) / Nvtms);
                /* WDLiu: Magic T1 in this form; different from BSIM4 beta. */
                instance.BSIM4gbs = SourceSatCurrent * (evbs + T1) / Nvtms + CKTgmin;
                instance.BSIM4cbs = SourceSatCurrent * (evbs + instance.BSIM4XExpBVS
                    - T1 - 1.0) + CKTgmin * vbs_jct;
                break;
            case 1:
                T2 = vbs_jct / Nvtms;
                if (T2 < -EXP_THRESHOLD)
                {   instance.BSIM4gbs = CKTgmin;
                    instance.BSIM4cbs = SourceSatCurrent * (MIN_EXP - 1.0) + CKTgmin * vbs_jct;
                }
                else if (vbs_jct <= instance.BSIM4vjsmFwd)
                {   evbs = std::exp(T2);
                    instance.BSIM4gbs = SourceSatCurrent * evbs / Nvtms + CKTgmin;
                    instance.BSIM4cbs = SourceSatCurrent * (evbs - 1.0) + CKTgmin * vbs_jct;
                }
                else
                {   T0 = instance.BSIM4IVjsmFwd / Nvtms;
                    instance.BSIM4gbs = T0 + CKTgmin;
                    instance.BSIM4cbs = instance.BSIM4IVjsmFwd - SourceSatCurrent + T0
                        * (vbs_jct - instance.BSIM4vjsmFwd) + CKTgmin * vbs_jct;
                }
                break;
            case 2:
                if (vbs_jct < instance.BSIM4vjsmRev)
                {   T0 = vbs_jct / Nvtms;
                    if (T0 < -EXP_THRESHOLD)
                    {   evbs = MIN_EXP;
                        devbs_dvb = 0.0;
                    }
                    else
                    {   evbs = std::exp(T0);
                        devbs_dvb = evbs / Nvtms;
                    }

                    T1 = evbs - 1.0;
                    T2 = instance.BSIM4IVjsmRev + instance.BSIM4SslpRev
                        * (vbs_jct - instance.BSIM4vjsmRev);
                    instance.BSIM4gbs = devbs_dvb * T2 + T1 * instance.BSIM4SslpRev + CKTgmin;
                    instance.BSIM4cbs = T1 * T2 + CKTgmin * vbs_jct;
                }
                else if (vbs_jct <= instance.BSIM4vjsmFwd)
                {   T0 = vbs_jct / Nvtms;
                    if (T0 < -EXP_THRESHOLD)
                    {   evbs = MIN_EXP;
                        devbs_dvb = 0.0;
                    }
                    else
                    {   evbs = std::exp(T0);
                        devbs_dvb = evbs / Nvtms;
                    }

                    T1 = (model.BSIM4bvs + vbs_jct) / Nvtms;
                    if (T1 > EXP_THRESHOLD)
                    {   T2 = MIN_EXP;
                        T3 = 0.0;
                    }
                    else
                    {   T2 = std::exp(-T1);
                        T3 = -T2 /Nvtms;
                    }
                    instance.BSIM4gbs = SourceSatCurrent * (devbs_dvb - model.BSIM4xjbvs * T3)
                        + CKTgmin;
                    instance.BSIM4cbs = SourceSatCurrent * (evbs + instance.BSIM4XExpBVS - 1.0
                        - model.BSIM4xjbvs * T2) + CKTgmin * vbs_jct;
                }
                else
                {     
                    instance.BSIM4gbs = instance.BSIM4SslpFwd + CKTgmin;
                    instance.BSIM4cbs = instance.BSIM4IVjsmFwd + instance.BSIM4SslpFwd * (vbs_jct
                        - instance.BSIM4vjsmFwd) + CKTgmin * vbs_jct;
                }
                break;
            default: break;
        }
    }

    Nvtmd = model.BSIM4vtm * model.BSIM4DjctEmissionCoeff;
/*    if ((instance.BSIM4Adeff <= 0.0) && (instance.BSIM4Pdeff <= 0.0))
      {   DrainSatCurrent = 1.0e-14;
      } v4.7 */
    if ((instance.BSIM4Adeff <= 0.0) && (instance.BSIM4Pdeff <= 0.0))
    {   DrainSatCurrent = 0.0;
    }
    else
    {   DrainSatCurrent = instance.BSIM4Adeff * model.BSIM4DjctTempSatCurDensity
                        + instance.BSIM4Pdeff * model.BSIM4DjctSidewallTempSatCurDensity
                        + pParam.BSIM4weffCJ * instance.BSIM4nf
                        * model.BSIM4DjctGateSidewallTempSatCurDensity;
    }

    if (DrainSatCurrent <= 0.0)
    {   instance.BSIM4gbd = CKTgmin;
        instance.BSIM4cbd = instance.BSIM4gbd * vbd_jct;
    }
    else
    {   
        switch(model.BSIM4dioMod)
        {   
            case 0:
                evbd = std::exp(vbd_jct / Nvtmd);
                T1 = model.BSIM4xjbvd * std::exp(-(model.BSIM4bvd + vbd_jct) / Nvtmd);
                /* WDLiu: Magic T1 in this form; different from BSIM4 beta. */
                instance.BSIM4gbd = DrainSatCurrent * (evbd + T1) / Nvtmd + CKTgmin;
                instance.BSIM4cbd = DrainSatCurrent * (evbd + instance.BSIM4XExpBVD
                                - T1 - 1.0) + CKTgmin * vbd_jct;
                break;
            case 1:
                T2 = vbd_jct / Nvtmd;
                if (T2 < -EXP_THRESHOLD)
                {   instance.BSIM4gbd = CKTgmin;
                    instance.BSIM4cbd = DrainSatCurrent * (MIN_EXP - 1.0) + CKTgmin * vbd_jct;
                }
                else if (vbd_jct <= instance.BSIM4vjdmFwd)
                {   evbd = std::exp(T2);
                    instance.BSIM4gbd = DrainSatCurrent * evbd / Nvtmd + CKTgmin;
                    instance.BSIM4cbd = DrainSatCurrent * (evbd - 1.0) + CKTgmin * vbd_jct;
                }
                else
                {   T0 = instance.BSIM4IVjdmFwd / Nvtmd;
                    instance.BSIM4gbd = T0 + CKTgmin;
                    instance.BSIM4cbd = instance.BSIM4IVjdmFwd - DrainSatCurrent + T0
                                    * (vbd_jct - instance.BSIM4vjdmFwd) + CKTgmin * vbd_jct;
                }
                break;
            case 2:
                if (vbd_jct < instance.BSIM4vjdmRev)
                {   T0 = vbd_jct / Nvtmd;
                    if (T0 < -EXP_THRESHOLD)
                    {   evbd = MIN_EXP;
                        devbd_dvb = 0.0;
                    }
                    else
                    {   evbd = std::exp(T0);
                        devbd_dvb = evbd / Nvtmd;
                    }

                    T1 = evbd - 1.0;
                    T2 = instance.BSIM4IVjdmRev + instance.BSIM4DslpRev
                        * (vbd_jct - instance.BSIM4vjdmRev);
                    instance.BSIM4gbd = devbd_dvb * T2 + T1 * instance.BSIM4DslpRev + CKTgmin;
                    instance.BSIM4cbd = T1 * T2 + CKTgmin * vbd_jct;
                }
                else if (vbd_jct <= instance.BSIM4vjdmFwd)
                {   T0 = vbd_jct / Nvtmd;
                    if (T0 < -EXP_THRESHOLD)
                    {   evbd = MIN_EXP;
                        devbd_dvb = 0.0;
                    }
                    else
                    {   evbd = std::exp(T0);
                        devbd_dvb = evbd / Nvtmd;
                    }

                    T1 = (model.BSIM4bvd + vbd_jct) / Nvtmd;
                    if (T1 > EXP_THRESHOLD)
                    {   T2 = MIN_EXP;
                        T3 = 0.0;
                    }
                    else
                    {   T2 = std::exp(-T1);
                        T3 = -T2 /Nvtmd;
                    }
                    instance.BSIM4gbd = DrainSatCurrent * (devbd_dvb - model.BSIM4xjbvd * T3)
                                    + CKTgmin;
                    instance.BSIM4cbd = DrainSatCurrent * (evbd + instance.BSIM4XExpBVD - 1.0
                                    - model.BSIM4xjbvd * T2) + CKTgmin * vbd_jct;
                }
                else
                {   instance.BSIM4gbd = instance.BSIM4DslpFwd + CKTgmin;
                    instance.BSIM4cbd = instance.BSIM4IVjdmFwd + instance.BSIM4DslpFwd * (vbd_jct
                                    - instance.BSIM4vjdmFwd) + CKTgmin * vbd_jct;
                }
                break;
            default: break;
        }
    }

    /* trap-assisted tunneling and recombination current for reverse bias  */
    Nvtmrssws = model.BSIM4vtm0 * model.BSIM4njtsswstemp;
    Nvtmrsswgs = model.BSIM4vtm0 * model.BSIM4njtsswgstemp;
    Nvtmrss = model.BSIM4vtm0 * model.BSIM4njtsstemp;
    Nvtmrsswd = model.BSIM4vtm0 * model.BSIM4njtsswdtemp;
    Nvtmrsswgd = model.BSIM4vtm0 * model.BSIM4njtsswgdtemp;
    Nvtmrsd = model.BSIM4vtm0 * model.BSIM4njtsdtemp;

    if ((model.BSIM4vtss - vbs_jct) < (model.BSIM4vtss * 1e-3))
    {   T9 = 1.0e3;
        T0 = - vbs_jct / Nvtmrss * T9;
        dexp(T0, T1, T10);
        dT1_dVb = T10 / Nvtmrss * T9;
    } else {
        T9 = 1.0 / (model.BSIM4vtss - vbs_jct);
        T0 = -vbs_jct / Nvtmrss * model.BSIM4vtss * T9;
        dT0_dVb = model.BSIM4vtss / Nvtmrss * (T9 + vbs_jct * T9 * T9) ;
        dexp(T0, T1, T10);
        dT1_dVb = T10 * dT0_dVb;
    }

    if ((model.BSIM4vtsd - vbd_jct) < (model.BSIM4vtsd * 1e-3) )
    {   T9 = 1.0e3;
        T0 = -vbd_jct / Nvtmrsd * T9;
        dexp(T0, T2, T10);
        dT2_dVb = T10 / Nvtmrsd * T9;
    } else {
        T9 = 1.0 / (model.BSIM4vtsd - vbd_jct);
        T0 = -vbd_jct / Nvtmrsd * model.BSIM4vtsd * T9;
        dT0_dVb = model.BSIM4vtsd / Nvtmrsd * (T9 + vbd_jct * T9 * T9) ;
        dexp(T0, T2, T10);
        dT2_dVb = T10 * dT0_dVb;
    }

    if ((model.BSIM4vtssws - vbs_jct) < (model.BSIM4vtssws * 1e-3) )
    {   T9 = 1.0e3;
        T0 = -vbs_jct / Nvtmrssws * T9;
        dexp(T0, T3, T10);
        dT3_dVb = T10 / Nvtmrssws * T9;
    } else {
        T9 = 1.0 / (model.BSIM4vtssws - vbs_jct);
        T0 = -vbs_jct / Nvtmrssws * model.BSIM4vtssws * T9;
        dT0_dVb = model.BSIM4vtssws / Nvtmrssws * (T9 + vbs_jct * T9 * T9) ;
        dexp(T0, T3, T10);
        dT3_dVb = T10 * dT0_dVb;
    }

    if ((model.BSIM4vtsswd - vbd_jct) < (model.BSIM4vtsswd * 1e-3) )
    {   T9 = 1.0e3;
        T0 = -vbd_jct / Nvtmrsswd * T9;
        dexp(T0, T4, T10);
        dT4_dVb = T10 / Nvtmrsswd * T9;
    } else {
        T9 = 1.0 / (model.BSIM4vtsswd - vbd_jct);
        T0 = -vbd_jct / Nvtmrsswd * model.BSIM4vtsswd * T9;
        dT0_dVb = model.BSIM4vtsswd / Nvtmrsswd * (T9 + vbd_jct * T9 * T9) ;
        dexp(T0, T4, T10);
        dT4_dVb = T10 * dT0_dVb;
    }

    if ((model.BSIM4vtsswgs - vbs_jct) < (model.BSIM4vtsswgs * 1e-3) )
    {   T9 = 1.0e3;
        T0 = -vbs_jct / Nvtmrsswgs * T9;
        dexp(T0, T5, T10);
        dT5_dVb = T10 / Nvtmrsswgs * T9;
    } else {
        T9 = 1.0 / (model.BSIM4vtsswgs - vbs_jct);
        T0 = -vbs_jct / Nvtmrsswgs * model.BSIM4vtsswgs * T9;
        dT0_dVb = model.BSIM4vtsswgs / Nvtmrsswgs * (T9 + vbs_jct * T9 * T9) ;
        dexp(T0, T5, T10);
        dT5_dVb = T10 * dT0_dVb;
    }

    if ((model.BSIM4vtsswgd - vbd_jct) < (model.BSIM4vtsswgd * 1e-3) )
    {   T9 = 1.0e3;
        T0 = -vbd_jct / Nvtmrsswgd * T9;
        dexp(T0, T6, T10);
        dT6_dVb = T10 / Nvtmrsswgd * T9;
    } else {
        T9 = 1.0 / (model.BSIM4vtsswgd - vbd_jct);
        T0 = -vbd_jct / Nvtmrsswgd * model.BSIM4vtsswgd * T9;
        dT0_dVb = model.BSIM4vtsswgd / Nvtmrsswgd * (T9 + vbd_jct * T9 * T9) ;
        dexp(T0, T6, T10);
        dT6_dVb = T10 * dT0_dVb;
    }

    instance.BSIM4gbs += instance.BSIM4SjctTempRevSatCur * dT1_dVb
            + instance.BSIM4SswTempRevSatCur * dT3_dVb
            + instance.BSIM4SswgTempRevSatCur * dT5_dVb;
    instance.BSIM4cbs -= instance.BSIM4SjctTempRevSatCur * (T1 - 1.0)
            + instance.BSIM4SswTempRevSatCur * (T3 - 1.0)
            + instance.BSIM4SswgTempRevSatCur * (T5 - 1.0);
    instance.BSIM4gbd += instance.BSIM4DjctTempRevSatCur * dT2_dVb
            + instance.BSIM4DswTempRevSatCur * dT4_dVb
            + instance.BSIM4DswgTempRevSatCur * dT6_dVb;
    instance.BSIM4cbd -= instance.BSIM4DjctTempRevSatCur * (T2 - 1.0)
            + instance.BSIM4DswTempRevSatCur * (T4 - 1.0)
            + instance.BSIM4DswgTempRevSatCur * (T6 - 1.0);

    /* End of diode DC model */

    if (vds >= 0.0)
    {   instance.BSIM4mode = 1;
        Vds = vds;
        Vgs = vgs;
        Vbs = vbs;
        Vdb = vds - vbs;  /* WDLiu: for GIDL */

    }
    else
    {   instance.BSIM4mode = -1;
        Vds = -vds;
        Vgs = vgd;
        Vbs = vbd;
        Vdb = -vbs;
    }


    /* dunga */
    if(model.BSIM4mtrlMod)
    {
        epsrox = 3.9;
        toxe = model.BSIM4eot;
        epssub = EPS0 * model.BSIM4epsrsub;
    }
    else
    {
        epsrox = model.BSIM4epsrox;
        toxe = model.BSIM4toxe;
        epssub = EPSSI;
    }


    T0 = Vbs - instance.BSIM4vbsc - 0.001;
    T1 = std::sqrt(T0 * T0 - 0.004 * instance.BSIM4vbsc);
    if (T0 >= 0.0)
    {   Vbseff = instance.BSIM4vbsc + 0.5 * (T0 + T1);
            dVbseff_dVb = 0.5 * (1.0 + T0 / T1);
    }
    else
    {   T2 = -0.002 / (T1 - T0);
        Vbseff = instance.BSIM4vbsc * (1.0 + T2);
        dVbseff_dVb = T2 * instance.BSIM4vbsc / T1;
    }

    /* JX: Correction to forward body bias  */
    T9 = 0.95 * pParam.BSIM4phi;
    T0 = T9 - Vbseff - 0.001;
    T1 = std::sqrt(T0 * T0 + 0.004 * T9);
    Vbseff = T9 - 0.5 * (T0 + T1);
    dVbseff_dVb *= 0.5 * (1.0 + T0 / T1);
    Phis = pParam.BSIM4phi - Vbseff;
    dPhis_dVb = -1.0;
    sqrtPhis = std::sqrt(Phis);
    dsqrtPhis_dVb = -0.5 / sqrtPhis;

    Xdep = pParam.BSIM4Xdep0 * sqrtPhis / pParam.BSIM4sqrtPhi;
    dXdep_dVb = (pParam.BSIM4Xdep0 / pParam.BSIM4sqrtPhi)
    * dsqrtPhis_dVb;

    Leff = pParam.BSIM4leff;
    Vtm = model.BSIM4vtm;
    Vtm0 = model.BSIM4vtm0;

    /* Vth Calculation */
    T3 = std::sqrt(Xdep);
    V0 = pParam.BSIM4vbi - pParam.BSIM4phi;

    T0 = pParam.BSIM4dvt2 * Vbseff;
    if (T0 >= - 0.5)
    {   T1 = 1.0 + T0;
        T2 = pParam.BSIM4dvt2;
    }
    else
    {   T4 = 1.0 / (3.0 + 8.0 * T0);
        T1 = (1.0 + 3.0 * T0) * T4;
        T2 = pParam.BSIM4dvt2 * T4 * T4;
    }
    lt1 = model.BSIM4factor1 * T3 * T1;
    dlt1_dVb = model.BSIM4factor1 * (0.5 / T3 * T1 * dXdep_dVb + T3 * T2);

    T0 = pParam.BSIM4dvt2w * Vbseff;
    if (T0 >= - 0.5)
    {   T1 = 1.0 + T0;
        T2 = pParam.BSIM4dvt2w;
    }
    else
    {   T4 = 1.0 / (3.0 + 8.0 * T0);
        T1 = (1.0 + 3.0 * T0) * T4;
        T2 = pParam.BSIM4dvt2w * T4 * T4;
    }
    ltw = model.BSIM4factor1 * T3 * T1;
    dltw_dVb = model.BSIM4factor1 * (0.5 / T3 * T1 * dXdep_dVb + T3 * T2);

    T0 = pParam.BSIM4dvt1 * Leff / lt1;
    if (T0 < EXP_THRESHOLD)
    {   T1 = std::exp(T0);
        T2 = T1 - 1.0;
        T3 = T2 * T2;
        T4 = T3 + 2.0 * T1 * MIN_EXP;
        Theta0 = T1 / T4;
        dT1_dVb = -T0 * T1 * dlt1_dVb / lt1;
        dTheta0_dVb = dT1_dVb * (T4 - 2.0 * T1 * (T2 + MIN_EXP)) / T4 / T4;
    }
    else
    {   Theta0 = 1.0 / (MAX_EXP - 2.0); /* 3.0 * MIN_EXP omitted */
        dTheta0_dVb = 0.0;
    }
    instance.BSIM4thetavth = pParam.BSIM4dvt0 * Theta0;
    Delt_vth = instance.BSIM4thetavth * V0;
    dDelt_vth_dVb = pParam.BSIM4dvt0 * dTheta0_dVb * V0;

    T0 = pParam.BSIM4dvt1w * pParam.BSIM4weff * Leff / ltw;
    if (T0 < EXP_THRESHOLD)
    {   T1 = std::exp(T0);
        T2 = T1 - 1.0;
        T3 = T2 * T2;
        T4 = T3 + 2.0 * T1 * MIN_EXP;
        T5 = T1 / T4;
        dT1_dVb = -T0 * T1 * dltw_dVb / ltw;
        dT5_dVb = dT1_dVb * (T4 - 2.0 * T1 * (T2 + MIN_EXP)) / T4 / T4;
    }
    else
    {   T5 = 1.0 / (MAX_EXP - 2.0); /* 3.0 * MIN_EXP omitted */
        dT5_dVb = 0.0;
    }
    T0 = pParam.BSIM4dvt0w * T5;
    T2 = T0 * V0;
    dT2_dVb = pParam.BSIM4dvt0w * dT5_dVb * V0;

    TempRatio =  CKTtemp / model.BSIM4tnom - 1.0;
    T0 = std::sqrt(1.0 + pParam.BSIM4lpe0 / Leff);
    T1 = pParam.BSIM4k1ox * (T0 - 1.0) * pParam.BSIM4sqrtPhi
        + (pParam.BSIM4kt1 + pParam.BSIM4kt1l / Leff
        + pParam.BSIM4kt2 * Vbseff) * TempRatio;
    Vth_NarrowW = toxe * pParam.BSIM4phi
            / (pParam.BSIM4weff + pParam.BSIM4w0);

    T3 = instance.BSIM4eta0 + pParam.BSIM4etab * Vbseff;
    if (T3 < 1.0e-4)
    {   T9 = 1.0 / (3.0 - 2.0e4 * T3);
        T3 = (2.0e-4 - T3) * T9;
        T4 = T9 * T9;
    }
    else
    {   T4 = 1.0;}
    dDIBL_Sft_dVd = T3 * pParam.BSIM4theta0vb0;
    DIBL_Sft = dDIBL_Sft_dVd * Vds;

    Lpe_Vb = std::sqrt(1.0 + pParam.BSIM4lpeb / Leff);

    Vth = model.BSIM4type * instance.BSIM4vth0 + (pParam.BSIM4k1ox * sqrtPhis
        - pParam.BSIM4k1 * pParam.BSIM4sqrtPhi) * Lpe_Vb
        - instance.BSIM4k2ox * Vbseff - Delt_vth - T2 + (pParam.BSIM4k3
        + pParam.BSIM4k3b * Vbseff) * Vth_NarrowW + T1 - DIBL_Sft;

    dVth_dVb = Lpe_Vb * pParam.BSIM4k1ox * dsqrtPhis_dVb - instance.BSIM4k2ox
            - dDelt_vth_dVb - dT2_dVb + pParam.BSIM4k3b * Vth_NarrowW
            - pParam.BSIM4etab * Vds * pParam.BSIM4theta0vb0 * T4
            + pParam.BSIM4kt2 * TempRatio;
    dVth_dVd = -dDIBL_Sft_dVd;


    /* Calculate n */
    tmp1 = epssub / Xdep;
    instance.BSIM4nstar = model.BSIM4vtm / Charge_q * (model.BSIM4coxe + tmp1 + pParam.BSIM4cit);
    tmp2 = pParam.BSIM4nfactor * tmp1;
    tmp3 = pParam.BSIM4cdsc + pParam.BSIM4cdscb * Vbseff
        + pParam.BSIM4cdscd * Vds;
    tmp4 = (tmp2 + tmp3 * Theta0 + pParam.BSIM4cit) / model.BSIM4coxe;
    if (tmp4 >= -0.5)
    {   n = 1.0 + tmp4;
        dn_dVb = (-tmp2 / Xdep * dXdep_dVb + tmp3 * dTheta0_dVb
                    + pParam.BSIM4cdscb * Theta0) / model.BSIM4coxe;
        dn_dVd = pParam.BSIM4cdscd * Theta0 / model.BSIM4coxe;
    }
    else
    {   T0 = 1.0 / (3.0 + 8.0 * tmp4);
        n = (1.0 + 3.0 * tmp4) * T0;
        T0 *= T0;
        dn_dVb = (-tmp2 / Xdep * dXdep_dVb + tmp3 * dTheta0_dVb
                    + pParam.BSIM4cdscb * Theta0) / model.BSIM4coxe * T0;
        dn_dVd = pParam.BSIM4cdscd * Theta0 / model.BSIM4coxe * T0;
    }


    /* Vth correction for Pocket implant */
    if (pParam.BSIM4dvtp0 > 0.0)
    {   T0 = -pParam.BSIM4dvtp1 * Vds;
        if (T0 < -EXP_THRESHOLD)
        {   T2 = MIN_EXP;
            dT2_dVd = 0.0;
        }
        else
        {   T2 = std::exp(T0);
            dT2_dVd = -pParam.BSIM4dvtp1 * T2;
        }

        T3 = Leff + pParam.BSIM4dvtp0 * (1.0 + T2);
        dT3_dVd = pParam.BSIM4dvtp0 * dT2_dVd;
        if (model.BSIM4tempMod < 2)
        {
            T4 = Vtm * std::log(Leff / T3);
            dT4_dVd = -Vtm * dT3_dVd / T3;
        }
        else
        {
            T4 = model.BSIM4vtm0 * std::log(Leff / T3);
            dT4_dVd = -model.BSIM4vtm0 * dT3_dVd / T3;
        }
        dDITS_Sft_dVd = dn_dVd * T4 + n * dT4_dVd;
        dDITS_Sft_dVb = T4 * dn_dVb;

        Vth -= n * T4;
        dVth_dVd -= dDITS_Sft_dVd;
        dVth_dVb -= dDITS_Sft_dVb;
    }

    /* v4.7 DITS_SFT2  */
    if ((pParam.BSIM4dvtp4  == 0.0) || (pParam.BSIM4dvtp2factor == 0.0)) {
      T0 = 0.0;
      DITS_Sft2 = 0.0;
    }
    else
    {
      //T0 = std::exp(2.0 * pParam.BSIM4dvtp4 * Vds);   /* beta code */
      T1 = 2.0 * pParam.BSIM4dvtp4 * Vds;
      dexp(T1, T0, T10);
      DITS_Sft2 = pParam.BSIM4dvtp2factor * (T0-1) / (T0+1);
      //dDITS_Sft2_dVd = pParam.BSIM4dvtp2factor * pParam.BSIM4dvtp4 * 4.0 * T0 / ((T0+1) * (T0+1));   /* beta code */
      dDITS_Sft2_dVd = pParam.BSIM4dvtp2factor * pParam.BSIM4dvtp4 * 4.0 * T10 / ((T0+1) * (T0+1));
      Vth -= DITS_Sft2;
      dVth_dVd -= dDITS_Sft2_dVd;
    }

    instance.BSIM4von = Vth;

    /* Poly Gate Si Depletion Effect */
    T0 = instance.BSIM4vfb + pParam.BSIM4phi;
    if(model.BSIM4mtrlMod == 0)
        T1 = EPSSI;
    else
        T1 = model.BSIM4epsrgate * EPS0;


    BSIM4polyDepletion(T0, pParam.BSIM4ngate, T1, model.BSIM4coxe, vgs, &vgs_eff, &dvgs_eff_dvg);

    BSIM4polyDepletion(T0, pParam.BSIM4ngate, T1, model.BSIM4coxe, vgd, &vgd_eff, &dvgd_eff_dvg);

    if(instance.BSIM4mode>0) {
        Vgs_eff = vgs_eff;
        dVgs_eff_dVg = dvgs_eff_dvg;
    } else {
        Vgs_eff = vgd_eff;
        dVgs_eff_dVg = dvgd_eff_dvg;
    }
    instance.BSIM4vgs_eff = vgs_eff;
    instance.BSIM4vgd_eff = vgd_eff;
    instance.BSIM4dvgs_eff_dvg = dvgs_eff_dvg;
    instance.BSIM4dvgd_eff_dvg = dvgd_eff_dvg;


    Vgst = Vgs_eff - Vth;

    /* Calculate Vgsteff */
    T0 = n * Vtm;
    T1 = pParam.BSIM4mstar * Vgst;
    T2 = T1 / T0;
    if (T2 > EXP_THRESHOLD)
    {   T10 = T1;
        dT10_dVg = pParam.BSIM4mstar * dVgs_eff_dVg;
        dT10_dVd = -dVth_dVd * pParam.BSIM4mstar;
        dT10_dVb = -dVth_dVb * pParam.BSIM4mstar;
    }
    else if (T2 < -EXP_THRESHOLD)
    {   T10 = Vtm * std::log(1.0 + MIN_EXP);
        dT10_dVg = 0.0;
        dT10_dVd = T10 * dn_dVd;
        dT10_dVb = T10 * dn_dVb;
        T10 *= n;
    }
    else
    {   ExpVgst = std::exp(T2);
        T3 = Vtm * std::log(1.0 + ExpVgst);
        T10 = n * T3;
        dT10_dVg = pParam.BSIM4mstar * ExpVgst / (1.0 + ExpVgst);
        dT10_dVb = T3 * dn_dVb - dT10_dVg * (dVth_dVb + Vgst * dn_dVb / n);
        dT10_dVd = T3 * dn_dVd - dT10_dVg * (dVth_dVd + Vgst * dn_dVd / n);
        dT10_dVg *= dVgs_eff_dVg;
    }

    T1 = pParam.BSIM4voffcbn - (1.0 - pParam.BSIM4mstar) * Vgst;
    T2 = T1 / T0;
    if (T2 < -EXP_THRESHOLD)
    {   T3 = model.BSIM4coxe * MIN_EXP / pParam.BSIM4cdep0;
        T9 = pParam.BSIM4mstar + T3 * n;
        dT9_dVg = 0.0;
        dT9_dVd = dn_dVd * T3;
        dT9_dVb = dn_dVb * T3;
    }
    else if (T2 > EXP_THRESHOLD)
    {   T3 = model.BSIM4coxe * MAX_EXP / pParam.BSIM4cdep0;
        T9 = pParam.BSIM4mstar + T3 * n;
        dT9_dVg = 0.0;
        dT9_dVd = dn_dVd * T3;
        dT9_dVb = dn_dVb * T3;
    }
    else
    {   ExpVgst = std::exp(T2);
        T3 = model.BSIM4coxe / pParam.BSIM4cdep0;
        T4 = T3 * ExpVgst;
        T5 = T1 * T4 / T0;
        T9 = pParam.BSIM4mstar + n * T4;
        dT9_dVg = T3 * (pParam.BSIM4mstar - 1.0) * ExpVgst / Vtm;
        dT9_dVb = T4 * dn_dVb - dT9_dVg * dVth_dVb - T5 * dn_dVb;
        dT9_dVd = T4 * dn_dVd - dT9_dVg * dVth_dVd - T5 * dn_dVd;
        dT9_dVg *= dVgs_eff_dVg;
    }
    instance.BSIM4Vgsteff = Vgsteff = T10 / T9;
    T11 = T9 * T9;
    dVgsteff_dVg = (T9 * dT10_dVg - T10 * dT9_dVg) / T11;
    dVgsteff_dVd = (T9 * dT10_dVd - T10 * dT9_dVd) / T11;
    dVgsteff_dVb = (T9 * dT10_dVb - T10 * dT9_dVb) / T11;

    /* Calculate Effective Channel Geometry */
    T9 = sqrtPhis - pParam.BSIM4sqrtPhi;
    Weff = pParam.BSIM4weff - 2.0 * (pParam.BSIM4dwg * Vgsteff
        + pParam.BSIM4dwb * T9);
    dWeff_dVg = -2.0 * pParam.BSIM4dwg;
    dWeff_dVb = -2.0 * pParam.BSIM4dwb * dsqrtPhis_dVb;

    if (Weff < 2.0e-8) /* to avoid the discontinuity problem due to Weff*/
    {   T0 = 1.0 / (6.0e-8 - 2.0 * Weff);
        Weff = 2.0e-8 * (4.0e-8 - Weff) * T0;
        T0 *= T0 * 4.0e-16;
        dWeff_dVg *= T0;
        dWeff_dVb *= T0;
    }

    if (model.BSIM4rdsMod == 1)
        Rds = dRds_dVg = dRds_dVb = 0.0;
    else
    {   T0 = 1.0 + pParam.BSIM4prwg * Vgsteff;
        dT0_dVg = -pParam.BSIM4prwg / T0 / T0;
        T1 = pParam.BSIM4prwb * T9;
        dT1_dVb = pParam.BSIM4prwb * dsqrtPhis_dVb;

        T2 = 1.0 / T0 + T1;
        T3 = T2 + std::sqrt(T2 * T2 + 0.01); /* 0.01 = 4.0 * 0.05 * 0.05 */
        dT3_dVg = 1.0 + T2 / (T3 - T2);
        dT3_dVb = dT3_dVg * dT1_dVb;
        dT3_dVg *= dT0_dVg;

        T4 = pParam.BSIM4rds0 * 0.5;
        Rds = pParam.BSIM4rdswmin + T3 * T4;
        dRds_dVg = T4 * dT3_dVg;
        dRds_dVb = T4 * dT3_dVb;

        if (Rds > 0.0)
            instance.BSIM4grdsw = 1.0 / Rds * instance.BSIM4nf; /*4.6.2*/
        else
            instance.BSIM4grdsw = 0.0;
    }

    /* Calculate Abulk */
    T9 = 0.5 * pParam.BSIM4k1ox * Lpe_Vb / sqrtPhis;
    T1 = T9 + instance.BSIM4k2ox - pParam.BSIM4k3b * Vth_NarrowW;
    dT1_dVb = -T9 / sqrtPhis * dsqrtPhis_dVb;

    T9 = std::sqrt(pParam.BSIM4xj * Xdep);
    tmp1 = Leff + 2.0 * T9;
    T5 = Leff / tmp1;
    tmp2 = pParam.BSIM4a0 * T5;
    tmp3 = pParam.BSIM4weff + pParam.BSIM4b1;
    tmp4 = pParam.BSIM4b0 / tmp3;
    T2 = tmp2 + tmp4;
    dT2_dVb = -T9 / tmp1 / Xdep * dXdep_dVb;
    T6 = T5 * T5;
    T7 = T5 * T6;

    Abulk0 = 1.0 + T1 * T2;
    dAbulk0_dVb = T1 * tmp2 * dT2_dVb + T2 * dT1_dVb;

    T8 = pParam.BSIM4ags * pParam.BSIM4a0 * T7;
    dAbulk_dVg = -T1 * T8;
    Abulk = Abulk0 + dAbulk_dVg * Vgsteff;
    dAbulk_dVb = dAbulk0_dVb - T8 * Vgsteff * (dT1_dVb
        + 3.0 * T1 * dT2_dVb);

    if (Abulk0 < 0.1) /* added to avoid the problems caused by Abulk0 */
    {   T9 = 1.0 / (3.0 - 20.0 * Abulk0);
        Abulk0 = (0.2 - Abulk0) * T9;
        dAbulk0_dVb *= T9 * T9;
    }

    if (Abulk < 0.1)
    {   T9 = 1.0 / (3.0 - 20.0 * Abulk);
        Abulk = (0.2 - Abulk) * T9;
        T10 = T9 * T9;
        dAbulk_dVb *= T10;
        dAbulk_dVg *= T10;
    }
    instance.BSIM4Abulk = Abulk;

    T2 = pParam.BSIM4keta * Vbseff;
    if (T2 >= -0.9)
    {   T0 = 1.0 / (1.0 + T2);
        dT0_dVb = -pParam.BSIM4keta * T0 * T0;
    }
    else
    {   T1 = 1.0 / (0.8 + T2);
        T0 = (17.0 + 20.0 * T2) * T1;
        dT0_dVb = -pParam.BSIM4keta * T1 * T1;
    }
    dAbulk_dVg *= T0;
    dAbulk_dVb = dAbulk_dVb * T0 + Abulk * dT0_dVb;
    dAbulk0_dVb = dAbulk0_dVb * T0 + Abulk0 * dT0_dVb;
    Abulk *= T0;
    Abulk0 *= T0;

    /* Mobility calculation */
    if (model.BSIM4mtrlMod && model.BSIM4mtrlCompatMod == 0)
        T14 = 2.0 * model.BSIM4type *(model.BSIM4phig - model.BSIM4easub - 0.5*model.BSIM4Eg0 + 0.45);
    else
        T14 = 0.0;

    if (model.BSIM4mobMod == 0)
    {   T0 = Vgsteff + Vth + Vth - T14;
        T2 = pParam.BSIM4ua + pParam.BSIM4uc * Vbseff;
        T3 = T0 / toxe;
        T12 = std::sqrt(Vth * Vth + 0.0001);
        T9 = 1.0/(Vgsteff + 2*T12);
        T10 = T9*toxe;
        T8 = pParam.BSIM4ud * T10 * T10 * Vth;
        T6 = T8 * Vth;
        T5 = T3 * (T2 + pParam.BSIM4ub * T3) + T6;
        T7 = - 2.0 * T6 * T9;
        T11 = T7 * Vth/T12;
        dDenomi_dVg = (T2 + 2.0 * pParam.BSIM4ub * T3) / toxe;
        T13 = 2.0 * (dDenomi_dVg + T11 + T8);
        dDenomi_dVd = T13 * dVth_dVd;
        dDenomi_dVb = T13 * dVth_dVb + pParam.BSIM4uc * T3;
        dDenomi_dVg+= T7;
    }
    else if (model.BSIM4mobMod == 1)
    {   T0 = Vgsteff + Vth + Vth - T14;
        T2 = 1.0 + pParam.BSIM4uc * Vbseff;
        T3 = T0 / toxe;
        T4 = T3 * (pParam.BSIM4ua + pParam.BSIM4ub * T3);
        T12 = std::sqrt(Vth * Vth + 0.0001);
        T9 = 1.0/(Vgsteff + 2*T12);
        T10 = T9*toxe;
        T8 = pParam.BSIM4ud * T10 * T10 * Vth;
        T6 = T8 * Vth;
        T5 = T4 * T2 + T6;
        T7 = - 2.0 * T6 * T9;
        T11 = T7 * Vth/T12;
        dDenomi_dVg = (pParam.BSIM4ua + 2.0 * pParam.BSIM4ub * T3) * T2 / toxe;
        T13 = 2.0 * (dDenomi_dVg + T11 + T8);
        dDenomi_dVd = T13 * dVth_dVd;
        dDenomi_dVb = T13 * dVth_dVb + pParam.BSIM4uc * T4;
        dDenomi_dVg+= T7;
    }
    else if (model.BSIM4mobMod == 2)
    {   T0 = (Vgsteff + instance.BSIM4vtfbphi1) / toxe;
        T1 = std::exp(pParam.BSIM4eu * std::log(T0));
        dT1_dVg = T1 * pParam.BSIM4eu / T0 / toxe;
        T2 = pParam.BSIM4ua + pParam.BSIM4uc * Vbseff;

        T12 = std::sqrt(Vth * Vth + 0.0001);
        T9 = 1.0/(Vgsteff + 2*T12);
        T10 = T9*toxe;
        T8 = pParam.BSIM4ud * T10 * T10 * Vth;
        T6 = T8 * Vth;
        T5 = T1 * T2 + T6;
        T7 = - 2.0 * T6 * T9;
        T11 = T7 * Vth/T12;
        dDenomi_dVg = T2 * dT1_dVg + T7;
        T13 = 2.0 * (T11 + T8);
        dDenomi_dVd = T13 * dVth_dVd;
        dDenomi_dVb = T13 * dVth_dVb + T1 * pParam.BSIM4uc;
    }
    else if (model.BSIM4mobMod == 4) /* Synopsys 08/30/2013 add */
    {
        T0 = Vgsteff + instance.BSIM4vtfbphi1 - T14;
        T2 = pParam.BSIM4ua + pParam.BSIM4uc * Vbseff;
        T3 = T0 / toxe;
        T12 = std::sqrt(instance.BSIM4vtfbphi1*instance.BSIM4vtfbphi1 + 0.0001);
        T9 = 1.0/(Vgsteff + 2*T12);
        T10 = T9*toxe;
        T8 = pParam.BSIM4ud * T10 * T10 * instance.BSIM4vtfbphi1;
        T6 = T8 * instance.BSIM4vtfbphi1;
        T5 = T3 * (T2 + pParam.BSIM4ub * T3) + T6;
        T7 = - 2.0 * T6 * T9;
        dDenomi_dVg = (T2 + 2.0 * pParam.BSIM4ub * T3) / toxe;
        dDenomi_dVd = 0.0;
        dDenomi_dVb = pParam.BSIM4uc * T3;
        dDenomi_dVg+= T7;
    }
    else if (model.BSIM4mobMod == 5) /* Synopsys 08/30/2013 add */
    {
        T0 = Vgsteff + instance.BSIM4vtfbphi1 - T14;
        T2 = 1.0 + pParam.BSIM4uc * Vbseff;
        T3 = T0 / toxe;
        T4 = T3 * (pParam.BSIM4ua + pParam.BSIM4ub * T3);
        T12 = std::sqrt(instance.BSIM4vtfbphi1 * instance.BSIM4vtfbphi1 + 0.0001);
        T9 = 1.0/(Vgsteff + 2*T12);
        T10 = T9*toxe;
        T8 = pParam.BSIM4ud * T10 * T10 * instance.BSIM4vtfbphi1;
        T6 = T8 * instance.BSIM4vtfbphi1;
        T5 = T4 * T2 + T6;
        T7 = - 2.0 * T6 * T9;
        dDenomi_dVg = (pParam.BSIM4ua + 2.0 * pParam.BSIM4ub * T3) * T2
                    / toxe;
        dDenomi_dVd = 0.0;
        dDenomi_dVb = pParam.BSIM4uc * T4;
        dDenomi_dVg+= T7;
    }
    else if (model.BSIM4mobMod == 6) /* Synopsys 08/30/2013 modify */
    {   T0 = (Vgsteff + instance.BSIM4vtfbphi1) / toxe;
        T1 = std::exp(pParam.BSIM4eu * std::log(T0));
        dT1_dVg = T1 * pParam.BSIM4eu / T0 / toxe;
        T2 = pParam.BSIM4ua + pParam.BSIM4uc * Vbseff;

        T12 = std::sqrt(instance.BSIM4vtfbphi1 * instance.BSIM4vtfbphi1 + 0.0001);
        T9 = 1.0/(Vgsteff + 2*T12);
        T10 = T9*toxe;
        T8 = pParam.BSIM4ud * T10 * T10 * instance.BSIM4vtfbphi1;
        T6 = T8 * instance.BSIM4vtfbphi1;
        T5 = T1 * T2 + T6;
        T7 = - 2.0 * T6 * T9;
        dDenomi_dVg = T2 * dT1_dVg + T7;
        dDenomi_dVd = 0;
        dDenomi_dVb = T1 * pParam.BSIM4uc;
    }

    /*high K mobility*/
    else
    {
        /*univsersal mobility*/
        T0 = (Vgsteff + instance.BSIM4vtfbphi1)* 1.0e-8 / toxe/6.0;
        T1 = std::exp(pParam.BSIM4eu * std::log(T0));
        dT1_dVg = T1 * pParam.BSIM4eu * 1.0e-8/ T0 / toxe/6.0;
        T2 = pParam.BSIM4ua + pParam.BSIM4uc * Vbseff;

        /*Coulombic*/
        VgsteffVth = pParam.BSIM4VgsteffVth;

        T10 = std::exp(pParam.BSIM4ucs * std::log(0.5 + 0.5 * Vgsteff/VgsteffVth));
        T11 =  pParam.BSIM4ud/T10;
        dT11_dVg = - 0.5 * pParam.BSIM4ucs * T11 /(0.5 + 0.5*Vgsteff/VgsteffVth)/VgsteffVth;

        dDenomi_dVg = T2 * dT1_dVg + dT11_dVg;
        dDenomi_dVd = 0.0;
        dDenomi_dVb = T1 * pParam.BSIM4uc;

        T5 = T1 * T2 + T11;
    }

    if (T5 >= -0.8)
    {   Denomi = 1.0 + T5; }
    else
    {   
        T9 = 1.0 / (7.0 + 10.0 * T5);
        Denomi = (0.6 + T5) * T9;
        T9 *= T9;
        dDenomi_dVg *= T9;
        dDenomi_dVd *= T9;
        dDenomi_dVb *= T9;
    }


    instance.BSIM4ueff = ueff = instance.BSIM4u0temp / Denomi;
    T9 = -ueff / Denomi;
    dueff_dVg = T9 * dDenomi_dVg;
    dueff_dVd = T9 * dDenomi_dVd;
    dueff_dVb = T9 * dDenomi_dVb;

    /* Saturation Drain Voltage  Vdsat */
    WVCox = Weff * instance.BSIM4vsattemp * model.BSIM4coxe;
    WVCoxRds = WVCox * Rds;

    Esat = 2.0 * instance.BSIM4vsattemp / ueff;
    instance.BSIM4EsatL = EsatL = Esat * Leff;
    T0 = -EsatL /ueff;
    dEsatL_dVg = T0 * dueff_dVg;
    dEsatL_dVd = T0 * dueff_dVd;
    dEsatL_dVb = T0 * dueff_dVb;

    /* Sqrt() */
    a1 = pParam.BSIM4a1;
    if (a1 == 0.0)
    {   Lambda = pParam.BSIM4a2;
        dLambda_dVg = 0.0;
    }
    else if (a1 > 0.0)
    {   T0 = 1.0 - pParam.BSIM4a2;
        T1 = T0 - pParam.BSIM4a1 * Vgsteff - 0.0001;
        T2 = std::sqrt(T1 * T1 + 0.0004 * T0);
        Lambda = pParam.BSIM4a2 + T0 - 0.5 * (T1 + T2);
        dLambda_dVg = 0.5 * pParam.BSIM4a1 * (1.0 + T1 / T2);
    }
    else
    {   T1 = pParam.BSIM4a2 + pParam.BSIM4a1 * Vgsteff - 0.0001;
        T2 = std::sqrt(T1 * T1 + 0.0004 * pParam.BSIM4a2);
        Lambda = 0.5 * (T1 + T2);
        dLambda_dVg = 0.5 * pParam.BSIM4a1 * (1.0 + T1 / T2);
    }

    Vgst2Vtm = Vgsteff + 2.0 * Vtm;
    if (Rds > 0)
    {   tmp2 = dRds_dVg / Rds + dWeff_dVg / Weff;
        tmp3 = dRds_dVb / Rds + dWeff_dVb / Weff;
    }
    else
    {   tmp2 = dWeff_dVg / Weff;
        tmp3 = dWeff_dVb / Weff;
    }

    if ((Rds == 0.0) && (Lambda == 1.0))
    {   T0 = 1.0 / (Abulk * EsatL + Vgst2Vtm);
        tmp1 = 0.0;
        T1 = T0 * T0;
        T2 = Vgst2Vtm * T0;
        T3 = EsatL * Vgst2Vtm;
        Vdsat = T3 * T0;

        dT0_dVg = -(Abulk * dEsatL_dVg + EsatL * dAbulk_dVg + 1.0) * T1;
        dT0_dVd = -(Abulk * dEsatL_dVd) * T1;
        dT0_dVb = -(Abulk * dEsatL_dVb + dAbulk_dVb * EsatL) * T1;

        dVdsat_dVg = T3 * dT0_dVg + T2 * dEsatL_dVg + EsatL * T0;
        dVdsat_dVd = T3 * dT0_dVd + T2 * dEsatL_dVd;
        dVdsat_dVb = T3 * dT0_dVb + T2 * dEsatL_dVb;
    }
    else
    {   tmp1 = dLambda_dVg / (Lambda * Lambda);
        T9 = Abulk * WVCoxRds;
        T8 = Abulk * T9;
        T7 = Vgst2Vtm * T9;
        T6 = Vgst2Vtm * WVCoxRds;
        T0 = 2.0 * Abulk * (T9 - 1.0 + 1.0 / Lambda);
        dT0_dVg = 2.0 * (T8 * tmp2 - Abulk * tmp1
            + (2.0 * T9 + 1.0 / Lambda - 1.0) * dAbulk_dVg);
        dT0_dVb = 2.0 * (T8 * (2.0 / Abulk * dAbulk_dVb + tmp3)
            + (1.0 / Lambda - 1.0) * dAbulk_dVb);
        dT0_dVd = 0.0;
        T1 = Vgst2Vtm * (2.0 / Lambda - 1.0) + Abulk * EsatL + 3.0 * T7;

        dT1_dVg = (2.0 / Lambda - 1.0) - 2.0 * Vgst2Vtm * tmp1
            + Abulk * dEsatL_dVg + EsatL * dAbulk_dVg + 3.0 * (T9
            + T7 * tmp2 + T6 * dAbulk_dVg);
        dT1_dVb = Abulk * dEsatL_dVb + EsatL * dAbulk_dVb
            + 3.0 * (T6 * dAbulk_dVb + T7 * tmp3);
        dT1_dVd = Abulk * dEsatL_dVd;

        T2 = Vgst2Vtm * (EsatL + 2.0 * T6);
        dT2_dVg = EsatL + Vgst2Vtm * dEsatL_dVg
            + T6 * (4.0 + 2.0 * Vgst2Vtm * tmp2);
        dT2_dVb = Vgst2Vtm * (dEsatL_dVb + 2.0 * T6 * tmp3);
        dT2_dVd = Vgst2Vtm * dEsatL_dVd;

        T3 = std::sqrt(T1 * T1 - 2.0 * T0 * T2);
        Vdsat = (T1 - T3) / T0;

        dT3_dVg = (T1 * dT1_dVg - 2.0 * (T0 * dT2_dVg + T2 * dT0_dVg))
            / T3;
        dT3_dVd = (T1 * dT1_dVd - 2.0 * (T0 * dT2_dVd + T2 * dT0_dVd))
            / T3;
        dT3_dVb = (T1 * dT1_dVb - 2.0 * (T0 * dT2_dVb + T2 * dT0_dVb))
            / T3;

        dVdsat_dVg = (dT1_dVg - (T1 * dT1_dVg - dT0_dVg * T2
            - T0 * dT2_dVg) / T3 - Vdsat * dT0_dVg) / T0;
        dVdsat_dVb = (dT1_dVb - (T1 * dT1_dVb - dT0_dVb * T2
            - T0 * dT2_dVb) / T3 - Vdsat * dT0_dVb) / T0;
        dVdsat_dVd = (dT1_dVd - (T1 * dT1_dVd - T0 * dT2_dVd) / T3) / T0;
    }
    instance.BSIM4vdsat = Vdsat;

    /* Calculate Vdseff */
    T1 = Vdsat - Vds - pParam.BSIM4delta;
    dT1_dVg = dVdsat_dVg;
    dT1_dVd = dVdsat_dVd - 1.0;
    dT1_dVb = dVdsat_dVb;

    T2 = std::sqrt(T1 * T1 + 4.0 * pParam.BSIM4delta * Vdsat);
    T0 = T1 / T2;
    T9 = 2.0 * pParam.BSIM4delta;
    T3 = T9 / T2;
    dT2_dVg = T0 * dT1_dVg + T3 * dVdsat_dVg;
    dT2_dVd = T0 * dT1_dVd + T3 * dVdsat_dVd;
    dT2_dVb = T0 * dT1_dVb + T3 * dVdsat_dVb;

    if (T1 >= 0.0)
    {   Vdseff = Vdsat - 0.5 * (T1 + T2);
        dVdseff_dVg = dVdsat_dVg - 0.5 * (dT1_dVg + dT2_dVg);
            dVdseff_dVd = dVdsat_dVd - 0.5 * (dT1_dVd + dT2_dVd);
            dVdseff_dVb = dVdsat_dVb - 0.5 * (dT1_dVb + dT2_dVb);
    }
    else
    {   T4 = T9 / (T2 - T1);
        T5 = 1.0 - T4;
        T6 = Vdsat * T4 / (T2 - T1);
        Vdseff = Vdsat * T5;
        dVdseff_dVg = dVdsat_dVg * T5 + T6 * (dT2_dVg - dT1_dVg);
        dVdseff_dVd = dVdsat_dVd * T5 + T6 * (dT2_dVd - dT1_dVd);
        dVdseff_dVb = dVdsat_dVb * T5 + T6 * (dT2_dVb - dT1_dVb);
    }

    if (Vds == 0.0)
    {  Vdseff = 0.0;
        dVdseff_dVg = 0.0;
        dVdseff_dVb = 0.0;
    }

    if (Vdseff > Vds)
        Vdseff = Vds;
    diffVds = Vds - Vdseff;
    instance.BSIM4Vdseff = Vdseff;

    /* Velocity Overshoot */
    if((model.BSIM4lambdaGiven) && (model.BSIM4lambda > 0.0) )
    {
        T1 =  Leff * ueff;
        T2 = pParam.BSIM4lambda / T1;
        T3 = -T2 / T1 * Leff;
        dT2_dVd = T3 * dueff_dVd;
        dT2_dVg = T3 * dueff_dVg;
        dT2_dVb = T3 * dueff_dVb;
        T5 = 1.0 / (Esat * pParam.BSIM4litl);
        T4 = -T5 / EsatL;
        dT5_dVg = dEsatL_dVg * T4;
        dT5_dVd = dEsatL_dVd * T4;
        dT5_dVb = dEsatL_dVb * T4;
        T6 = 1.0 + diffVds  * T5;
        dT6_dVg = dT5_dVg * diffVds - dVdseff_dVg * T5;
        dT6_dVd = dT5_dVd * diffVds + (1.0 - dVdseff_dVd) * T5;
        dT6_dVb = dT5_dVb * diffVds - dVdseff_dVb * T5;
        T7 = 2.0 / (T6 * T6 + 1.0);
        T8 = 1.0 - T7;
        T9 = T6 * T7 * T7;
        dT8_dVg = T9 * dT6_dVg;
        dT8_dVd = T9 * dT6_dVd;
        dT8_dVb = T9 * dT6_dVb;
        T10 = 1.0 + T2 * T8;
        dT10_dVg = dT2_dVg * T8 + T2 * dT8_dVg;
        dT10_dVd = dT2_dVd * T8 + T2 * dT8_dVd;
        dT10_dVb = dT2_dVb * T8 + T2 * dT8_dVb;
        if(T10 == 1.0)
            dT10_dVg = dT10_dVd = dT10_dVb = 0.0;

        dEsatL_dVg *= T10;
        dEsatL_dVg += EsatL * dT10_dVg;
        dEsatL_dVd *= T10;
        dEsatL_dVd += EsatL * dT10_dVd;
        dEsatL_dVb *= T10;
        dEsatL_dVb += EsatL * dT10_dVb;
        EsatL *= T10;
        Esat = EsatL / Leff;  /* bugfix by Wenwei Yang (4.6.4) */
        instance.BSIM4EsatL = EsatL;
    }

    /* Calculate Vasat */
    tmp4 = 1.0 - 0.5 * Abulk * Vdsat / Vgst2Vtm;
    T9 = WVCoxRds * Vgsteff;
    T8 = T9 / Vgst2Vtm;
    T0 = EsatL + Vdsat + 2.0 * T9 * tmp4;

    T7 = 2.0 * WVCoxRds * tmp4;
    dT0_dVg = dEsatL_dVg + dVdsat_dVg + T7 * (1.0 + tmp2 * Vgsteff)
        - T8 * (Abulk * dVdsat_dVg - Abulk * Vdsat / Vgst2Vtm
        + Vdsat * dAbulk_dVg);

    dT0_dVb = dEsatL_dVb + dVdsat_dVb + T7 * tmp3 * Vgsteff
        - T8 * (dAbulk_dVb * Vdsat + Abulk * dVdsat_dVb);
    dT0_dVd = dEsatL_dVd + dVdsat_dVd - T8 * Abulk * dVdsat_dVd;

    T9 = WVCoxRds * Abulk;
    T1 = 2.0 / Lambda - 1.0 + T9;
    dT1_dVg = -2.0 * tmp1 +  WVCoxRds * (Abulk * tmp2 + dAbulk_dVg);
    dT1_dVb = dAbulk_dVb * WVCoxRds + T9 * tmp3;

    Vasat = T0 / T1;
    dVasat_dVg = (dT0_dVg - Vasat * dT1_dVg) / T1;
    dVasat_dVb = (dT0_dVb - Vasat * dT1_dVb) / T1;
    dVasat_dVd = dT0_dVd / T1;

    /* Calculate Idl first */

    tmp1 = instance.BSIM4vtfbphi2;
    tmp2 = 2.0e8 * instance.BSIM4toxp;
    dT0_dVg = 1.0 / tmp2;
    T0 = (Vgsteff + tmp1) * dT0_dVg;

    tmp3 = std::exp(model.BSIM4bdos * 0.7 * std::log(T0));
    T1 = 1.0 + tmp3;
    T2 = model.BSIM4bdos * 0.7 * tmp3 / T0;
    Tcen = model.BSIM4ados * 1.9e-9 / T1;
    dTcen_dVg = -Tcen * T2 * dT0_dVg / T1;

    Coxeff = epssub * instance.BSIM4coxp
            / (epssub + instance.BSIM4coxp * Tcen);
    instance.BSIM4Coxeff = Coxeff;
    dCoxeff_dVg = -Coxeff * Coxeff * dTcen_dVg / epssub;

    CoxeffWovL = Coxeff * Weff / Leff;
    beta = ueff * CoxeffWovL;
    T3 = ueff / Leff;
    dbeta_dVg = CoxeffWovL * dueff_dVg + T3
            * (Weff * dCoxeff_dVg + Coxeff * dWeff_dVg);
    dbeta_dVd = CoxeffWovL * dueff_dVd;
    dbeta_dVb = CoxeffWovL * dueff_dVb + T3 * Coxeff * dWeff_dVb;

    instance.BSIM4AbovVgst2Vtm = Abulk / Vgst2Vtm;
    T0 = 1.0 - 0.5 * Vdseff * instance.BSIM4AbovVgst2Vtm;
    dT0_dVg = -0.5 * (Abulk * dVdseff_dVg
            - Abulk * Vdseff / Vgst2Vtm + Vdseff * dAbulk_dVg) / Vgst2Vtm;
    dT0_dVd = -0.5 * Abulk * dVdseff_dVd / Vgst2Vtm;
    dT0_dVb = -0.5 * (Abulk * dVdseff_dVb + dAbulk_dVb * Vdseff)
            / Vgst2Vtm;

    fgche1 = Vgsteff * T0;
    dfgche1_dVg = Vgsteff * dT0_dVg + T0;
    dfgche1_dVd = Vgsteff * dT0_dVd;
    dfgche1_dVb = Vgsteff * dT0_dVb;

    T9 = Vdseff / EsatL;
    fgche2 = 1.0 + T9;
    dfgche2_dVg = (dVdseff_dVg - T9 * dEsatL_dVg) / EsatL;
    dfgche2_dVd = (dVdseff_dVd - T9 * dEsatL_dVd) / EsatL;
    dfgche2_dVb = (dVdseff_dVb - T9 * dEsatL_dVb) / EsatL;

    gche = beta * fgche1 / fgche2;
    dgche_dVg = (beta * dfgche1_dVg + fgche1 * dbeta_dVg
            - gche * dfgche2_dVg) / fgche2;
    dgche_dVd = (beta * dfgche1_dVd + fgche1 * dbeta_dVd
            - gche * dfgche2_dVd) / fgche2;
    dgche_dVb = (beta * dfgche1_dVb + fgche1 * dbeta_dVb
            - gche * dfgche2_dVb) / fgche2;

    T0 = 1.0 + gche * Rds;
    Idl = gche / T0;
    T1 = (1.0 - Idl * Rds) / T0;
    T2 = Idl * Idl;
    dIdl_dVg = T1 * dgche_dVg - T2 * dRds_dVg;
    dIdl_dVd = T1 * dgche_dVd;
    dIdl_dVb = T1 * dgche_dVb - T2 * dRds_dVb;

    /* Calculate degradation factor due to pocket implant */

    if (pParam.BSIM4fprout <= 0.0)
    {   FP = 1.0;
        dFP_dVg = 0.0;
    }
    else
    {   T9 = pParam.BSIM4fprout * std::sqrt(Leff) / Vgst2Vtm;
            FP = 1.0 / (1.0 + T9);
            dFP_dVg = FP * FP * T9 / Vgst2Vtm;
    }

    /* Calculate VACLM */
    T8 = pParam.BSIM4pvag / EsatL;
    T9 = T8 * Vgsteff;
    if (T9 > -0.9)
    {   PvagTerm = 1.0 + T9;
        dPvagTerm_dVg = T8 * (1.0 - Vgsteff * dEsatL_dVg / EsatL);
        dPvagTerm_dVb = -T9 * dEsatL_dVb / EsatL;
        dPvagTerm_dVd = -T9 * dEsatL_dVd / EsatL;
    }
    else
    {   T4 = 1.0 / (17.0 + 20.0 * T9);
        PvagTerm = (0.8 + T9) * T4;
        T4 *= T4;
        dPvagTerm_dVg = T8 * (1.0 - Vgsteff * dEsatL_dVg / EsatL) * T4;
        T9 *= T4 / EsatL;
        dPvagTerm_dVb = -T9 * dEsatL_dVb;
        dPvagTerm_dVd = -T9 * dEsatL_dVd;
    }

    if ((pParam.BSIM4pclm > MIN_EXP) && (diffVds > 1.0e-10))
    {   T0 = 1.0 + Rds * Idl;
        dT0_dVg = dRds_dVg * Idl + Rds * dIdl_dVg;
        dT0_dVd = Rds * dIdl_dVd;
        dT0_dVb = dRds_dVb * Idl + Rds * dIdl_dVb;

        T2 = Vdsat / Esat;
        T1 = Leff + T2;
        dT1_dVg = (dVdsat_dVg - T2 * dEsatL_dVg / Leff) / Esat;
        dT1_dVd = (dVdsat_dVd - T2 * dEsatL_dVd / Leff) / Esat;
        dT1_dVb = (dVdsat_dVb - T2 * dEsatL_dVb / Leff) / Esat;

        Cclm = FP * PvagTerm * T0 * T1 / (pParam.BSIM4pclm * pParam.BSIM4litl);
        dCclm_dVg = Cclm * (dFP_dVg / FP + dPvagTerm_dVg / PvagTerm
                + dT0_dVg / T0 + dT1_dVg / T1);
        dCclm_dVb = Cclm * (dPvagTerm_dVb / PvagTerm + dT0_dVb / T0
                + dT1_dVb / T1);
        dCclm_dVd = Cclm * (dPvagTerm_dVd / PvagTerm + dT0_dVd / T0
                + dT1_dVd / T1);
        VACLM = Cclm * diffVds;

        dVACLM_dVg = dCclm_dVg * diffVds - dVdseff_dVg * Cclm;
        dVACLM_dVb = dCclm_dVb * diffVds - dVdseff_dVb * Cclm;
        dVACLM_dVd = dCclm_dVd * diffVds + (1.0 - dVdseff_dVd) * Cclm;
    }
    else
    {   VACLM = Cclm = MAX_EXP;
        dVACLM_dVd = dVACLM_dVg = dVACLM_dVb = 0.0;
        dCclm_dVd = dCclm_dVg = dCclm_dVb = 0.0;
    }

    /* Calculate VADIBL */
    if (pParam.BSIM4thetaRout > MIN_EXP)
    {   
        T8 = Abulk * Vdsat;
        T0 = Vgst2Vtm * T8;
        dT0_dVg = Vgst2Vtm * Abulk * dVdsat_dVg + T8
            + Vgst2Vtm * Vdsat * dAbulk_dVg;
        dT0_dVb = Vgst2Vtm * (dAbulk_dVb * Vdsat + Abulk * dVdsat_dVb);
        dT0_dVd = Vgst2Vtm * Abulk * dVdsat_dVd;

        T1 = Vgst2Vtm + T8;
        dT1_dVg = 1.0 + Abulk * dVdsat_dVg + Vdsat * dAbulk_dVg;
        dT1_dVb = Abulk * dVdsat_dVb + dAbulk_dVb * Vdsat;
        dT1_dVd = Abulk * dVdsat_dVd;

        T9 = T1 * T1;
        T2 = pParam.BSIM4thetaRout;
        VADIBL = (Vgst2Vtm - T0 / T1) / T2;
        dVADIBL_dVg = (1.0 - dT0_dVg / T1 + T0 * dT1_dVg / T9) / T2;
        dVADIBL_dVb = (-dT0_dVb / T1 + T0 * dT1_dVb / T9) / T2;
        dVADIBL_dVd = (-dT0_dVd / T1 + T0 * dT1_dVd / T9) / T2;

        T7 = pParam.BSIM4pdiblb * Vbseff;
        if (T7 >= -0.9)
        {   T3 = 1.0 / (1.0 + T7);
            VADIBL *= T3;
            dVADIBL_dVg *= T3;
            dVADIBL_dVb = (dVADIBL_dVb - VADIBL * pParam.BSIM4pdiblb) * T3;
            dVADIBL_dVd *= T3;
        }
        else
        {   T4 = 1.0 / (0.8 + T7);
            T3 = (17.0 + 20.0 * T7) * T4;
            dVADIBL_dVg *= T3;
            dVADIBL_dVb = dVADIBL_dVb * T3
                - VADIBL * pParam.BSIM4pdiblb * T4 * T4;
            dVADIBL_dVd *= T3;
            VADIBL *= T3;
        }

        dVADIBL_dVg = dVADIBL_dVg * PvagTerm + VADIBL * dPvagTerm_dVg;
        dVADIBL_dVb = dVADIBL_dVb * PvagTerm + VADIBL * dPvagTerm_dVb;
        dVADIBL_dVd = dVADIBL_dVd * PvagTerm + VADIBL * dPvagTerm_dVd;
        VADIBL *= PvagTerm;
    }
    else
    {   VADIBL = MAX_EXP;
        dVADIBL_dVd = dVADIBL_dVg = dVADIBL_dVb = 0.0;
    }

    /* Calculate Va */
    Va = Vasat + VACLM;
    dVa_dVg = dVasat_dVg + dVACLM_dVg;
    dVa_dVb = dVasat_dVb + dVACLM_dVb;
    dVa_dVd = dVasat_dVd + dVACLM_dVd;

    /* Calculate VADITS */
    T0 = pParam.BSIM4pditsd * Vds;
    if (T0 > EXP_THRESHOLD)
    {   T1 = MAX_EXP;
        dT1_dVd = 0;
    }
    else
    {   T1 = std::exp(T0);
        dT1_dVd = T1 * pParam.BSIM4pditsd;
    }

    if (pParam.BSIM4pdits > MIN_EXP)
    {   T2 = 1.0 + model.BSIM4pditsl * Leff;
        VADITS = (1.0 + T2 * T1) / pParam.BSIM4pdits;
        dVADITS_dVg = VADITS * dFP_dVg;
        dVADITS_dVd = FP * T2 * dT1_dVd / pParam.BSIM4pdits;
        VADITS *= FP;
    }
    else
    {   VADITS = MAX_EXP;
        dVADITS_dVg = dVADITS_dVd = 0;
    }

    /* Calculate VASCBE */
    if ((pParam.BSIM4pscbe2 > 0.0)&&(pParam.BSIM4pscbe1>=0.0)) /*4.6.2*/
    {   
        if (diffVds > pParam.BSIM4pscbe1 * pParam.BSIM4litl / EXP_THRESHOLD)
        {   T0 =  pParam.BSIM4pscbe1 * pParam.BSIM4litl / diffVds;
            VASCBE = Leff * std::exp(T0) / pParam.BSIM4pscbe2;
            T1 = T0 * VASCBE / diffVds;
            dVASCBE_dVg = T1 * dVdseff_dVg;
            dVASCBE_dVd = -T1 * (1.0 - dVdseff_dVd);
            dVASCBE_dVb = T1 * dVdseff_dVb;
        }
        else
        {   VASCBE = MAX_EXP * Leff/pParam.BSIM4pscbe2;
            dVASCBE_dVg = dVASCBE_dVd = dVASCBE_dVb = 0.0;
        }
    }
    else
    {   VASCBE = MAX_EXP;
        dVASCBE_dVg = dVASCBE_dVd = dVASCBE_dVb = 0.0;
    }

    /* Add DIBL to Ids */
    T9 = diffVds / VADIBL;
    T0 = 1.0 + T9;
    Idsa = Idl * T0;
    dIdsa_dVg = T0 * dIdl_dVg - Idl * (dVdseff_dVg + T9 * dVADIBL_dVg) / VADIBL;
    dIdsa_dVd = T0 * dIdl_dVd + Idl
            * (1.0 - dVdseff_dVd - T9 * dVADIBL_dVd) / VADIBL;
    dIdsa_dVb = T0 * dIdl_dVb - Idl * (dVdseff_dVb + T9 * dVADIBL_dVb) / VADIBL;

    /* Add DITS to Ids */
    T9 = diffVds / VADITS;
    T0 = 1.0 + T9;
    dIdsa_dVg = T0 * dIdsa_dVg - Idsa * (dVdseff_dVg + T9 * dVADITS_dVg) / VADITS;
    dIdsa_dVd = T0 * dIdsa_dVd + Idsa
    * (1.0 - dVdseff_dVd - T9 * dVADITS_dVd) / VADITS;
    dIdsa_dVb = T0 * dIdsa_dVb - Idsa * dVdseff_dVb / VADITS;
    Idsa *= T0;

    /* Add CLM to Ids */
    T0 = std::log(Va / Vasat);
    dT0_dVg = dVa_dVg / Va - dVasat_dVg / Vasat;
    dT0_dVb = dVa_dVb / Va - dVasat_dVb / Vasat;
    dT0_dVd = dVa_dVd / Va - dVasat_dVd / Vasat;
    T1 = T0 / Cclm;
    T9 = 1.0 + T1;
    dT9_dVg = (dT0_dVg - T1 * dCclm_dVg) / Cclm;
    dT9_dVb = (dT0_dVb - T1 * dCclm_dVb) / Cclm;
    dT9_dVd = (dT0_dVd - T1 * dCclm_dVd) / Cclm;

    dIdsa_dVg = dIdsa_dVg * T9 + Idsa * dT9_dVg;
    dIdsa_dVb = dIdsa_dVb * T9 + Idsa * dT9_dVb;
    dIdsa_dVd = dIdsa_dVd * T9 + Idsa * dT9_dVd;
    Idsa *= T9;

    /* Substrate current begins */
    tmp = pParam.BSIM4alpha0 + pParam.BSIM4alpha1 * Leff;
    if ((tmp <= 0.0) || (pParam.BSIM4beta0 <= 0.0))
    {   Isub = Gbd = Gbb = Gbg = 0.0;
    }
    else
    {   
        T2 = tmp / Leff;
        if (diffVds > pParam.BSIM4beta0 / EXP_THRESHOLD)
        {   
            T0 = -pParam.BSIM4beta0 / diffVds;
            T1 = T2 * diffVds * std::exp(T0);
            T3 = T1 / diffVds * (T0 - 1.0);
            dT1_dVg = T3 * dVdseff_dVg;
            dT1_dVd = T3 * (dVdseff_dVd - 1.0);
            dT1_dVb = T3 * dVdseff_dVb;
        }
        else
        {   
            T3 = T2 * MIN_EXP;
            T1 = T3 * diffVds;
            dT1_dVg = -T3 * dVdseff_dVg;
            dT1_dVd = T3 * (1.0 - dVdseff_dVd);
            dT1_dVb = -T3 * dVdseff_dVb;
        }
        T4 = Idsa * Vdseff;
        Isub = T1 * T4;
        Gbg = T1 * (dIdsa_dVg * Vdseff + Idsa * dVdseff_dVg)
            + T4 * dT1_dVg;
        Gbd = T1 * (dIdsa_dVd * Vdseff + Idsa * dVdseff_dVd)
            + T4 * dT1_dVd;
        Gbb = T1 * (dIdsa_dVb * Vdseff + Idsa * dVdseff_dVb)
            + T4 * dT1_dVb;

        Gbd += Gbg * dVgsteff_dVd;
        Gbb += Gbg * dVgsteff_dVb;
        Gbg *= dVgsteff_dVg;
        Gbb *= dVbseff_dVb;
    }

    instance.BSIM4csub = Isub;
    instance.BSIM4gbbs = Gbb;
    instance.BSIM4gbgs = Gbg;
    instance.BSIM4gbds = Gbd;

    /* Add SCBE to Ids */
    T9 = diffVds / VASCBE;
    T0 = 1.0 + T9;
    Ids = Idsa * T0;

    Gm = T0 * dIdsa_dVg - Idsa
    * (dVdseff_dVg + T9 * dVASCBE_dVg) / VASCBE;
    Gds = T0 * dIdsa_dVd + Idsa
    * (1.0 - dVdseff_dVd - T9 * dVASCBE_dVd) / VASCBE;
    Gmb = T0 * dIdsa_dVb - Idsa
    * (dVdseff_dVb + T9 * dVASCBE_dVb) / VASCBE;


    tmp1 = Gds + Gm * dVgsteff_dVd;
    tmp2 = Gmb + Gm * dVgsteff_dVb;
    tmp3 = Gm;

    Gm = (Ids * dVdseff_dVg + Vdseff * tmp3) * dVgsteff_dVg;
    Gds = Ids * (dVdseff_dVd + dVdseff_dVg * dVgsteff_dVd)
        + Vdseff * tmp1;
    Gmb = (Ids * (dVdseff_dVb + dVdseff_dVg * dVgsteff_dVb)
        + Vdseff * tmp2) * dVbseff_dVb;

    cdrain = Ids * Vdseff;

    /* Source End Velocity Limit  */
    if((model.BSIM4vtlGiven) && (model.BSIM4vtl > 0.0) ) {
        T12 = 1.0 / Leff / CoxeffWovL;
        T11 = T12 / Vgsteff;
        T10 = -T11 / Vgsteff;
        vs = cdrain * T11; /* vs */
        dvs_dVg = Gm * T11 + cdrain * T10 * dVgsteff_dVg;
        dvs_dVd = Gds * T11 + cdrain * T10 * dVgsteff_dVd;
        dvs_dVb = Gmb * T11 + cdrain * T10 * dVgsteff_dVb;
        T0 = 2 * MM;
        T1 = vs / (pParam.BSIM4vtl * pParam.BSIM4tfactor);
        if(T1 > 0.0)
        { T2 = 1.0 + std::exp(T0 * std::log(T1));
        T3 = (T2 - 1.0) * T0 / vs;
        Fsevl = 1.0 / std::exp(log(T2)/ T0);
        dT2_dVg = T3 * dvs_dVg;
        dT2_dVd = T3 * dvs_dVd;
        dT2_dVb = T3 * dvs_dVb;
        T4 = -1.0 / T0 * Fsevl / T2;
        dFsevl_dVg = T4 * dT2_dVg;
        dFsevl_dVd = T4 * dT2_dVd;
        dFsevl_dVb = T4 * dT2_dVb;
        } else {
        Fsevl = 1.0;
        dFsevl_dVg = 0.0;
        dFsevl_dVd = 0.0;
        dFsevl_dVb = 0.0;
        }
        Gm *=Fsevl;
        Gm += cdrain * dFsevl_dVg;
        Gmb *=Fsevl;
        Gmb += cdrain * dFsevl_dVb;
        Gds *=Fsevl;
        Gds += cdrain * dFsevl_dVd;

        cdrain *= Fsevl;
    }

    instance.BSIM4gds = Gds;
    instance.BSIM4gm = Gm;
    instance.BSIM4gmbs = Gmb;
    instance.BSIM4IdovVds = Ids;
    if( instance.BSIM4IdovVds <= model.BSIM4idovvdsc) instance.BSIM4IdovVds = model.BSIM4idovvdsc; 

    /* Calculate Rg */
    if ((instance.BSIM4rgateMod > 1) ||
        (instance.BSIM4trnqsMod != 0) || (instance.BSIM4acnqsMod != 0))
    {   
        T9 = pParam.BSIM4xrcrg2 * model.BSIM4vtm;
        T0 = T9 * beta;
        dT0_dVd = (dbeta_dVd + dbeta_dVg * dVgsteff_dVd) * T9;
        dT0_dVb = (dbeta_dVb + dbeta_dVg * dVgsteff_dVb) * T9;
        dT0_dVg = dbeta_dVg * T9;

        instance.BSIM4gcrg = pParam.BSIM4xrcrg1 * ( T0 + Ids);
        instance.BSIM4gcrgd = pParam.BSIM4xrcrg1 * (dT0_dVd + tmp1);
        instance.BSIM4gcrgb = pParam.BSIM4xrcrg1 * (dT0_dVb + tmp2)
            * dVbseff_dVb;
        instance.BSIM4gcrgg = pParam.BSIM4xrcrg1 * (dT0_dVg + tmp3)
            * dVgsteff_dVg;

        if (instance.BSIM4nf != 1.0)
        {   
            instance.BSIM4gcrg *= instance.BSIM4nf;
            instance.BSIM4gcrgg *= instance.BSIM4nf;
            instance.BSIM4gcrgd *= instance.BSIM4nf;
            instance.BSIM4gcrgb *= instance.BSIM4nf;
        }

        if (instance.BSIM4rgateMod == 2)
        {   
            T10 = instance.BSIM4grgeltd * instance.BSIM4grgeltd;
            T11 = instance.BSIM4grgeltd + instance.BSIM4gcrg;
            instance.BSIM4gcrg = instance.BSIM4grgeltd * instance.BSIM4gcrg / T11;
            T12 = T10 / T11 / T11;
            instance.BSIM4gcrgg *= T12;
            instance.BSIM4gcrgd *= T12;
            instance.BSIM4gcrgb *= T12;
        }
        instance.BSIM4gcrgs = -(instance.BSIM4gcrgg + instance.BSIM4gcrgd + instance.BSIM4gcrgb);
    }


    /* Calculate bias-dependent external S/D resistance */
    if (model.BSIM4rdsMod)
    {   
        /* Rs(V) */
        T0 = vgs - pParam.BSIM4vfbsd;
        T1 = std::sqrt(T0 * T0 + 1.0e-4);
        vgs_eff = 0.5 * (T0 + T1);
        dvgs_eff_dvg = vgs_eff / T1;

        T0 = 1.0 + pParam.BSIM4prwg * vgs_eff;
        dT0_dvg = -pParam.BSIM4prwg / T0 / T0 * dvgs_eff_dvg;
        T1 = -pParam.BSIM4prwb * vbs;
        dT1_dvb = -pParam.BSIM4prwb;

        T2 = 1.0 / T0 + T1;
        T3 = T2 + std::sqrt(T2 * T2 + 0.01);
        dT3_dvg = T3 / (T3 - T2);
        dT3_dvb = dT3_dvg * dT1_dvb;
        dT3_dvg *= dT0_dvg;

        T4 = pParam.BSIM4rs0 * 0.5;
        Rs = pParam.BSIM4rswmin + T3 * T4;
        dRs_dvg = T4 * dT3_dvg;
        dRs_dvb = T4 * dT3_dvb;

        T0 = 1.0 + instance.BSIM4sourceConductance * Rs;
        instance.BSIM4gstot = instance.BSIM4sourceConductance / T0;
        T0 = -instance.BSIM4gstot * instance.BSIM4gstot;
        dgstot_dvd = 0.0; /* place holder */
        dgstot_dvg = T0 * dRs_dvg;
        dgstot_dvb = T0 * dRs_dvb;
        dgstot_dvs = -(dgstot_dvg + dgstot_dvb + dgstot_dvd);

        /* Rd(V) */
        T0 = vgd - pParam.BSIM4vfbsd;
        T1 = std::sqrt(T0 * T0 + 1.0e-4);
        vgd_eff = 0.5 * (T0 + T1);
        dvgd_eff_dvg = vgd_eff / T1;

        T0 = 1.0 + pParam.BSIM4prwg * vgd_eff;
        dT0_dvg = -pParam.BSIM4prwg / T0 / T0 * dvgd_eff_dvg;
        T1 = -pParam.BSIM4prwb * vbd;
        dT1_dvb = -pParam.BSIM4prwb;

        T2 = 1.0 / T0 + T1;
        T3 = T2 + std::sqrt(T2 * T2 + 0.01);
        dT3_dvg = T3 / (T3 - T2);
        dT3_dvb = dT3_dvg * dT1_dvb;
        dT3_dvg *= dT0_dvg;

        T4 = pParam.BSIM4rd0 * 0.5;
        Rd = pParam.BSIM4rdwmin + T3 * T4;
        dRd_dvg = T4 * dT3_dvg;
        dRd_dvb = T4 * dT3_dvb;

        T0 = 1.0 + instance.BSIM4drainConductance * Rd;
        instance.BSIM4gdtot = instance.BSIM4drainConductance / T0;
        T0 = -instance.BSIM4gdtot * instance.BSIM4gdtot;
        dgdtot_dvs = 0.0;
        dgdtot_dvg = T0 * dRd_dvg;
        dgdtot_dvb = T0 * dRd_dvb;
        dgdtot_dvd = -(dgdtot_dvg + dgdtot_dvb + dgdtot_dvs);

        instance.BSIM4gstotd = vses * dgstot_dvd;
        instance.BSIM4gstotg = vses * dgstot_dvg;
        instance.BSIM4gstots = vses * dgstot_dvs;
        instance.BSIM4gstotb = vses * dgstot_dvb;

        T2 = vdes - vds;
        instance.BSIM4gdtotd = T2 * dgdtot_dvd;
        instance.BSIM4gdtotg = T2 * dgdtot_dvg;
        instance.BSIM4gdtots = T2 * dgdtot_dvs;
        instance.BSIM4gdtotb = T2 * dgdtot_dvb;
    }
    else /* WDLiu: for bypass */
    {   
        instance.BSIM4gstot = instance.BSIM4gstotd = instance.BSIM4gstotg = 0.0;
        instance.BSIM4gstots = instance.BSIM4gstotb = 0.0;
        instance.BSIM4gdtot = instance.BSIM4gdtotd = instance.BSIM4gdtotg = 0.0;
        instance.BSIM4gdtots = instance.BSIM4gdtotb = 0.0;
    }

    /* GIDL/GISL Models */

    if(model.BSIM4mtrlMod == 0)
    T0 = 3.0 * toxe;
    else
    T0 = model.BSIM4epsrsub * toxe / epsrox;

    /* Calculate GIDL current */

    vgs_eff = instance.BSIM4vgs_eff;
    dvgs_eff_dvg = instance.BSIM4dvgs_eff_dvg;
    vgd_eff = instance.BSIM4vgd_eff;
    dvgd_eff_dvg = instance.BSIM4dvgd_eff_dvg;

    if (model.BSIM4gidlMod==0){

        if(model.BSIM4mtrlMod ==0)
            T1 = (vds - vgs_eff - pParam.BSIM4egidl ) / T0;
        else
            T1 = (vds - vgs_eff - pParam.BSIM4egidl + pParam.BSIM4vfbsd) / T0;

        if ((pParam.BSIM4agidl <= 0.0) || (pParam.BSIM4bgidl <= 0.0)
            || (T1 <= 0.0) || (pParam.BSIM4cgidl <= 0.0) || (vbd > 0.0))
            Igidl = Ggidld = Ggidlg = Ggidlb = 0.0;
        else {
            dT1_dVd = 1.0 / T0;
            dT1_dVg = -dvgs_eff_dvg * dT1_dVd;
            T2 = pParam.BSIM4bgidl / T1;
            if (T2 < 100.0)
            {   Igidl = pParam.BSIM4agidl * pParam.BSIM4weffCJ * T1 * std::exp(-T2);
                T3 = Igidl * (1.0 + T2) / T1;
                Ggidld = T3 * dT1_dVd;
                Ggidlg = T3 * dT1_dVg;
            }
            else
            {   Igidl = pParam.BSIM4agidl * pParam.BSIM4weffCJ * 3.720075976e-44;
                Ggidld = Igidl * dT1_dVd;
                Ggidlg = Igidl * dT1_dVg;
                Igidl *= T1;
            }

            T4 = vbd * vbd;
            T5 = -vbd * T4;
            T6 = pParam.BSIM4cgidl + T5;
            T7 = T5 / T6;
            T8 = 3.0 * pParam.BSIM4cgidl * T4 / T6 / T6;
            Ggidld = Ggidld * T7 + Igidl * T8;
            Ggidlg = Ggidlg * T7;
            Ggidlb = -Igidl * T8;
            Igidl *= T7;
        }
        instance.BSIM4Igidl = Igidl;
        instance.BSIM4ggidld = Ggidld;
        instance.BSIM4ggidlg = Ggidlg;
        instance.BSIM4ggidlb = Ggidlb;
        /* Calculate GISL current  */

        if(model.BSIM4mtrlMod ==0)
            T1 = (-vds - vgd_eff - pParam.BSIM4egisl ) / T0;
        else
            T1 = (-vds - vgd_eff - pParam.BSIM4egisl + pParam.BSIM4vfbsd ) / T0;

        if ((pParam.BSIM4agisl <= 0.0) || (pParam.BSIM4bgisl <= 0.0)
            || (T1 <= 0.0) || (pParam.BSIM4cgisl <= 0.0) || (vbs > 0.0))
            Igisl = Ggisls = Ggislg = Ggislb = 0.0;
        else {
            dT1_dVd = 1.0 / T0;
            dT1_dVg = -dvgd_eff_dvg * dT1_dVd;
            T2 = pParam.BSIM4bgisl / T1;
            if (T2 < 100.0)
            {   Igisl = pParam.BSIM4agisl * pParam.BSIM4weffCJ * T1 * std::exp(-T2);
                T3 = Igisl * (1.0 + T2) / T1;
                Ggisls = T3 * dT1_dVd;
                Ggislg = T3 * dT1_dVg;
            }
            else
            {   Igisl = pParam.BSIM4agisl * pParam.BSIM4weffCJ * 3.720075976e-44;
                Ggisls = Igisl * dT1_dVd;
                Ggislg = Igisl * dT1_dVg;
                Igisl *= T1;
            }

            T4 = vbs * vbs;
            T5 = -vbs * T4;
            T6 = pParam.BSIM4cgisl + T5;
            T7 = T5 / T6;
            T8 = 3.0 * pParam.BSIM4cgisl * T4 / T6 / T6;
            Ggisls = Ggisls * T7 + Igisl * T8;
            Ggislg = Ggislg * T7;
            Ggislb = -Igisl * T8;
            Igisl *= T7;
        }
        instance.BSIM4Igisl = Igisl;
        instance.BSIM4ggisls = Ggisls;
        instance.BSIM4ggislg = Ggislg;
        instance.BSIM4ggislb = Ggislb;
    }
    else{
	/* v4.7 New Gidl/GISL model */

		/* GISL */
                    if (model.BSIM4mtrlMod==0)
                       T1 = (-vds - pParam.BSIM4rgisl * vgd_eff - pParam.BSIM4egisl) / T0;
                    else
                       T1 = (-vds - pParam.BSIM4rgisl * vgd_eff - pParam.BSIM4egisl + pParam.BSIM4vfbsd) / T0;

            if ((pParam.BSIM4agisl <= 0.0) ||
                            (pParam.BSIM4bgisl <= 0.0) || (T1 <= 0.0) ||
                            (pParam.BSIM4cgisl < 0.0)  )
                        Igisl = Ggisls = Ggislg = Ggislb = 0.0;
                    else
                    {
                        dT1_dVd = 1 / T0;
                        dT1_dVg = - pParam.BSIM4rgisl * dT1_dVd * dvgd_eff_dvg;
                        T2 = pParam.BSIM4bgisl / T1;
                        if (T2 < EXPL_THRESHOLD)
                        {
                            Igisl = pParam.BSIM4weffCJ * pParam.BSIM4agisl * T1 * std::exp(-T2);
                            T3 = Igisl / T1 * (T2 + 1);
                            Ggisls = T3 * dT1_dVd;
                            Ggislg = T3 * dT1_dVg;
                        }
            else
                        {
                            T3 = pParam.BSIM4weffCJ * pParam.BSIM4agisl * MIN_EXPL;
                            Igisl = T3 * T1 ;
                            Ggisls  = T3 * dT1_dVd;
                            Ggislg  = T3 * dT1_dVg;

                        }
                        T4 = vbs - pParam.BSIM4fgisl;
		/*--chetan dabhi solution for clamping T4-*/
			if(T4 > model.BSIM4gidlclamp)					
			    T4=model.BSIM4gidlclamp; 

            if (T4==0)
                            T5 = EXPL_THRESHOLD;
                        else
                            T5 = pParam.BSIM4kgisl / T4;
                        if (T5<EXPL_THRESHOLD)
                        {T6 = std::exp(T5);
                            Ggislb = -Igisl * T6 * T5 / T4;
                        }
                        else
                        {T6 = MAX_EXPL;
                            Ggislb=0.0;
                        }
                        Ggisls*=T6;
                        Ggislg*=T6;
                        Igisl*=T6;

                    }
                instance.BSIM4Igisl = Igisl;
                instance.BSIM4ggisls = Ggisls;
                instance.BSIM4ggislg = Ggislg;
                instance.BSIM4ggislb = Ggislb;
                    /* End of GISL */

                /* GIDL */
                    if (model.BSIM4mtrlMod==0)
                        T1 = (vds - pParam.BSIM4rgidl * vgs_eff - pParam.BSIM4egidl) /  T0;
                    else
                        T1 = (vds - pParam.BSIM4rgidl * vgs_eff - pParam.BSIM4egidl + pParam.BSIM4vfbsd) / T0;



                    if ((pParam.BSIM4agidl <= 0.0) ||
                            (pParam.BSIM4bgidl <= 0.0) || (T1 <= 0.0) ||
                            (pParam.BSIM4cgidl < 0.0)  )
                        Igidl = Ggidld = Ggidlg = Ggidlb = 0.0;
                    else
                    {
                        dT1_dVd = 1 / T0;
                        dT1_dVg = - pParam.BSIM4rgidl * dT1_dVd * dvgs_eff_dvg;
                        T2 = pParam.BSIM4bgidl / T1;
                        if (T2 < EXPL_THRESHOLD)
                        {
                            Igidl = pParam.BSIM4weffCJ * pParam.BSIM4agidl * T1 * std::exp(-T2);
                            T3 = Igidl / T1 * (T2 + 1);
                            Ggidld = T3 * dT1_dVd;
                            Ggidlg = T3 * dT1_dVg;

                        } else
                        {
                            T3 = pParam.BSIM4weffCJ * pParam.BSIM4agidl * MIN_EXPL;
                            Igidl = T3 * T1 ;
                            Ggidld  = T3 * dT1_dVd;
                            Ggidlg  = T3 * dT1_dVg;
                        }
                        T4 = vbd - pParam.BSIM4fgidl;
		/*--chetan dabhi solution for clamping T4-*/
			if(T4 > model.BSIM4gidlclamp)					
			    T4=model.BSIM4gidlclamp; 
                        if (T4==0)
                            T5 = EXPL_THRESHOLD;
                        else
                            T5 = pParam.BSIM4kgidl / T4;
                        if (T5<EXPL_THRESHOLD)
                        {T6 = std::exp(T5);
                            Ggidlb = -Igidl * T6 * T5 / T4;
                        }
                        else
                        {T6 = MAX_EXPL;
                            Ggidlb=0.0;
                        }
                        Ggidld *= T6;
                        Ggidlg *= T6;
                        Igidl *= T6;
                    }
        instance.BSIM4Igidl = Igidl;
                instance.BSIM4ggidld = Ggidld;
            instance.BSIM4ggidlg = Ggidlg;
            instance.BSIM4ggidlb = Ggidlb;
                /* End of New GIDL */
        }
        /*End of Gidl*/



          /* Calculate gate tunneling current */
          if ((model.BSIM4igcMod != 0) || (model.BSIM4igbMod != 0))
          {   Vfb = instance.BSIM4vfbzb;
              V3 = Vfb - Vgs_eff + Vbseff - DELTA_3;
              if (Vfb <= 0.0)
                  T0 = std::sqrt(V3 * V3 - 4.0 * DELTA_3 * Vfb);
              else
                  T0 = std::sqrt(V3 * V3 + 4.0 * DELTA_3 * Vfb);
              T1 = 0.5 * (1.0 + V3 / T0);
              Vfbeff = Vfb - 0.5 * (V3 + T0);
              dVfbeff_dVg = T1 * dVgs_eff_dVg;
              dVfbeff_dVb = -T1; /* WDLiu: -No surprise? No. -Good! */

              Voxacc = Vfb - Vfbeff;
              dVoxacc_dVg = -dVfbeff_dVg;
              dVoxacc_dVb = -dVfbeff_dVb;
              if (Voxacc < 0.0) /* WDLiu: Avoiding numerical instability. */
                  Voxacc = dVoxacc_dVg = dVoxacc_dVb = 0.0;

              T0 = 0.5 * pParam.BSIM4k1ox;
              T3 = Vgs_eff - Vfbeff - Vbseff - Vgsteff;
              if (pParam.BSIM4k1ox == 0.0)
                  Voxdepinv = dVoxdepinv_dVg = dVoxdepinv_dVd
                = dVoxdepinv_dVb = 0.0;
              else if (T3 < 0.0)
              {   Voxdepinv = -T3;
                  dVoxdepinv_dVg = -dVgs_eff_dVg + dVfbeff_dVg
                                 + dVgsteff_dVg;
                  dVoxdepinv_dVd = dVgsteff_dVd;
                  dVoxdepinv_dVb = dVfbeff_dVb + 1.0 + dVgsteff_dVb;
          }
              else
              {   T1 = std::sqrt(T0 * T0 + T3);
                  T2 = T0 / T1;
                  Voxdepinv = pParam.BSIM4k1ox * (T1 - T0);
                  dVoxdepinv_dVg = T2 * (dVgs_eff_dVg - dVfbeff_dVg
                 - dVgsteff_dVg);
                  dVoxdepinv_dVd = -T2 * dVgsteff_dVd;
                  dVoxdepinv_dVb = -T2 * (dVfbeff_dVb + 1.0 + dVgsteff_dVb);
              }

              Voxdepinv += Vgsteff;
              dVoxdepinv_dVg += dVgsteff_dVg;
              dVoxdepinv_dVd += dVgsteff_dVd;
              dVoxdepinv_dVb += dVgsteff_dVb;
          }

          if(model.BSIM4tempMod < 2)
            tmp = Vtm;
          else /* model.BSIM4tempMod = 2 , 3*/
            tmp = Vtm0;
          if (model.BSIM4igcMod)
          {   T0 = tmp * pParam.BSIM4nigc;
          if(model.BSIM4igcMod == 1) {
                VxNVt = (Vgs_eff - model.BSIM4type * instance.BSIM4vth0) / T0;
                if (VxNVt > EXP_THRESHOLD)
                {   Vaux = Vgs_eff - model.BSIM4type * instance.BSIM4vth0;
                    dVaux_dVg = dVgs_eff_dVg;
            dVaux_dVd = 0.0;
                    dVaux_dVb = 0.0;
                }
          } else if (model.BSIM4igcMod == 2) {
                VxNVt = (Vgs_eff - instance.BSIM4von) / T0;
                if (VxNVt > EXP_THRESHOLD)
                {   Vaux = Vgs_eff - instance.BSIM4von;
                    dVaux_dVg = dVgs_eff_dVg;
                    dVaux_dVd = -dVth_dVd;
                    dVaux_dVb = -dVth_dVb;
                }
              }
              if (VxNVt < -EXP_THRESHOLD)
              {   Vaux = T0 * std::log(1.0 + MIN_EXP);
                  dVaux_dVg = dVaux_dVd = dVaux_dVb = 0.0;
              }
              else if ((VxNVt >= -EXP_THRESHOLD) && (VxNVt <= EXP_THRESHOLD))
              {   ExpVxNVt = std::exp(VxNVt);
                  Vaux = T0 * std::log(1.0 + ExpVxNVt);
                  dVaux_dVg = ExpVxNVt / (1.0 + ExpVxNVt);
          if(model.BSIM4igcMod == 1) {
            dVaux_dVd = 0.0;
                    dVaux_dVb = 0.0;
                  } else if (model.BSIM4igcMod == 2) {
                        dVaux_dVd = -dVaux_dVg* dVth_dVd;  /* Synopsys 08/30/2013 modify */
                        dVaux_dVb = -dVaux_dVg* dVth_dVb;  /* Synopsys 08/30/2013 modify */
          }
          dVaux_dVg *= dVgs_eff_dVg;
              }

              T2 = Vgs_eff * Vaux;
              dT2_dVg = dVgs_eff_dVg * Vaux + Vgs_eff * dVaux_dVg;
              dT2_dVd = Vgs_eff * dVaux_dVd;
              dT2_dVb = Vgs_eff * dVaux_dVb;

              T11 = pParam.BSIM4Aechvb;
              T12 = pParam.BSIM4Bechvb;
              T3 = pParam.BSIM4aigc * pParam.BSIM4cigc
                 - pParam.BSIM4bigc;
              T4 = pParam.BSIM4bigc * pParam.BSIM4cigc;
              T5 = T12 * (pParam.BSIM4aigc + T3 * Voxdepinv
                 - T4 * Voxdepinv * Voxdepinv);

              if (T5 > EXP_THRESHOLD)
              {   T6 = MAX_EXP;
                  dT6_dVg = dT6_dVd = dT6_dVb = 0.0;
              }
              else if (T5 < -EXP_THRESHOLD)
              {   T6 = MIN_EXP;
                  dT6_dVg = dT6_dVd = dT6_dVb = 0.0;
              }
              else
              {   T6 = std::exp(T5);
                  dT6_dVg = T6 * T12 * (T3 - 2.0 * T4 * Voxdepinv);
                  dT6_dVd = dT6_dVg * dVoxdepinv_dVd;
                  dT6_dVb = dT6_dVg * dVoxdepinv_dVb;
                  dT6_dVg *= dVoxdepinv_dVg;
              }

              Igc = T11 * T2 * T6;
              dIgc_dVg = T11 * (T2 * dT6_dVg + T6 * dT2_dVg);
              dIgc_dVd = T11 * (T2 * dT6_dVd + T6 * dT2_dVd);
              dIgc_dVb = T11 * (T2 * dT6_dVb + T6 * dT2_dVb);

              if (model.BSIM4pigcdGiven)
              {   Pigcd = pParam.BSIM4pigcd;
                  dPigcd_dVg = dPigcd_dVd = dPigcd_dVb = 0.0;
              }
              else
              {  /* T11 = pParam.BSIM4Bechvb * toxe; v4.7 */
          T11 = -pParam.BSIM4Bechvb;
                  T12 = Vgsteff + 1.0e-20;
                  T13 = T11 / T12 / T12;
                  T14 = -T13 / T12;
                  Pigcd = T13 * (1.0 - 0.5 * Vdseff / T12);
                  dPigcd_dVg = T14 * (2.0 + 0.5 * (dVdseff_dVg
                              - 3.0 * Vdseff / T12));
                  dPigcd_dVd = 0.5 * T14 * dVdseff_dVd;
                  dPigcd_dVb = 0.5 * T14 * dVdseff_dVb;
              }

              T7 = -Pigcd * Vdseff; /* bugfix */
              dT7_dVg = -Vdseff * dPigcd_dVg - Pigcd * dVdseff_dVg;
              dT7_dVd = -Vdseff * dPigcd_dVd - Pigcd * dVdseff_dVd + dT7_dVg * dVgsteff_dVd;
              dT7_dVb = -Vdseff * dPigcd_dVb - Pigcd * dVdseff_dVb + dT7_dVg * dVgsteff_dVb;
              dT7_dVg *= dVgsteff_dVg;
              /*dT7_dVb *= dVbseff_dVb;*/ /* Synopsys, 2013/08/30 */
              T8 = T7 * T7 + 2.0e-4;
              dT8_dVg = 2.0 * T7;
              dT8_dVd = dT8_dVg * dT7_dVd;
              dT8_dVb = dT8_dVg * dT7_dVb;
              dT8_dVg *= dT7_dVg;

              if (T7 > EXP_THRESHOLD)
              {   T9 = MAX_EXP;
                  dT9_dVg = dT9_dVd = dT9_dVb = 0.0;
              }
              else if (T7 < -EXP_THRESHOLD)
              {   T9 = MIN_EXP;
                  dT9_dVg = dT9_dVd = dT9_dVb = 0.0;
              }
              else
              {   T9 = std::exp(T7);
                  dT9_dVg = T9 * dT7_dVg;
                  dT9_dVd = T9 * dT7_dVd;
                  dT9_dVb = T9 * dT7_dVb;
              }

              T0 = T8 * T8;
              T1 = T9 - 1.0 + 1.0e-4;
              T10 = (T1 - T7) / T8;
              dT10_dVg = (dT9_dVg - dT7_dVg - T10 * dT8_dVg) / T8;
              dT10_dVd = (dT9_dVd - dT7_dVd - T10 * dT8_dVd) / T8;
              dT10_dVb = (dT9_dVb - dT7_dVb - T10 * dT8_dVb) / T8;

              Igcs = Igc * T10;
              dIgcs_dVg = dIgc_dVg * T10 + Igc * dT10_dVg;
              dIgcs_dVd = dIgc_dVd * T10 + Igc * dT10_dVd;
              dIgcs_dVb = dIgc_dVb * T10 + Igc * dT10_dVb;

              T1 = T9 - 1.0 - 1.0e-4;
              T10 = (T7 * T9 - T1) / T8;
              dT10_dVg = (dT7_dVg * T9 + (T7 - 1.0) * dT9_dVg
                       - T10 * dT8_dVg) / T8;
              dT10_dVd = (dT7_dVd * T9 + (T7 - 1.0) * dT9_dVd
                       - T10 * dT8_dVd) / T8;
              dT10_dVb = (dT7_dVb * T9 + (T7 - 1.0) * dT9_dVb
                       - T10 * dT8_dVb) / T8;
              Igcd = Igc * T10;
              dIgcd_dVg = dIgc_dVg * T10 + Igc * dT10_dVg;
              dIgcd_dVd = dIgc_dVd * T10 + Igc * dT10_dVd;
              dIgcd_dVb = dIgc_dVb * T10 + Igc * dT10_dVb;

              instance.BSIM4Igcs = Igcs;
              instance.BSIM4gIgcsg = dIgcs_dVg;
              instance.BSIM4gIgcsd = dIgcs_dVd;
              instance.BSIM4gIgcsb =  dIgcs_dVb * dVbseff_dVb;
              instance.BSIM4Igcd = Igcd;
              instance.BSIM4gIgcdg = dIgcd_dVg;
              instance.BSIM4gIgcdd = dIgcd_dVd;
              instance.BSIM4gIgcdb = dIgcd_dVb * dVbseff_dVb;

              T0 = vgs - (pParam.BSIM4vfbsd + pParam.BSIM4vfbsdoff);
              vgs_eff = std::sqrt(T0 * T0 + 1.0e-4);
              dvgs_eff_dvg = T0 / vgs_eff;

              T2 = vgs * vgs_eff;
              dT2_dVg = vgs * dvgs_eff_dvg + vgs_eff;
              T11 = pParam.BSIM4AechvbEdgeS;
              T12 = pParam.BSIM4BechvbEdge;
              T3 = pParam.BSIM4aigs * pParam.BSIM4cigs
                 - pParam.BSIM4bigs;
              T4 = pParam.BSIM4bigs * pParam.BSIM4cigs;
              T5 = T12 * (pParam.BSIM4aigs + T3 * vgs_eff
                 - T4 * vgs_eff * vgs_eff);
              if (T5 > EXP_THRESHOLD)
              {   T6 = MAX_EXP;
                  dT6_dVg = 0.0;
              }
              else if (T5 < -EXP_THRESHOLD)
              {   T6 = MIN_EXP;
                  dT6_dVg = 0.0;
              }
              else
              {   T6 = std::exp(T5);
                  dT6_dVg = T6 * T12 * (T3 - 2.0 * T4 * vgs_eff)
              * dvgs_eff_dvg;
              }
              Igs = T11 * T2 * T6;
              dIgs_dVg = T11 * (T2 * dT6_dVg + T6 * dT2_dVg);
              dIgs_dVs = -dIgs_dVg;


              T0 = vgd - (pParam.BSIM4vfbsd + pParam.BSIM4vfbsdoff);
              vgd_eff = std::sqrt(T0 * T0 + 1.0e-4);
              dvgd_eff_dvg = T0 / vgd_eff;

              T2 = vgd * vgd_eff;
              dT2_dVg = vgd * dvgd_eff_dvg + vgd_eff;
              T11 = pParam.BSIM4AechvbEdgeD;
              T3 = pParam.BSIM4aigd * pParam.BSIM4cigd
                 - pParam.BSIM4bigd;
              T4 = pParam.BSIM4bigd * pParam.BSIM4cigd;
              T5 = T12 * (pParam.BSIM4aigd + T3 * vgd_eff
                 - T4 * vgd_eff * vgd_eff);
              if (T5 > EXP_THRESHOLD)
              {   T6 = MAX_EXP;
                  dT6_dVg = 0.0;
              }
              else if (T5 < -EXP_THRESHOLD)
              {   T6 = MIN_EXP;
                  dT6_dVg = 0.0;
              }
              else
              {   T6 = std::exp(T5);
                  dT6_dVg = T6 * T12 * (T3 - 2.0 * T4 * vgd_eff)
                          * dvgd_eff_dvg;
              }
              Igd = T11 * T2 * T6;
              dIgd_dVg = T11 * (T2 * dT6_dVg + T6 * dT2_dVg);
              dIgd_dVd = -dIgd_dVg;

              instance.BSIM4Igs = Igs;
              instance.BSIM4gIgsg = dIgs_dVg;
              instance.BSIM4gIgss = dIgs_dVs;
              instance.BSIM4Igd = Igd;
              instance.BSIM4gIgdg = dIgd_dVg;
              instance.BSIM4gIgdd = dIgd_dVd;
          }
          else
          {   instance.BSIM4Igcs = instance.BSIM4gIgcsg = instance.BSIM4gIgcsd
                  = instance.BSIM4gIgcsb = 0.0;
              instance.BSIM4Igcd = instance.BSIM4gIgcdg = instance.BSIM4gIgcdd
                          = instance.BSIM4gIgcdb = 0.0;
              instance.BSIM4Igs = instance.BSIM4gIgsg = instance.BSIM4gIgss = 0.0;
              instance.BSIM4Igd = instance.BSIM4gIgdg = instance.BSIM4gIgdd = 0.0;
          }

          if (model.BSIM4igbMod)
          {   T0 = tmp * pParam.BSIM4nigbacc;
          T1 = -Vgs_eff + Vbseff + Vfb;
              VxNVt = T1 / T0;
              if (VxNVt > EXP_THRESHOLD)
              {   Vaux = T1;
                  dVaux_dVg = -dVgs_eff_dVg;
                  dVaux_dVb = 1.0;
              }
              else if (VxNVt < -EXP_THRESHOLD)
              {   Vaux = T0 * std::log(1.0 + MIN_EXP);
                  dVaux_dVg = dVaux_dVb = 0.0;
              }
              else
              {   ExpVxNVt = std::exp(VxNVt);
                  Vaux = T0 * std::log(1.0 + ExpVxNVt);
                  dVaux_dVb = ExpVxNVt / (1.0 + ExpVxNVt);
                  dVaux_dVg = -dVaux_dVb * dVgs_eff_dVg;
              }

          T2 = (Vgs_eff - Vbseff) * Vaux;
              dT2_dVg = dVgs_eff_dVg * Vaux + (Vgs_eff - Vbseff) * dVaux_dVg;
              dT2_dVb = -Vaux + (Vgs_eff - Vbseff) * dVaux_dVb;

              T11 = 4.97232e-7 * pParam.BSIM4weff
          * pParam.BSIM4leff * pParam.BSIM4ToxRatio;
              T12 = -7.45669e11 * toxe;
          T3 = pParam.BSIM4aigbacc * pParam.BSIM4cigbacc
                 - pParam.BSIM4bigbacc;
              T4 = pParam.BSIM4bigbacc * pParam.BSIM4cigbacc;
          T5 = T12 * (pParam.BSIM4aigbacc + T3 * Voxacc
         - T4 * Voxacc * Voxacc);

              if (T5 > EXP_THRESHOLD)
              {   T6 = MAX_EXP;
                  dT6_dVg = dT6_dVb = 0.0;
              }
              else if (T5 < -EXP_THRESHOLD)
              {   T6 = MIN_EXP;
                  dT6_dVg = dT6_dVb = 0.0;
              }
              else
              {   T6 = std::exp(T5);
                  dT6_dVg = T6 * T12 * (T3 - 2.0 * T4 * Voxacc);
                  dT6_dVb = dT6_dVg * dVoxacc_dVb;
                  dT6_dVg *= dVoxacc_dVg;
              }

              Igbacc = T11 * T2 * T6;
              dIgbacc_dVg = T11 * (T2 * dT6_dVg + T6 * dT2_dVg);
              dIgbacc_dVb = T11 * (T2 * dT6_dVb + T6 * dT2_dVb);


              T0 = tmp * pParam.BSIM4nigbinv;
              T1 = Voxdepinv - pParam.BSIM4eigbinv;
              VxNVt = T1 / T0;
              if (VxNVt > EXP_THRESHOLD)
              {   Vaux = T1;
                  dVaux_dVg = dVoxdepinv_dVg;
                  dVaux_dVd = dVoxdepinv_dVd;
                  dVaux_dVb = dVoxdepinv_dVb;
              }
              else if (VxNVt < -EXP_THRESHOLD)
              {   Vaux = T0 * std::log(1.0 + MIN_EXP);
                  dVaux_dVg = dVaux_dVd = dVaux_dVb = 0.0;
              }
              else
              {   ExpVxNVt = std::exp(VxNVt);
                  Vaux = T0 * std::log(1.0 + ExpVxNVt);
          dVaux_dVg = ExpVxNVt / (1.0 + ExpVxNVt);
                  dVaux_dVd = dVaux_dVg * dVoxdepinv_dVd;
                  dVaux_dVb = dVaux_dVg * dVoxdepinv_dVb;
                  dVaux_dVg *= dVoxdepinv_dVg;
              }

              T2 = (Vgs_eff - Vbseff) * Vaux;
              dT2_dVg = dVgs_eff_dVg * Vaux + (Vgs_eff - Vbseff) * dVaux_dVg;
              dT2_dVd = (Vgs_eff - Vbseff) * dVaux_dVd;
              dT2_dVb = -Vaux + (Vgs_eff - Vbseff) * dVaux_dVb;

              T11 *= 0.75610;
              T12 *= 1.31724;
              T3 = pParam.BSIM4aigbinv * pParam.BSIM4cigbinv
                 - pParam.BSIM4bigbinv;
              T4 = pParam.BSIM4bigbinv * pParam.BSIM4cigbinv;
              T5 = T12 * (pParam.BSIM4aigbinv + T3 * Voxdepinv
                 - T4 * Voxdepinv * Voxdepinv);

              if (T5 > EXP_THRESHOLD)
              {   T6 = MAX_EXP;
                  dT6_dVg = dT6_dVd = dT6_dVb = 0.0;
              }
              else if (T5 < -EXP_THRESHOLD)
              {   T6 = MIN_EXP;
                  dT6_dVg = dT6_dVd = dT6_dVb = 0.0;
              }
              else
              {   T6 = std::exp(T5);
                  dT6_dVg = T6 * T12 * (T3 - 2.0 * T4 * Voxdepinv);
                  dT6_dVd = dT6_dVg * dVoxdepinv_dVd;
                  dT6_dVb = dT6_dVg * dVoxdepinv_dVb;
                  dT6_dVg *= dVoxdepinv_dVg;
              }

              Igbinv = T11 * T2 * T6;
              dIgbinv_dVg = T11 * (T2 * dT6_dVg + T6 * dT2_dVg);
              dIgbinv_dVd = T11 * (T2 * dT6_dVd + T6 * dT2_dVd);
              dIgbinv_dVb = T11 * (T2 * dT6_dVb + T6 * dT2_dVb);

              instance.BSIM4Igb = Igbinv + Igbacc;
              instance.BSIM4gIgbg = dIgbinv_dVg + dIgbacc_dVg;
              instance.BSIM4gIgbd = dIgbinv_dVd;
              instance.BSIM4gIgbb = (dIgbinv_dVb + dIgbacc_dVb) * dVbseff_dVb;
          }
          else
          {  instance.BSIM4Igb = instance.BSIM4gIgbg = instance.BSIM4gIgbd
                = instance.BSIM4gIgbs = instance.BSIM4gIgbb = 0.0;
          } /* End of Gate current */

      if (instance.BSIM4nf != 1.0)
          {   cdrain *= instance.BSIM4nf;
              instance.BSIM4gds *= instance.BSIM4nf;
              instance.BSIM4gm *= instance.BSIM4nf;
              instance.BSIM4gmbs *= instance.BSIM4nf;
          instance.BSIM4IdovVds *= instance.BSIM4nf;

              instance.BSIM4gbbs *= instance.BSIM4nf;
              instance.BSIM4gbgs *= instance.BSIM4nf;
              instance.BSIM4gbds *= instance.BSIM4nf;
              instance.BSIM4csub *= instance.BSIM4nf;

              instance.BSIM4Igidl *= instance.BSIM4nf;
              instance.BSIM4ggidld *= instance.BSIM4nf;
              instance.BSIM4ggidlg *= instance.BSIM4nf;
          instance.BSIM4ggidlb *= instance.BSIM4nf;

              instance.BSIM4Igisl *= instance.BSIM4nf;
              instance.BSIM4ggisls *= instance.BSIM4nf;
              instance.BSIM4ggislg *= instance.BSIM4nf;
              instance.BSIM4ggislb *= instance.BSIM4nf;

              instance.BSIM4Igcs *= instance.BSIM4nf;
              instance.BSIM4gIgcsg *= instance.BSIM4nf;
              instance.BSIM4gIgcsd *= instance.BSIM4nf;
              instance.BSIM4gIgcsb *= instance.BSIM4nf;
              instance.BSIM4Igcd *= instance.BSIM4nf;
              instance.BSIM4gIgcdg *= instance.BSIM4nf;
              instance.BSIM4gIgcdd *= instance.BSIM4nf;
              instance.BSIM4gIgcdb *= instance.BSIM4nf;

              instance.BSIM4Igs *= instance.BSIM4nf;
              instance.BSIM4gIgsg *= instance.BSIM4nf;
              instance.BSIM4gIgss *= instance.BSIM4nf;
              instance.BSIM4Igd *= instance.BSIM4nf;
              instance.BSIM4gIgdg *= instance.BSIM4nf;
              instance.BSIM4gIgdd *= instance.BSIM4nf;

              instance.BSIM4Igb *= instance.BSIM4nf;
              instance.BSIM4gIgbg *= instance.BSIM4nf;
              instance.BSIM4gIgbd *= instance.BSIM4nf;
              instance.BSIM4gIgbb *= instance.BSIM4nf;
      }

          instance.BSIM4ggidls = -(instance.BSIM4ggidld + instance.BSIM4ggidlg
                + instance.BSIM4ggidlb);
          instance.BSIM4ggisld = -(instance.BSIM4ggisls + instance.BSIM4ggislg
                + instance.BSIM4ggislb);
          instance.BSIM4gIgbs = -(instance.BSIM4gIgbg + instance.BSIM4gIgbd
                           + instance.BSIM4gIgbb);
          instance.BSIM4gIgcss = -(instance.BSIM4gIgcsg + instance.BSIM4gIgcsd
                            + instance.BSIM4gIgcsb);
          instance.BSIM4gIgcds = -(instance.BSIM4gIgcdg + instance.BSIM4gIgcdd
                            + instance.BSIM4gIgcdb);
      instance.BSIM4cd = cdrain;


      /* Calculations for noise analysis */

          if (model.BSIM4tnoiMod == 0)
          {   Abulk = Abulk0 * pParam.BSIM4abulkCVfactor;
              Vdsat = Vgsteff / Abulk;
              T0 = Vdsat - Vds - DELTA_4;
              T1 = std::sqrt(T0 * T0 + 4.0 * DELTA_4 * Vdsat);
              if (T0 >= 0.0)
                  Vdseff = Vdsat - 0.5 * (T0 + T1);
              else
              {   T3 = (DELTA_4 + DELTA_4) / (T1 - T0);
                  T4 = 1.0 - T3;
                  T5 = Vdsat * T3 / (T1 - T0);
                  Vdseff = Vdsat * T4;
              }
              if (Vds == 0.0)
                  Vdseff = 0.0;

              T0 = Abulk * Vdseff;
              T1 = 12.0 * (Vgsteff - 0.5 * T0 + 1.0e-20);
              T2 = Vdseff / T1;
              T3 = T0 * T2;
              instance.BSIM4qinv = Coxeff * pParam.BSIM4weffCV * instance.BSIM4nf
                              * pParam.BSIM4leffCV
                              * (Vgsteff - 0.5 * T0 + Abulk * T3);
          }
      else if(model.BSIM4tnoiMod == 2)
      {
              instance.BSIM4noiGd0 = instance.BSIM4nf * beta * Vgsteff / (1.0 + gche * Rds);
      }

      /*
           *  BSIM4 C-V begins
       */

          if ((model.BSIM4xpart < 0) || (!ChargeComputationNeeded))
      {   qgate  = qdrn = qsrc = qbulk = 0.0;
              instance.BSIM4cggb = instance.BSIM4cgsb = instance.BSIM4cgdb = 0.0;
              instance.BSIM4cdgb = instance.BSIM4cdsb = instance.BSIM4cddb = 0.0;
              instance.BSIM4cbgb = instance.BSIM4cbsb = instance.BSIM4cbdb = 0.0;
              instance.BSIM4csgb = instance.BSIM4cssb = instance.BSIM4csdb = 0.0;
              instance.BSIM4cgbb = instance.BSIM4csbb = instance.BSIM4cdbb = instance.BSIM4cbbb = 0.0;
              instance.BSIM4cqdb = instance.BSIM4cqsb = instance.BSIM4cqgb
                              = instance.BSIM4cqbb = 0.0;
              instance.BSIM4gtau = 0.0;
              goto finished;
          }
      else if (model.BSIM4capMod == 0)
      {
              if (Vbseff < 0.0)
          {   VbseffCV = Vbs; /*4.6.2*/
                  dVbseffCV_dVb = 1.0;
              }
          else
          {   VbseffCV = pParam.BSIM4phi - Phis;
                  dVbseffCV_dVb = -dPhis_dVb * dVbseff_dVb; /*4.6.2*/
              }

              Vfb = pParam.BSIM4vfbcv;
              Vth = Vfb + pParam.BSIM4phi + pParam.BSIM4k1ox * sqrtPhis;
              Vgst = Vgs_eff - Vth;
              dVth_dVb = pParam.BSIM4k1ox * dsqrtPhis_dVb *dVbseff_dVb; /*4.6.2*/
              dVgst_dVb = -dVth_dVb;
              dVgst_dVg = dVgs_eff_dVg;

              CoxWL = model.BSIM4coxe * pParam.BSIM4weffCV
                    * pParam.BSIM4leffCV * instance.BSIM4nf;
              Arg1 = Vgs_eff - VbseffCV - Vfb;

              if (Arg1 <= 0.0)
          {   qgate = CoxWL * Arg1;
                  qbulk = -qgate;
                  qdrn = 0.0;

                  instance.BSIM4cggb = CoxWL * dVgs_eff_dVg;
                  instance.BSIM4cgdb = 0.0;
                  instance.BSIM4cgsb = CoxWL * (dVbseffCV_dVb - dVgs_eff_dVg);

                  instance.BSIM4cdgb = 0.0;
                  instance.BSIM4cddb = 0.0;
                  instance.BSIM4cdsb = 0.0;

                  instance.BSIM4cbgb = -CoxWL * dVgs_eff_dVg;
                  instance.BSIM4cbdb = 0.0;
                  instance.BSIM4cbsb = -instance.BSIM4cgsb;
              } /* Arg1 <= 0.0, end of accumulation */
          else if (Vgst <= 0.0)
          {   T1 = 0.5 * pParam.BSIM4k1ox;
              T2 = std::sqrt(T1 * T1 + Arg1);
              qgate = CoxWL * pParam.BSIM4k1ox * (T2 - T1);
                  qbulk = -qgate;
                  qdrn = 0.0;

              T0 = CoxWL * T1 / T2;
              instance.BSIM4cggb = T0 * dVgs_eff_dVg;
              instance.BSIM4cgdb = 0.0;
                  instance.BSIM4cgsb = T0 * (dVbseffCV_dVb - dVgs_eff_dVg);

                  instance.BSIM4cdgb = 0.0;
                  instance.BSIM4cddb = 0.0;
                  instance.BSIM4cdsb = 0.0;

                  instance.BSIM4cbgb = -instance.BSIM4cggb;
                  instance.BSIM4cbdb = 0.0;
                  instance.BSIM4cbsb = -instance.BSIM4cgsb;
              } /* Vgst <= 0.0, end of depletion */
          else
          {   One_Third_CoxWL = CoxWL / 3.0;
                  Two_Third_CoxWL = 2.0 * One_Third_CoxWL;

                  AbulkCV = Abulk0 * pParam.BSIM4abulkCVfactor;
                  dAbulkCV_dVb = pParam.BSIM4abulkCVfactor * dAbulk0_dVb*dVbseff_dVb;

              dVdsat_dVg = 1.0 / AbulkCV;  /*4.6.2*/
              Vdsat = Vgst * dVdsat_dVg;
              dVdsat_dVb = - (Vdsat * dAbulkCV_dVb + dVth_dVb)* dVdsat_dVg;

                  if (model.BSIM4xpart > 0.5)
          {   /* 0/100 Charge partition model */
              if (Vdsat <= Vds)
              {   /* saturation region */
                      T1 = Vdsat / 3.0;
                      qgate = CoxWL * (Vgs_eff - Vfb
                    - pParam.BSIM4phi - T1);
                      T2 = -Two_Third_CoxWL * Vgst;
                      qbulk = -(qgate + T2);
                      qdrn = 0.0;

                      instance.BSIM4cggb = One_Third_CoxWL * (3.0
                      - dVdsat_dVg) * dVgs_eff_dVg;
                      T2 = -One_Third_CoxWL * dVdsat_dVb;
                      instance.BSIM4cgsb = -(instance.BSIM4cggb + T2);
                          instance.BSIM4cgdb = 0.0;

                          instance.BSIM4cdgb = 0.0;
                          instance.BSIM4cddb = 0.0;
                          instance.BSIM4cdsb = 0.0;

                      instance.BSIM4cbgb = -(instance.BSIM4cggb
                      - Two_Third_CoxWL * dVgs_eff_dVg);
                      T3 = -(T2 + Two_Third_CoxWL * dVth_dVb);
                      instance.BSIM4cbsb = -(instance.BSIM4cbgb + T3);
                          instance.BSIM4cbdb = 0.0;
                  }
              else
              {   /* linear region */
              Alphaz = Vgst / Vdsat;
                      T1 = 2.0 * Vdsat - Vds;
                      T2 = Vds / (3.0 * T1);
                      T3 = T2 * Vds;
                      T9 = 0.25 * CoxWL;
                      T4 = T9 * Alphaz;
                      T7 = 2.0 * Vds - T1 - 3.0 * T3;
                      T8 = T3 - T1 - 2.0 * Vds;
                      qgate = CoxWL * (Vgs_eff - Vfb
                    - pParam.BSIM4phi - 0.5 * (Vds - T3));
                      T10 = T4 * T8;
                      qdrn = T4 * T7;
                      qbulk = -(qgate + qdrn + T10);

                          T5 = T3 / T1;
                          instance.BSIM4cggb = CoxWL * (1.0 - T5 * dVdsat_dVg)
                      * dVgs_eff_dVg;
                          T11 = -CoxWL * T5 * dVdsat_dVb;
                          instance.BSIM4cgdb = CoxWL * (T2 - 0.5 + 0.5 * T5);
                          instance.BSIM4cgsb = -(instance.BSIM4cggb + T11
                                          + instance.BSIM4cgdb);
                          T6 = 1.0 / Vdsat;
                          dAlphaz_dVg = T6 * (1.0 - Alphaz * dVdsat_dVg);
                          dAlphaz_dVb = -T6 * (dVth_dVb + Alphaz * dVdsat_dVb);
                          T7 = T9 * T7;
                          T8 = T9 * T8;
                          T9 = 2.0 * T4 * (1.0 - 3.0 * T5);
                          instance.BSIM4cdgb = (T7 * dAlphaz_dVg - T9
                      * dVdsat_dVg) * dVgs_eff_dVg;
                          T12 = T7 * dAlphaz_dVb - T9 * dVdsat_dVb;
                          instance.BSIM4cddb = T4 * (3.0 - 6.0 * T2 - 3.0 * T5);
                          instance.BSIM4cdsb = -(instance.BSIM4cdgb + T12
                                          + instance.BSIM4cddb);

                          T9 = 2.0 * T4 * (1.0 + T5);
                          T10 = (T8 * dAlphaz_dVg - T9 * dVdsat_dVg)
                  * dVgs_eff_dVg;
                          T11 = T8 * dAlphaz_dVb - T9 * dVdsat_dVb;
                          T12 = T4 * (2.0 * T2 + T5 - 1.0);
                          T0 = -(T10 + T11 + T12);

                          instance.BSIM4cbgb = -(instance.BSIM4cggb
                      + instance.BSIM4cdgb + T10);
                          instance.BSIM4cbdb = -(instance.BSIM4cgdb
                      + instance.BSIM4cddb + T12);
                          instance.BSIM4cbsb = -(instance.BSIM4cgsb
                      + instance.BSIM4cdsb + T0);
                      }
                  }
          else if (model.BSIM4xpart < 0.5)
          {   /* 40/60 Charge partition model */
              if (Vds >= Vdsat)
              {   /* saturation region */
                      T1 = Vdsat / 3.0;
                      qgate = CoxWL * (Vgs_eff - Vfb
                    - pParam.BSIM4phi - T1);
                      T2 = -Two_Third_CoxWL * Vgst;
                      qbulk = -(qgate + T2);
                      qdrn = 0.4 * T2;

                      instance.BSIM4cggb = One_Third_CoxWL * (3.0
                      - dVdsat_dVg) * dVgs_eff_dVg;
                      T2 = -One_Third_CoxWL * dVdsat_dVb;
                      instance.BSIM4cgsb = -(instance.BSIM4cggb + T2);
                      instance.BSIM4cgdb = 0.0;

              T3 = 0.4 * Two_Third_CoxWL;
                          instance.BSIM4cdgb = -T3 * dVgs_eff_dVg;
                          instance.BSIM4cddb = 0.0;
              T4 = T3 * dVth_dVb;
                          instance.BSIM4cdsb = -(T4 + instance.BSIM4cdgb);

                      instance.BSIM4cbgb = -(instance.BSIM4cggb
                      - Two_Third_CoxWL * dVgs_eff_dVg);
                      T3 = -(T2 + Two_Third_CoxWL * dVth_dVb);
                      instance.BSIM4cbsb = -(instance.BSIM4cbgb + T3);
                          instance.BSIM4cbdb = 0.0;
                  }
              else
              {   /* linear region  */
              Alphaz = Vgst / Vdsat;
              T1 = 2.0 * Vdsat - Vds;
              T2 = Vds / (3.0 * T1);
              T3 = T2 * Vds;
              T9 = 0.25 * CoxWL;
              T4 = T9 * Alphaz;
              qgate = CoxWL * (Vgs_eff - Vfb - pParam.BSIM4phi
                - 0.5 * (Vds - T3));

              T5 = T3 / T1;
                          instance.BSIM4cggb = CoxWL * (1.0 - T5 * dVdsat_dVg)
                      * dVgs_eff_dVg;
                          tmp = -CoxWL * T5 * dVdsat_dVb;
                          instance.BSIM4cgdb = CoxWL * (T2 - 0.5 + 0.5 * T5);
                          instance.BSIM4cgsb = -(instance.BSIM4cggb
                      + instance.BSIM4cgdb + tmp);

              T6 = 1.0 / Vdsat;
                          dAlphaz_dVg = T6 * (1.0 - Alphaz * dVdsat_dVg);
                          dAlphaz_dVb = -T6 * (dVth_dVb + Alphaz * dVdsat_dVb);

              T6 = 8.0 * Vdsat * Vdsat - 6.0 * Vdsat * Vds
                 + 1.2 * Vds * Vds;
              T8 = T2 / T1;
              T7 = Vds - T1 - T8 * T6;
              qdrn = T4 * T7;
              T7 *= T9;
              tmp = T8 / T1;
              tmp1 = T4 * (2.0 - 4.0 * tmp * T6
                   + T8 * (16.0 * Vdsat - 6.0 * Vds));

                          instance.BSIM4cdgb = (T7 * dAlphaz_dVg - tmp1
                      * dVdsat_dVg) * dVgs_eff_dVg;
                          T10 = T7 * dAlphaz_dVb - tmp1 * dVdsat_dVb;
                          instance.BSIM4cddb = T4 * (2.0 - (1.0 / (3.0 * T1
                      * T1) + 2.0 * tmp) * T6 + T8
                      * (6.0 * Vdsat - 2.4 * Vds));
                          instance.BSIM4cdsb = -(instance.BSIM4cdgb
                      + T10 + instance.BSIM4cddb);

              T7 = 2.0 * (T1 + T3);
              qbulk = -(qgate - T4 * T7);
              T7 *= T9;
              T0 = 4.0 * T4 * (1.0 - T5);
              T12 = (-T7 * dAlphaz_dVg - T0 * dVdsat_dVg) * dVgs_eff_dVg
                      - instance.BSIM4cdgb;  /*4.6.2*/
              T11 = -T7 * dAlphaz_dVb - T10 - T0 * dVdsat_dVb;
              T10 = -4.0 * T4 * (T2 - 0.5 + 0.5 * T5)
                  - instance.BSIM4cddb;
                          tmp = -(T10 + T11 + T12);

                          instance.BSIM4cbgb = -(instance.BSIM4cggb
                      + instance.BSIM4cdgb + T12);
                          instance.BSIM4cbdb = -(instance.BSIM4cgdb
                      + instance.BSIM4cddb + T10);
                          instance.BSIM4cbsb = -(instance.BSIM4cgsb
                      + instance.BSIM4cdsb + tmp);
                      }
                  }
          else
          {   /* 50/50 partitioning */
              if (Vds >= Vdsat)
              {   /* saturation region */
                      T1 = Vdsat / 3.0;
                      qgate = CoxWL * (Vgs_eff - Vfb
                    - pParam.BSIM4phi - T1);
                      T2 = -Two_Third_CoxWL * Vgst;
                      qbulk = -(qgate + T2);
                      qdrn = 0.5 * T2;

                      instance.BSIM4cggb = One_Third_CoxWL * (3.0
                      - dVdsat_dVg) * dVgs_eff_dVg;
                      T2 = -One_Third_CoxWL * dVdsat_dVb;
                      instance.BSIM4cgsb = -(instance.BSIM4cggb + T2);
                      instance.BSIM4cgdb = 0.0;

                          instance.BSIM4cdgb = -One_Third_CoxWL * dVgs_eff_dVg;
                          instance.BSIM4cddb = 0.0;
              T4 = One_Third_CoxWL * dVth_dVb;
                          instance.BSIM4cdsb = -(T4 + instance.BSIM4cdgb);

                      instance.BSIM4cbgb = -(instance.BSIM4cggb
                      - Two_Third_CoxWL * dVgs_eff_dVg);
                      T3 = -(T2 + Two_Third_CoxWL * dVth_dVb);
                      instance.BSIM4cbsb = -(instance.BSIM4cbgb + T3);
                          instance.BSIM4cbdb = 0.0;
                  }
              else
              {   /* linear region */
              Alphaz = Vgst / Vdsat;
              T1 = 2.0 * Vdsat - Vds;
              T2 = Vds / (3.0 * T1);
              T3 = T2 * Vds;
              T9 = 0.25 * CoxWL;
              T4 = T9 * Alphaz;
              qgate = CoxWL * (Vgs_eff - Vfb - pParam.BSIM4phi
                - 0.5 * (Vds - T3));

              T5 = T3 / T1;
                          instance.BSIM4cggb = CoxWL * (1.0 - T5 * dVdsat_dVg)
                      * dVgs_eff_dVg;
                          tmp = -CoxWL * T5 * dVdsat_dVb;
                          instance.BSIM4cgdb = CoxWL * (T2 - 0.5 + 0.5 * T5);
                          instance.BSIM4cgsb = -(instance.BSIM4cggb
                      + instance.BSIM4cgdb + tmp);

              T6 = 1.0 / Vdsat;
                          dAlphaz_dVg = T6 * (1.0 - Alphaz * dVdsat_dVg);
                          dAlphaz_dVb = -T6 * (dVth_dVb + Alphaz * dVdsat_dVb);

              T7 = T1 + T3;
              qdrn = -T4 * T7;
              qbulk = - (qgate + qdrn + qdrn);
              T7 *= T9;
              T0 = T4 * (2.0 * T5 - 2.0);

                          instance.BSIM4cdgb = (T0 * dVdsat_dVg - T7
                      * dAlphaz_dVg) * dVgs_eff_dVg;
              T12 = T0 * dVdsat_dVb - T7 * dAlphaz_dVb;
                          instance.BSIM4cddb = T4 * (1.0 - 2.0 * T2 - T5);
                          instance.BSIM4cdsb = -(instance.BSIM4cdgb + T12
                                          + instance.BSIM4cddb);

                          instance.BSIM4cbgb = -(instance.BSIM4cggb
                      + 2.0 * instance.BSIM4cdgb);
                          instance.BSIM4cbdb = -(instance.BSIM4cgdb
                      + 2.0 * instance.BSIM4cddb);
                          instance.BSIM4cbsb = -(instance.BSIM4cgsb
                      + 2.0 * instance.BSIM4cdsb);
                      } /* end of linear region */
              } /* end of 50/50 partition */
          } /* end of inversion */
          } /* end of capMod=0 */
      else
      {   if (Vbseff < 0.0)
          {   VbseffCV = Vbseff;
                  dVbseffCV_dVb = 1.0;
              }
          else
          {   VbseffCV = pParam.BSIM4phi - Phis;
                  dVbseffCV_dVb = -dPhis_dVb;
              }

              CoxWL = model.BSIM4coxe * pParam.BSIM4weffCV
            * pParam.BSIM4leffCV * instance.BSIM4nf;

          if(model.BSIM4cvchargeMod == 0)
        {
          /* Seperate VgsteffCV with noff and voffcv */
          noff = n * pParam.BSIM4noff;
          dnoff_dVd = pParam.BSIM4noff * dn_dVd;
          dnoff_dVb = pParam.BSIM4noff * dn_dVb;
          T0 = Vtm * noff;
          voffcv = pParam.BSIM4voffcv;
          VgstNVt = (Vgst - voffcv) / T0;

          if (VgstNVt > EXP_THRESHOLD)
            {
              Vgsteff = Vgst - voffcv;
              dVgsteff_dVg = dVgs_eff_dVg;
              dVgsteff_dVd = -dVth_dVd;
              dVgsteff_dVb = -dVth_dVb;
            }
          else if (VgstNVt < -EXP_THRESHOLD)
            {
              Vgsteff = T0 * std::log(1.0 + MIN_EXP);
              dVgsteff_dVg = 0.0;
              dVgsteff_dVd = Vgsteff / noff;
              dVgsteff_dVb = dVgsteff_dVd * dnoff_dVb;
              dVgsteff_dVd *= dnoff_dVd;
            }
          else
            {
              ExpVgst = std::exp(VgstNVt);
              Vgsteff = T0 * std::log(1.0 + ExpVgst);
              dVgsteff_dVg = ExpVgst / (1.0 + ExpVgst);
              dVgsteff_dVd = -dVgsteff_dVg * (dVth_dVd + (Vgst - voffcv)
                   / noff * dnoff_dVd) + Vgsteff / noff * dnoff_dVd;
              dVgsteff_dVb = -dVgsteff_dVg * (dVth_dVb + (Vgst - voffcv)
                   / noff * dnoff_dVb) + Vgsteff / noff * dnoff_dVb;
              dVgsteff_dVg *= dVgs_eff_dVg;
            }
          /* End of VgsteffCV for cvchargeMod = 0 */
        }
          else
        {
          T0 = n * Vtm;
          T1 = pParam.BSIM4mstarcv * Vgst;
          T2 = T1 / T0;
          if (T2 > EXP_THRESHOLD)
            {
              T10 = T1;
              dT10_dVg = pParam.BSIM4mstarcv * dVgs_eff_dVg;
              dT10_dVd = -dVth_dVd * pParam.BSIM4mstarcv;
              dT10_dVb = -dVth_dVb * pParam.BSIM4mstarcv;
            }
          else if (T2 < -EXP_THRESHOLD)
            {
              T10 = Vtm * std::log(1.0 + MIN_EXP);
              dT10_dVg = 0.0;
              dT10_dVd = T10 * dn_dVd;
              dT10_dVb = T10 * dn_dVb;
              T10 *= n;
            }
          else
            {
              ExpVgst = std::exp(T2);
              T3 = Vtm * std::log(1.0 + ExpVgst);
              T10 = n * T3;
              dT10_dVg = pParam.BSIM4mstarcv * ExpVgst / (1.0 + ExpVgst);
              dT10_dVb = T3 * dn_dVb - dT10_dVg * (dVth_dVb + Vgst * dn_dVb / n);
              dT10_dVd = T3 * dn_dVd - dT10_dVg * (dVth_dVd + Vgst * dn_dVd / n);
              dT10_dVg *= dVgs_eff_dVg;
            }

          T1 = pParam.BSIM4voffcbncv - (1.0 - pParam.BSIM4mstarcv) * Vgst;
          T2 = T1 / T0;
          if (T2 < -EXP_THRESHOLD)
            {
              T3 = model.BSIM4coxe * MIN_EXP / pParam.BSIM4cdep0;
              T9 = pParam.BSIM4mstarcv + T3 * n;
              dT9_dVg = 0.0;
              dT9_dVd = dn_dVd * T3;
              dT9_dVb = dn_dVb * T3;
            }
          else if (T2 > EXP_THRESHOLD)
            {
              T3 = model.BSIM4coxe * MAX_EXP / pParam.BSIM4cdep0;
              T9 = pParam.BSIM4mstarcv + T3 * n;
              dT9_dVg = 0.0;
              dT9_dVd = dn_dVd * T3;
              dT9_dVb = dn_dVb * T3;
            }
          else
            {
              ExpVgst = std::exp(T2);
              T3 = model.BSIM4coxe / pParam.BSIM4cdep0;
              T4 = T3 * ExpVgst;
              T5 = T1 * T4 / T0;
              T9 = pParam.BSIM4mstarcv + n * T4;
              dT9_dVg = T3 * (pParam.BSIM4mstarcv - 1.0) * ExpVgst / Vtm;
              dT9_dVb = T4 * dn_dVb - dT9_dVg * dVth_dVb - T5 * dn_dVb;
              dT9_dVd = T4 * dn_dVd - dT9_dVg * dVth_dVd - T5 * dn_dVd;
              dT9_dVg *= dVgs_eff_dVg;
            }

          Vgsteff = T10 / T9;
          T11 = T9 * T9;
          dVgsteff_dVg = (T9 * dT10_dVg - T10 * dT9_dVg) / T11;
          dVgsteff_dVd = (T9 * dT10_dVd - T10 * dT9_dVd) / T11;
          dVgsteff_dVb = (T9 * dT10_dVb - T10 * dT9_dVb) / T11;
          /* End of VgsteffCV for cvchargeMod = 1 */
        }


          if (model.BSIM4capMod == 1)
          {   Vfb = instance.BSIM4vfbzb;
                  V3 = Vfb - Vgs_eff + VbseffCV - DELTA_3;
          if (Vfb <= 0.0)
              T0 = std::sqrt(V3 * V3 - 4.0 * DELTA_3 * Vfb);
          else
              T0 = std::sqrt(V3 * V3 + 4.0 * DELTA_3 * Vfb);

          T1 = 0.5 * (1.0 + V3 / T0);
          Vfbeff = Vfb - 0.5 * (V3 + T0);
          dVfbeff_dVg = T1 * dVgs_eff_dVg;
          dVfbeff_dVb = -T1 * dVbseffCV_dVb;
          Qac0 = CoxWL * (Vfbeff - Vfb);
          dQac0_dVg = CoxWL * dVfbeff_dVg;
          dQac0_dVb = CoxWL * dVfbeff_dVb;

                  T0 = 0.5 * pParam.BSIM4k1ox;
          T3 = Vgs_eff - Vfbeff - VbseffCV - Vgsteff;
                  if (pParam.BSIM4k1ox == 0.0)
                  {   T1 = 0.0;
                      T2 = 0.0;
                  }
          else if (T3 < 0.0)
          {   T1 = T0 + T3 / pParam.BSIM4k1ox;
                      T2 = CoxWL;
          }
          else
          {   T1 = std::sqrt(T0 * T0 + T3);
                      T2 = CoxWL * T0 / T1;
          }

          Qsub0 = CoxWL * pParam.BSIM4k1ox * (T1 - T0);

                  dQsub0_dVg = T2 * (dVgs_eff_dVg - dVfbeff_dVg - dVgsteff_dVg);
                  dQsub0_dVd = -T2 * dVgsteff_dVd;
                  dQsub0_dVb = -T2 * (dVfbeff_dVb + dVbseffCV_dVb
                             + dVgsteff_dVb);

                  AbulkCV = Abulk0 * pParam.BSIM4abulkCVfactor;
                  dAbulkCV_dVb = pParam.BSIM4abulkCVfactor * dAbulk0_dVb;
              VdsatCV = Vgsteff / AbulkCV;

          T0 = VdsatCV - Vds - DELTA_4;
          dT0_dVg = 1.0 / AbulkCV;
          dT0_dVb = -VdsatCV * dAbulkCV_dVb / AbulkCV;
          T1 = std::sqrt(T0 * T0 + 4.0 * DELTA_4 * VdsatCV);
                  dT1_dVg = (T0 + DELTA_4 + DELTA_4) / T1;
                  dT1_dVd = -T0 / T1;
                  dT1_dVb = dT1_dVg * dT0_dVb;
          dT1_dVg *= dT0_dVg;
          if (T0 >= 0.0)
          {   VdseffCV = VdsatCV - 0.5 * (T0 + T1);
                      dVdseffCV_dVg = 0.5 * (dT0_dVg - dT1_dVg);
                      dVdseffCV_dVd = 0.5 * (1.0 - dT1_dVd);
                      dVdseffCV_dVb = 0.5 * (dT0_dVb - dT1_dVb);
          }
                  else
                  {   T3 = (DELTA_4 + DELTA_4) / (T1 - T0);
              T4 = 1.0 - T3;
              T5 = VdsatCV * T3 / (T1 - T0);
              VdseffCV = VdsatCV * T4;
                      dVdseffCV_dVg = dT0_dVg * T4 + T5 * (dT1_dVg - dT0_dVg);
                      dVdseffCV_dVd = T5 * (dT1_dVd + 1.0);
                      dVdseffCV_dVb = dT0_dVb * (T4 - T5) + T5 * dT1_dVb;
                  }

                  if (Vds == 0.0)
                  {  VdseffCV = 0.0;
                     dVdseffCV_dVg = 0.0;
                     dVdseffCV_dVb = 0.0;
                  }

          T0 = AbulkCV * VdseffCV;
          T1 = 12.0 * (Vgsteff - 0.5 * T0 + 1.0e-20);
          T2 = VdseffCV / T1;
          T3 = T0 * T2;

          T4 = (1.0 - 12.0 * T2 * T2 * AbulkCV);
          T5 = (6.0 * T0 * (4.0 * Vgsteff - T0) / (T1 * T1) - 0.5);
          T6 = 12.0 * T2 * T2 * Vgsteff;

          qgate = CoxWL * (Vgsteff - 0.5 * VdseffCV + T3);
          Cgg1 = CoxWL * (T4 + T5 * dVdseffCV_dVg);
          Cgd1 = CoxWL * T5 * dVdseffCV_dVd + Cgg1 * dVgsteff_dVd;
          Cgb1 = CoxWL * (T5 * dVdseffCV_dVb + T6 * dAbulkCV_dVb)
            + Cgg1 * dVgsteff_dVb;
          Cgg1 *= dVgsteff_dVg;

          T7 = 1.0 - AbulkCV;
          qbulk = CoxWL * T7 * (0.5 * VdseffCV - T3);
          T4 = -T7 * (T4 - 1.0);
          T5 = -T7 * T5;
          T6 = -(T7 * T6 + (0.5 * VdseffCV - T3));
          Cbg1 = CoxWL * (T4 + T5 * dVdseffCV_dVg);
          Cbd1 = CoxWL * T5 * dVdseffCV_dVd + Cbg1 * dVgsteff_dVd;
          Cbb1 = CoxWL * (T5 * dVdseffCV_dVb + T6 * dAbulkCV_dVb)
            + Cbg1 * dVgsteff_dVb;
          Cbg1 *= dVgsteff_dVg;

          if (model.BSIM4xpart > 0.5)
            {   /* 0/100 Charge petition model */
              T1 = T1 + T1;
              qsrc = -CoxWL * (0.5 * Vgsteff + 0.25 * T0
                       - T0 * T0 / T1);
              T7 = (4.0 * Vgsteff - T0) / (T1 * T1);
              T4 = -(0.5 + 24.0 * T0 * T0 / (T1 * T1));
              T5 = -(0.25 * AbulkCV - 12.0 * AbulkCV * T0 * T7);
              T6 = -(0.25 * VdseffCV - 12.0 * T0 * VdseffCV * T7);
              Csg = CoxWL * (T4 + T5 * dVdseffCV_dVg);
              Csd = CoxWL * T5 * dVdseffCV_dVd + Csg * dVgsteff_dVd;
              Csb = CoxWL * (T5 * dVdseffCV_dVb + T6 * dAbulkCV_dVb)
            + Csg * dVgsteff_dVb;
              Csg *= dVgsteff_dVg;
            }
          else if (model.BSIM4xpart < 0.5)
            {   /* 40/60 Charge petition model */
              T1 = T1 / 12.0;
              T2 = 0.5 * CoxWL / (T1 * T1);
              T3 = Vgsteff * (2.0 * T0 * T0 / 3.0 + Vgsteff
                      * (Vgsteff - 4.0 * T0 / 3.0))
            - 2.0 * T0 * T0 * T0 / 15.0;
              qsrc = -T2 * T3;
              T7 = 4.0 / 3.0 * Vgsteff * (Vgsteff - T0)
            + 0.4 * T0 * T0;
              T4 = -2.0 * qsrc / T1 - T2 * (Vgsteff * (3.0
                                   * Vgsteff - 8.0 * T0 / 3.0)
                            + 2.0 * T0 * T0 / 3.0);
              T5 = (qsrc / T1 + T2 * T7) * AbulkCV;
              T6 = (qsrc / T1 * VdseffCV + T2 * T7 * VdseffCV);
              Csg = (T4 + T5 * dVdseffCV_dVg);
              Csd = T5 * dVdseffCV_dVd + Csg * dVgsteff_dVd;
              Csb = (T5 * dVdseffCV_dVb + T6 * dAbulkCV_dVb)
            + Csg * dVgsteff_dVb;
              Csg *= dVgsteff_dVg;
            }
          else
            {   /* 50/50 Charge petition model */
              qsrc = -0.5 * (qgate + qbulk);
              Csg = -0.5 * (Cgg1 + Cbg1);
              Csb = -0.5 * (Cgb1 + Cbb1);
              Csd = -0.5 * (Cgd1 + Cbd1);
            }

          qgate += Qac0 + Qsub0;
          qbulk -= (Qac0 + Qsub0);
          qdrn = -(qgate + qbulk + qsrc);

          Cgg = dQac0_dVg + dQsub0_dVg + Cgg1;
          Cgd = dQsub0_dVd + Cgd1;
          Cgb = dQac0_dVb + dQsub0_dVb + Cgb1;

          Cbg = Cbg1 - dQac0_dVg - dQsub0_dVg;
          Cbd = Cbd1 - dQsub0_dVd;
          Cbb = Cbb1 - dQac0_dVb - dQsub0_dVb;

          Cgb *= dVbseff_dVb;
          Cbb *= dVbseff_dVb;
          Csb *= dVbseff_dVb;

          instance.BSIM4cggb = Cgg;
          instance.BSIM4cgsb = -(Cgg + Cgd + Cgb);
          instance.BSIM4cgdb = Cgd;
          instance.BSIM4cdgb = -(Cgg + Cbg + Csg);
          instance.BSIM4cdsb = (Cgg + Cgd + Cgb + Cbg + Cbd + Cbb
                     + Csg + Csd + Csb);
          instance.BSIM4cddb = -(Cgd + Cbd + Csd);
          instance.BSIM4cbgb = Cbg;
          instance.BSIM4cbsb = -(Cbg + Cbd + Cbb);
          instance.BSIM4cbdb = Cbd;
          }

              /* Charge-Thickness capMod (CTM) begins */
          else if (model.BSIM4capMod == 2)
          {   V3 = instance.BSIM4vfbzb - Vgs_eff + VbseffCV - DELTA_3;
          if (instance.BSIM4vfbzb <= 0.0)
              T0 = std::sqrt(V3 * V3 - 4.0 * DELTA_3 * instance.BSIM4vfbzb);
          else
              T0 = std::sqrt(V3 * V3 + 4.0 * DELTA_3 * instance.BSIM4vfbzb);

          T1 = 0.5 * (1.0 + V3 / T0);
          Vfbeff = instance.BSIM4vfbzb - 0.5 * (V3 + T0);
          dVfbeff_dVg = T1 * dVgs_eff_dVg;
          dVfbeff_dVb = -T1 * dVbseffCV_dVb;

                  Cox = instance.BSIM4coxp;
                  Tox = 1.0e8 * instance.BSIM4toxp;
                  T0 = (Vgs_eff - VbseffCV - instance.BSIM4vfbzb) / Tox;
                  dT0_dVg = dVgs_eff_dVg / Tox;
                  dT0_dVb = -dVbseffCV_dVb / Tox;

                  tmp = T0 * pParam.BSIM4acde;
                  if ((-EXP_THRESHOLD < tmp) && (tmp < EXP_THRESHOLD))
                  {   Tcen = pParam.BSIM4ldeb * std::exp(tmp);
                      dTcen_dVg = pParam.BSIM4acde * Tcen;
                      dTcen_dVb = dTcen_dVg * dT0_dVb;
                      dTcen_dVg *= dT0_dVg;
                  }
                  else if (tmp <= -EXP_THRESHOLD)
                  {   Tcen = pParam.BSIM4ldeb * MIN_EXP;
                      dTcen_dVg = dTcen_dVb = 0.0;
                  }
                  else
                  {   Tcen = pParam.BSIM4ldeb * MAX_EXP;
                      dTcen_dVg = dTcen_dVb = 0.0;
                  }

                  LINK = 1.0e-3 * instance.BSIM4toxp;
                  V3 = pParam.BSIM4ldeb - Tcen - LINK;
                  V4 = std::sqrt(V3 * V3 + 4.0 * LINK * pParam.BSIM4ldeb);
                  Tcen = pParam.BSIM4ldeb - 0.5 * (V3 + V4);
                  T1 = 0.5 * (1.0 + V3 / V4);
                  dTcen_dVg *= T1;
                  dTcen_dVb *= T1;

                  Ccen = epssub / Tcen;
                  T2 = Cox / (Cox + Ccen);
                  Coxeff = T2 * Ccen;
                  T3 = -Ccen / Tcen;
                  dCoxeff_dVg = T2 * T2 * T3;
                  dCoxeff_dVb = dCoxeff_dVg * dTcen_dVb;
                  dCoxeff_dVg *= dTcen_dVg;
                  CoxWLcen = CoxWL * Coxeff / model.BSIM4coxe;

                  Qac0 = CoxWLcen * (Vfbeff - instance.BSIM4vfbzb);
                  QovCox = Qac0 / Coxeff;
                  dQac0_dVg = CoxWLcen * dVfbeff_dVg
                            + QovCox * dCoxeff_dVg;
                  dQac0_dVb = CoxWLcen * dVfbeff_dVb
                + QovCox * dCoxeff_dVb;

                  T0 = 0.5 * pParam.BSIM4k1ox;
                  T3 = Vgs_eff - Vfbeff - VbseffCV - Vgsteff;
                  if (pParam.BSIM4k1ox == 0.0)
                  {   T1 = 0.0;
                      T2 = 0.0;
                  }
                  else if (T3 < 0.0)
                  {   T1 = T0 + T3 / pParam.BSIM4k1ox;
                      T2 = CoxWLcen;
                  }
                  else
                  {   T1 = std::sqrt(T0 * T0 + T3);
                      T2 = CoxWLcen * T0 / T1;
                  }

                  Qsub0 = CoxWLcen * pParam.BSIM4k1ox * (T1 - T0);
                  QovCox = Qsub0 / Coxeff;
                  dQsub0_dVg = T2 * (dVgs_eff_dVg - dVfbeff_dVg - dVgsteff_dVg)
                             + QovCox * dCoxeff_dVg;
                  dQsub0_dVd = -T2 * dVgsteff_dVd;
                  dQsub0_dVb = -T2 * (dVfbeff_dVb + dVbseffCV_dVb + dVgsteff_dVb)
                             + QovCox * dCoxeff_dVb;

          /* Gate-bias dependent delta Phis begins */
          if (pParam.BSIM4k1ox <= 0.0)
          {   Denomi = 0.25 * pParam.BSIM4moin * Vtm;
                      T0 = 0.5 * pParam.BSIM4sqrtPhi;
          }
          else
          {   Denomi = pParam.BSIM4moin * Vtm
                 * pParam.BSIM4k1ox * pParam.BSIM4k1ox;
                      T0 = pParam.BSIM4k1ox * pParam.BSIM4sqrtPhi;
          }
                  T1 = 2.0 * T0 + Vgsteff;

          DeltaPhi = Vtm * std::log(1.0 + T1 * Vgsteff / Denomi);
          dDeltaPhi_dVg = 2.0 * Vtm * (T1 -T0) / (Denomi + T1 * Vgsteff);
          /* End of delta Phis */

          /* VgDP = Vgsteff - DeltaPhi */
          T0 = Vgsteff - DeltaPhi - 0.001;
          dT0_dVg = 1.0 - dDeltaPhi_dVg;
          T1 = std::sqrt(T0 * T0 + Vgsteff * 0.004);
          VgDP = 0.5 * (T0 + T1);
          dVgDP_dVg = 0.5 * (dT0_dVg + (T0 * dT0_dVg + 0.002) / T1);

                  Tox += Tox; /* WDLiu: Tcen reevaluated below due to different Vgsteff */
                  T0 = (Vgsteff + instance.BSIM4vtfbphi2) / Tox;
                  tmp = std::exp(model.BSIM4bdos * 0.7 * std::log(T0));
                  T1 = 1.0 + tmp;
                  T2 = model.BSIM4bdos * 0.7 * tmp / (T0 * Tox);
                  Tcen = model.BSIM4ados * 1.9e-9 / T1;
                  dTcen_dVg = -Tcen * T2 / T1;
                  dTcen_dVd = dTcen_dVg * dVgsteff_dVd;
                  dTcen_dVb = dTcen_dVg * dVgsteff_dVb;
                  dTcen_dVg *= dVgsteff_dVg;

          Ccen = epssub / Tcen;
          T0 = Cox / (Cox + Ccen);
          Coxeff = T0 * Ccen;
          T1 = -Ccen / Tcen;
          dCoxeff_dVg = T0 * T0 * T1;
          dCoxeff_dVd = dCoxeff_dVg * dTcen_dVd;
          dCoxeff_dVb = dCoxeff_dVg * dTcen_dVb;
          dCoxeff_dVg *= dTcen_dVg;
          CoxWLcen = CoxWL * Coxeff / model.BSIM4coxe;

                  AbulkCV = Abulk0 * pParam.BSIM4abulkCVfactor;
                  dAbulkCV_dVb = pParam.BSIM4abulkCVfactor * dAbulk0_dVb;
                  VdsatCV = VgDP / AbulkCV;

                  T0 = VdsatCV - Vds - DELTA_4;
                  dT0_dVg = dVgDP_dVg / AbulkCV;
                  dT0_dVb = -VdsatCV * dAbulkCV_dVb / AbulkCV;
                  T1 = std::sqrt(T0 * T0 + 4.0 * DELTA_4 * VdsatCV);
                  dT1_dVg = (T0 + DELTA_4 + DELTA_4) / T1;
                  dT1_dVd = -T0 / T1;
                  dT1_dVb = dT1_dVg * dT0_dVb;
                  dT1_dVg *= dT0_dVg;
                  if (T0 >= 0.0)
                  {   VdseffCV = VdsatCV - 0.5 * (T0 + T1);
                      dVdseffCV_dVg = 0.5 * (dT0_dVg - dT1_dVg);
                      dVdseffCV_dVd = 0.5 * (1.0 - dT1_dVd);
                      dVdseffCV_dVb = 0.5 * (dT0_dVb - dT1_dVb);
                  }
                  else
                  {   T3 = (DELTA_4 + DELTA_4) / (T1 - T0);
                      T4 = 1.0 - T3;
                      T5 = VdsatCV * T3 / (T1 - T0);
                      VdseffCV = VdsatCV * T4;
                      dVdseffCV_dVg = dT0_dVg * T4 + T5 * (dT1_dVg - dT0_dVg);
                      dVdseffCV_dVd = T5 * (dT1_dVd + 1.0);
                      dVdseffCV_dVb = dT0_dVb * (T4 - T5) + T5 * dT1_dVb;
                  }

                  if (Vds == 0.0)
                  {  VdseffCV = 0.0;
                     dVdseffCV_dVg = 0.0;
                     dVdseffCV_dVb = 0.0;
                  }

                  T0 = AbulkCV * VdseffCV;
          T1 = VgDP;
                  T2 = 12.0 * (T1 - 0.5 * T0 + 1.0e-20);
                  T3 = T0 / T2;
                  T4 = 1.0 - 12.0 * T3 * T3;
                  T5 = AbulkCV * (6.0 * T0 * (4.0 * T1 - T0) / (T2 * T2) - 0.5);
          T6 = T5 * VdseffCV / AbulkCV;

                  qgate = CoxWLcen * (T1 - T0 * (0.5 - T3));
          QovCox = qgate / Coxeff;
          Cgg1 = CoxWLcen * (T4 * dVgDP_dVg
               + T5 * dVdseffCV_dVg);
          Cgd1 = CoxWLcen * T5 * dVdseffCV_dVd + Cgg1
               * dVgsteff_dVd + QovCox * dCoxeff_dVd;
          Cgb1 = CoxWLcen * (T5 * dVdseffCV_dVb + T6 * dAbulkCV_dVb)
               + Cgg1 * dVgsteff_dVb + QovCox * dCoxeff_dVb;
          Cgg1 = Cgg1 * dVgsteff_dVg + QovCox * dCoxeff_dVg;


                  T7 = 1.0 - AbulkCV;
                  T8 = T2 * T2;
                  T9 = 12.0 * T7 * T0 * T0 / (T8 * AbulkCV);
                  T10 = T9 * dVgDP_dVg;
                  T11 = -T7 * T5 / AbulkCV;
                  T12 = -(T9 * T1 / AbulkCV + VdseffCV * (0.5 - T0 / T2));

          qbulk = CoxWLcen * T7 * (0.5 * VdseffCV - T0 * VdseffCV / T2);
          QovCox = qbulk / Coxeff;
          Cbg1 = CoxWLcen * (T10 + T11 * dVdseffCV_dVg);
          Cbd1 = CoxWLcen * T11 * dVdseffCV_dVd + Cbg1
               * dVgsteff_dVd + QovCox * dCoxeff_dVd;
          Cbb1 = CoxWLcen * (T11 * dVdseffCV_dVb + T12 * dAbulkCV_dVb)
               + Cbg1 * dVgsteff_dVb + QovCox * dCoxeff_dVb;
          Cbg1 = Cbg1 * dVgsteff_dVg + QovCox * dCoxeff_dVg;

                  if (model.BSIM4xpart > 0.5)
          {   /* 0/100 partition */
              qsrc = -CoxWLcen * (T1 / 2.0 + T0 / 4.0
               - 0.5 * T0 * T0 / T2);
              QovCox = qsrc / Coxeff;
              T2 += T2;
              T3 = T2 * T2;
              T7 = -(0.25 - 12.0 * T0 * (4.0 * T1 - T0) / T3);
              T4 = -(0.5 + 24.0 * T0 * T0 / T3) * dVgDP_dVg;
              T5 = T7 * AbulkCV;
              T6 = T7 * VdseffCV;

              Csg = CoxWLcen * (T4 + T5 * dVdseffCV_dVg);
              Csd = CoxWLcen * T5 * dVdseffCV_dVd + Csg * dVgsteff_dVd
              + QovCox * dCoxeff_dVd;
              Csb = CoxWLcen * (T5 * dVdseffCV_dVb + T6 * dAbulkCV_dVb)
              + Csg * dVgsteff_dVb + QovCox * dCoxeff_dVb;
              Csg = Csg * dVgsteff_dVg + QovCox * dCoxeff_dVg;
                  }
          else if (model.BSIM4xpart < 0.5)
          {   /* 40/60 partition */
              T2 = T2 / 12.0;
              T3 = 0.5 * CoxWLcen / (T2 * T2);
                      T4 = T1 * (2.0 * T0 * T0 / 3.0 + T1 * (T1 - 4.0
                         * T0 / 3.0)) - 2.0 * T0 * T0 * T0 / 15.0;
              qsrc = -T3 * T4;
              QovCox = qsrc / Coxeff;
              T8 = 4.0 / 3.0 * T1 * (T1 - T0) + 0.4 * T0 * T0;
              T5 = -2.0 * qsrc / T2 - T3 * (T1 * (3.0 * T1 - 8.0
             * T0 / 3.0) + 2.0 * T0 * T0 / 3.0);
              T6 = AbulkCV * (qsrc / T2 + T3 * T8);
              T7 = T6 * VdseffCV / AbulkCV;

              Csg = T5 * dVgDP_dVg + T6 * dVdseffCV_dVg;
              Csd = Csg * dVgsteff_dVd + T6 * dVdseffCV_dVd
              + QovCox * dCoxeff_dVd;
              Csb = Csg * dVgsteff_dVb + T6 * dVdseffCV_dVb
                  + T7 * dAbulkCV_dVb + QovCox * dCoxeff_dVb;
              Csg = Csg * dVgsteff_dVg + QovCox * dCoxeff_dVg;
                  }
          else
          {   /* 50/50 partition */
                      qsrc = -0.5 * qgate;
                      Csg = -0.5 * Cgg1;
                      Csd = -0.5 * Cgd1;
                      Csb = -0.5 * Cgb1;
                  }

          qgate += Qac0 + Qsub0 - qbulk;
          qbulk -= (Qac0 + Qsub0);
                  qdrn = -(qgate + qbulk + qsrc);

          Cbg = Cbg1 - dQac0_dVg - dQsub0_dVg;
          Cbd = Cbd1 - dQsub0_dVd;
          Cbb = Cbb1 - dQac0_dVb - dQsub0_dVb;

                  Cgg = Cgg1 - Cbg;
                  Cgd = Cgd1 - Cbd;
                  Cgb = Cgb1 - Cbb;

          Cgb *= dVbseff_dVb;
          Cbb *= dVbseff_dVb;
          Csb *= dVbseff_dVb;

                  instance.BSIM4cggb = Cgg;
              instance.BSIM4cgsb = -(Cgg + Cgd + Cgb);
              instance.BSIM4cgdb = Cgd;
                  instance.BSIM4cdgb = -(Cgg + Cbg + Csg);
              instance.BSIM4cdsb = (Cgg + Cgd + Cgb + Cbg + Cbd + Cbb
                      + Csg + Csd + Csb);
              instance.BSIM4cddb = -(Cgd + Cbd + Csd);
                  instance.BSIM4cbgb = Cbg;
              instance.BSIM4cbsb = -(Cbg + Cbd + Cbb);
              instance.BSIM4cbdb = Cbd;
          }  /* End of CTM */
          }

          instance.BSIM4csgb = - instance.BSIM4cggb - instance.BSIM4cdgb - instance.BSIM4cbgb;
          instance.BSIM4csdb = - instance.BSIM4cgdb - instance.BSIM4cddb - instance.BSIM4cbdb;
          instance.BSIM4cssb = - instance.BSIM4cgsb - instance.BSIM4cdsb - instance.BSIM4cbsb;
          instance.BSIM4cgbb = - instance.BSIM4cgdb - instance.BSIM4cggb - instance.BSIM4cgsb;
          instance.BSIM4cdbb = - instance.BSIM4cddb - instance.BSIM4cdgb - instance.BSIM4cdsb;
          instance.BSIM4cbbb = - instance.BSIM4cbgb - instance.BSIM4cbdb - instance.BSIM4cbsb;
          instance.BSIM4csbb = - instance.BSIM4cgbb - instance.BSIM4cdbb - instance.BSIM4cbbb;
          instance.BSIM4qgate = qgate;
          instance.BSIM4qbulk = qbulk;
          instance.BSIM4qdrn = qdrn;
          instance.BSIM4qsrc = -(qgate + qbulk + qdrn);

          /* NQS begins */
          if ((instance.BSIM4trnqsMod) || (instance.BSIM4acnqsMod))
          {   instance.BSIM4qchqs = qcheq = -(qbulk + qgate);
              instance.BSIM4cqgb = -(instance.BSIM4cggb + instance.BSIM4cbgb);
              instance.BSIM4cqdb = -(instance.BSIM4cgdb + instance.BSIM4cbdb);
              instance.BSIM4cqsb = -(instance.BSIM4cgsb + instance.BSIM4cbsb);
              instance.BSIM4cqbb = -(instance.BSIM4cqgb + instance.BSIM4cqdb
                              + instance.BSIM4cqsb);

              CoxWL = model.BSIM4coxe * pParam.BSIM4weffCV * instance.BSIM4nf
                    * pParam.BSIM4leffCV;
              T1 = instance.BSIM4gcrg / CoxWL; /* 1 / tau */
              instance.BSIM4gtau = T1 * ScalingFactor;

          if (instance.BSIM4acnqsMod)
                  instance.BSIM4taunet = 1.0 / T1;

              instance.BSIM4states0[BSIM4qcheq] = qcheq;
              if (CKTmode & MODEINITTRAN)
                  instance.BSIM4states1[BSIM4qcheq] =
                                   instance.BSIM4states0[BSIM4qcheq];
        //   if (instance.BSIM4trnqsMod)
        //       {   error = NIintegrate(ckt, &geq, &ceq, 0.0, instance.BSIM4qcheq);
        //           if (error)
        //               return(error);
        //   }
          }


finished:

    /* Calculate junction C-V */
    if (ChargeComputationNeeded)
    {   czbd = model.BSIM4DunitAreaTempJctCap * instance.BSIM4Adeff; /* bug fix */
              czbs = model.BSIM4SunitAreaTempJctCap * instance.BSIM4Aseff;
              czbdsw = model.BSIM4DunitLengthSidewallTempJctCap * instance.BSIM4Pdeff;
              czbdswg = model.BSIM4DunitLengthGateSidewallTempJctCap
                      * pParam.BSIM4weffCJ * instance.BSIM4nf;
              czbssw = model.BSIM4SunitLengthSidewallTempJctCap * instance.BSIM4Pseff;
              czbsswg = model.BSIM4SunitLengthGateSidewallTempJctCap
                      * pParam.BSIM4weffCJ * instance.BSIM4nf;

              MJS = model.BSIM4SbulkJctBotGradingCoeff;
              MJSWS = model.BSIM4SbulkJctSideGradingCoeff;
          MJSWGS = model.BSIM4SbulkJctGateSideGradingCoeff;

              MJD = model.BSIM4DbulkJctBotGradingCoeff;
              MJSWD = model.BSIM4DbulkJctSideGradingCoeff;
              MJSWGD = model.BSIM4DbulkJctGateSideGradingCoeff;

              /* Source Bulk Junction */
          if (vbs_jct == 0.0)
          {     instance.BSIM4states0[BSIM4qbs] = 0.0;
                instance.BSIM4capbs = czbs + czbssw + czbsswg;
          }
          else if (vbs_jct < 0.0)
          {   if (czbs > 0.0)
          {   arg = 1.0 - vbs_jct / model.BSIM4PhiBS;
              if (MJS == 0.5)
                          sarg = 1.0 / std::sqrt(arg);
              else
                          sarg = std::exp(-MJS * std::log(arg));
                      instance.BSIM4states0[BSIM4qbs] = model.BSIM4PhiBS * czbs
                           * (1.0 - arg * sarg) / (1.0 - MJS);
              instance.BSIM4capbs = czbs * sarg;
          }
          else
          {   instance.BSIM4states0[BSIM4qbs] = 0.0;
              instance.BSIM4capbs = 0.0;
          }
          if (czbssw > 0.0)
          {   arg = 1.0 - vbs_jct / model.BSIM4PhiBSWS;
              if (MJSWS == 0.5)
                          sarg = 1.0 / std::sqrt(arg);
              else
                          sarg = std::exp(-MJSWS * std::log(arg));
                      instance.BSIM4states0[BSIM4qbs] += model.BSIM4PhiBSWS * czbssw
                       * (1.0 - arg * sarg) / (1.0 - MJSWS);
                      instance.BSIM4capbs += czbssw * sarg;
          }
          if (czbsswg > 0.0)
          {   arg = 1.0 - vbs_jct / model.BSIM4PhiBSWGS;
              if (MJSWGS == 0.5)
                          sarg = 1.0 / std::sqrt(arg);
              else
                          sarg = std::exp(-MJSWGS * std::log(arg));
                      instance.BSIM4states0[BSIM4qbs] += model.BSIM4PhiBSWGS * czbsswg
                       * (1.0 - arg * sarg) / (1.0 - MJSWGS);
                      instance.BSIM4capbs += czbsswg * sarg;
          }

              }
          else
          {   T0 = czbs + czbssw + czbsswg;
                  T1 = vbs_jct * (czbs * MJS / model.BSIM4PhiBS + czbssw * MJSWS
                     / model.BSIM4PhiBSWS + czbsswg * MJSWGS / model.BSIM4PhiBSWGS);
                  instance.BSIM4states0[BSIM4qbs] = vbs_jct * (T0 + 0.5 * T1);
                  instance.BSIM4capbs = T0 + T1;
              }

              /* Drain Bulk Junction */
          if (vbd_jct == 0.0)
          {   instance.BSIM4states0[BSIM4qbd] = 0.0;
                  instance.BSIM4capbd = czbd + czbdsw + czbdswg;
          }
          else if (vbd_jct < 0.0)
          {   if (czbd > 0.0)
          {   arg = 1.0 - vbd_jct / model.BSIM4PhiBD;
              if (MJD == 0.5)
                          sarg = 1.0 / std::sqrt(arg);
              else
                          sarg = std::exp(-MJD * std::log(arg));
                      instance.BSIM4states0[BSIM4qbd] = model.BSIM4PhiBD* czbd
                           * (1.0 - arg * sarg) / (1.0 - MJD);
                      instance.BSIM4capbd = czbd * sarg;
          }
          else
          {     instance.BSIM4states0[BSIM4qbd] = 0.0;
                instance.BSIM4capbd = 0.0;
          }
          if (czbdsw > 0.0)
          {   arg = 1.0 - vbd_jct / model.BSIM4PhiBSWD;
              if (MJSWD == 0.5)
                          sarg = 1.0 / std::sqrt(arg);
              else
                          sarg = std::exp(-MJSWD * std::log(arg));
                      instance.BSIM4states0[BSIM4qbd] += model.BSIM4PhiBSWD * czbdsw
                           * (1.0 - arg * sarg) / (1.0 - MJSWD);
                      instance.BSIM4capbd += czbdsw * sarg;
          }
          if (czbdswg > 0.0)
          {   arg = 1.0 - vbd_jct / model.BSIM4PhiBSWGD;
              if (MJSWGD == 0.5)
                          sarg = 1.0 / std::sqrt(arg);
              else
                          sarg = std::exp(-MJSWGD * std::log(arg));
                      instance.BSIM4states0[BSIM4qbd] += model.BSIM4PhiBSWGD * czbdswg
                       * (1.0 - arg * sarg) / (1.0 - MJSWGD);
                      instance.BSIM4capbd += czbdswg * sarg;
          }
              }
          else
          {   T0 = czbd + czbdsw + czbdswg;
                  T1 = vbd_jct * (czbd * MJD / model.BSIM4PhiBD + czbdsw * MJSWD
                     / model.BSIM4PhiBSWD + czbdswg * MJSWGD / model.BSIM4PhiBSWGD);
                  instance.BSIM4states0[BSIM4qbd] = vbd_jct * (T0 + 0.5 * T1);
                  instance.BSIM4capbd = T0 + T1;
              }
    }


    /*
    *  check convergence
    */

    // if ((instance.BSIM4off == 0) || (!(CKTmode & MODEINITFIX)))
    // {   
    //     if (Check == 1)
    //     {                   
    //         // {   ckt->CKTnoncon++;
    //     }
    //     else
    //     {   
    //         if (instance.BSIM4mode >= 0)
    //         {   Idtot = instance.BSIM4cd + instance.BSIM4csub
    //                 + instance.BSIM4Igidl - instance.BSIM4cbd;
    //         }
    //         else
    //         {   Idtot = instance.BSIM4cd + instance.BSIM4cbd - instance.BSIM4Igidl; /* bugfix */
    //         }
    //         tol0 = ckt->CKTreltol * std::max( std::abs(cdhat), std::abs(Idtot))
    //                 + ckt->CKTabstol;
    //         tol1 = ckt->CKTreltol * std::max( std::abs(cseshat), std::abs(Isestot))
    //                 + ckt->CKTabstol;
    //         tol2 = ckt->CKTreltol * std::max( std::abs(cdedhat), std::abs(Idedtot))
    //             + ckt->CKTabstol;
    //         tol3 = ckt->CKTreltol * std::max( std::abs(cgshat), std::abs(Igstot))
    //             + ckt->CKTabstol;
    //         tol4 = ckt->CKTreltol * std::max( std::abs(cgdhat), std::abs(Igdtot))
    //             + ckt->CKTabstol;
    //         tol5 = ckt->CKTreltol * std::max( std::abs(cgbhat), std::abs(Igbtot))
    //             + ckt->CKTabstol;
    //         if (( std::abs(cdhat - Idtot) >= tol0) || ( std::abs(cseshat - Isestot) >= tol1)
    //                 || ( std::abs(cdedhat - Idedtot) >= tol2))
    //         {   
    //             // ckt->CKTnoncon++;
    //         }
    //         else if (( std::abs(cgshat - Igstot) >= tol3) || ( std::abs(cgdhat - Igdtot) >= tol4)
    //                 || ( std::abs(cgbhat - Igbtot) >= tol5))
    //             {   
    //                 // ckt->CKTnoncon++;
    //             }
    //         else
    //         {   Ibtot = instance.BSIM4cbs + instance.BSIM4cbd
    //             - instance.BSIM4Igidl - instance.BSIM4Igisl - instance.BSIM4csub;
    //                 tol6 = ckt->CKTreltol * std::max( std::abs(cbhat), std::abs(Ibtot))
    //                     + ckt->CKTabstol;
    //             if ( std::abs(cbhat - Ibtot) > tol6)
    //             {   ckt->CKTnoncon++;
    //             }
    //         }
    //     }
    // }
    
    instance.BSIM4states0[BSIM4vds] = vds;
    instance.BSIM4states0[BSIM4vgs] = vgs;
    instance.BSIM4states0[BSIM4vbs] = vbs;
    instance.BSIM4states0[BSIM4vbd] = vbd;
    instance.BSIM4states0[BSIM4vges] = vges;
    instance.BSIM4states0[BSIM4vgms] = vgms;
    instance.BSIM4states0[BSIM4vdbs] = vdbs;
    instance.BSIM4states0[BSIM4vdbd] = vdbd;
    instance.BSIM4states0[BSIM4vsbs] = vsbs;
    instance.BSIM4states0[BSIM4vses] = vses;
    instance.BSIM4states0[BSIM4vdes] = vdes;
    instance.BSIM4states0[BSIM4qdef] = qdef;


    if (!ChargeComputationNeeded)
        goto line850;

    if (instance.BSIM4rgateMod == 3)
    {
        vgdx = vgmd;
        vgsx = vgms;
    }
    else  /* For rgateMod == 0, 1 and 2 */
    {
        vgdx = vgd;
        vgsx = vgs;
    }
    if (model.BSIM4capMod == 0)
    {
        cgdo = pParam.BSIM4cgdo;
        qgdo = pParam.BSIM4cgdo * vgdx;
        cgso = pParam.BSIM4cgso;
        qgso = pParam.BSIM4cgso * vgsx;
    }
    else /* For both capMod == 1 and 2 */
    {   T0 = vgdx + DELTA_1;
        T1 = std::sqrt(T0 * T0 + 4.0 * DELTA_1);
        T2 = 0.5 * (T0 - T1);

        T3 = pParam.BSIM4weffCV * pParam.BSIM4cgdl;
        T4 = std::sqrt(1.0 - 4.0 * T2 / pParam.BSIM4ckappad);
        cgdo = pParam.BSIM4cgdo + T3 - T3 * (1.0 - 1.0 / T4)
        * (0.5 - 0.5 * T0 / T1);
        qgdo = (pParam.BSIM4cgdo + T3) * vgdx - T3 * (T2
        + 0.5 * pParam.BSIM4ckappad * (T4 - 1.0));

        T0 = vgsx + DELTA_1;
        T1 = std::sqrt(T0 * T0 + 4.0 * DELTA_1);
        T2 = 0.5 * (T0 - T1);
        T3 = pParam.BSIM4weffCV * pParam.BSIM4cgsl;
        T4 = std::sqrt(1.0 - 4.0 * T2 / pParam.BSIM4ckappas);
        cgso = pParam.BSIM4cgso + T3 - T3 * (1.0 - 1.0 / T4)
        * (0.5 - 0.5 * T0 / T1);
        qgso = (pParam.BSIM4cgso + T3) * vgsx - T3 * (T2
        + 0.5 * pParam.BSIM4ckappas * (T4 - 1.0));
    }

    if (instance.BSIM4nf != 1.0)
    {   cgdo *= instance.BSIM4nf;
        cgso *= instance.BSIM4nf;
        qgdo *= instance.BSIM4nf;
        qgso *= instance.BSIM4nf;
    }
    instance.BSIM4cgdo = cgdo;
    instance.BSIM4qgdo = qgdo;
    instance.BSIM4cgso = cgso;
    instance.BSIM4qgso = qgso;


line755:
    // ag0 = ckt->CKTag[0];
    ag0 = 0.0;
    if (instance.BSIM4mode > 0)
    {   
        if (instance.BSIM4trnqsMod == 0)
        {   
            qdrn -= qgdo;
            if (instance.BSIM4rgateMod == 3)
            {   
                gcgmgmb = (cgdo + cgso + pParam.BSIM4cgbo) * ag0;
                gcgmdb = -cgdo * ag0;
                gcgmsb = -cgso * ag0;
                gcgmbb = -pParam.BSIM4cgbo * ag0;

                gcdgmb = gcgmdb;
                gcsgmb = gcgmsb;
                gcbgmb = gcgmbb;

                gcggb = instance.BSIM4cggb * ag0;
                gcgdb = instance.BSIM4cgdb * ag0;
                gcgsb = instance.BSIM4cgsb * ag0;
                gcgbb = -(gcggb + gcgdb + gcgsb);

                gcdgb = instance.BSIM4cdgb * ag0;
                gcsgb = -(instance.BSIM4cggb + instance.BSIM4cbgb
                                + instance.BSIM4cdgb) * ag0;
                gcbgb = instance.BSIM4cbgb * ag0;

                qgmb = pParam.BSIM4cgbo * vgmb;
                qgmid = qgdo + qgso + qgmb;
                qbulk -= qgmb;
                qsrc = -(qgate + qgmid + qbulk + qdrn);
            }
            else
            {   gcggb = (instance.BSIM4cggb + cgdo + cgso
                        + pParam.BSIM4cgbo ) * ag0;
                gcgdb = (instance.BSIM4cgdb - cgdo) * ag0;
                gcgsb = (instance.BSIM4cgsb - cgso) * ag0;
                gcgbb = -(gcggb + gcgdb + gcgsb);

                gcdgb = (instance.BSIM4cdgb - cgdo) * ag0;
                gcsgb = -(instance.BSIM4cggb + instance.BSIM4cbgb
                        + instance.BSIM4cdgb + cgso) * ag0;
                gcbgb = (instance.BSIM4cbgb - pParam.BSIM4cgbo) * ag0;

                gcdgmb = gcsgmb = gcbgmb = 0.0;

                qgb = pParam.BSIM4cgbo * vgb;
                qgate += qgdo + qgso + qgb;
                qbulk -= qgb;
                qsrc = -(qgate + qbulk + qdrn);
            }
            gcddb = (instance.BSIM4cddb + instance.BSIM4capbd + cgdo) * ag0;
            gcdsb = instance.BSIM4cdsb * ag0;

            gcsdb = -(instance.BSIM4cgdb + instance.BSIM4cbdb
                + instance.BSIM4cddb) * ag0;
            gcssb = (instance.BSIM4capbs + cgso - (instance.BSIM4cgsb
                + instance.BSIM4cbsb + instance.BSIM4cdsb)) * ag0;

            if (!instance.BSIM4rbodyMod)
            {   gcdbb = -(gcdgb + gcddb + gcdsb + gcdgmb);
                gcsbb = -(gcsgb + gcsdb + gcssb + gcsgmb);
                gcbdb = (instance.BSIM4cbdb - instance.BSIM4capbd) * ag0;
                gcbsb = (instance.BSIM4cbsb - instance.BSIM4capbs) * ag0;
                gcdbdb = 0.0; gcsbsb = 0.0;
            }
            else
            {   gcdbb  = -(instance.BSIM4cddb + instance.BSIM4cdgb
                        + instance.BSIM4cdsb) * ag0;
                gcsbb = -(gcsgb + gcsdb + gcssb + gcsgmb)
                    + instance.BSIM4capbs * ag0;
                gcbdb = instance.BSIM4cbdb * ag0;
                gcbsb = instance.BSIM4cbsb * ag0;

                gcdbdb = -instance.BSIM4capbd * ag0;
                gcsbsb = -instance.BSIM4capbs * ag0;
            }
            gcbbb = -(gcbdb + gcbgb + gcbsb + gcbgmb);

            ggtg = ggtd = ggtb = ggts = 0.0;
            sxpart = 0.6;
            dxpart = 0.4;
            ddxpart_dVd = ddxpart_dVg = ddxpart_dVb = ddxpart_dVs = 0.0;
            dsxpart_dVd = dsxpart_dVg = dsxpart_dVb = dsxpart_dVs = 0.0;
        }
        else
        {   qcheq = instance.BSIM4qchqs;
            CoxWL = model.BSIM4coxe * pParam.BSIM4weffCV * instance.BSIM4nf
                * pParam.BSIM4leffCV;
            T0 = qdef * ScalingFactor / CoxWL;

            ggtg = instance.BSIM4gtg = T0 * instance.BSIM4gcrgg;
            ggtd = instance.BSIM4gtd = T0 * instance.BSIM4gcrgd;
            ggts = instance.BSIM4gts = T0 * instance.BSIM4gcrgs;
            ggtb = instance.BSIM4gtb = T0 * instance.BSIM4gcrgb;
            gqdef = ScalingFactor * ag0;

            gcqgb = instance.BSIM4cqgb * ag0;
            gcqdb = instance.BSIM4cqdb * ag0;
            gcqsb = instance.BSIM4cqsb * ag0;
            gcqbb = instance.BSIM4cqbb * ag0;

            if (fabs(qcheq) <= 1.0e-5 * CoxWL)
            {   if (model.BSIM4xpart < 0.5)
                {   dxpart = 0.4;
                }
                else if (model.BSIM4xpart > 0.5)
                {   dxpart = 0.0;
                }
                else
                {   dxpart = 0.5;
                }
                ddxpart_dVd = ddxpart_dVg = ddxpart_dVb
                            = ddxpart_dVs = 0.0;
            }
            else
            {   dxpart = qdrn / qcheq;
                Cdd = instance.BSIM4cddb;
                Csd = -(instance.BSIM4cgdb + instance.BSIM4cddb
                    + instance.BSIM4cbdb);
                ddxpart_dVd = (Cdd - dxpart * (Cdd + Csd)) / qcheq;
                Cdg = instance.BSIM4cdgb;
                Csg = -(instance.BSIM4cggb + instance.BSIM4cdgb
                    + instance.BSIM4cbgb);
                ddxpart_dVg = (Cdg - dxpart * (Cdg + Csg)) / qcheq;

                Cds = instance.BSIM4cdsb;
                Css = -(instance.BSIM4cgsb + instance.BSIM4cdsb
                    + instance.BSIM4cbsb);
                ddxpart_dVs = (Cds - dxpart * (Cds + Css)) / qcheq;

                ddxpart_dVb = -(ddxpart_dVd + ddxpart_dVg + ddxpart_dVs);
            }
            sxpart = 1.0 - dxpart;
            dsxpart_dVd = -ddxpart_dVd;
            dsxpart_dVg = -ddxpart_dVg;
            dsxpart_dVs = -ddxpart_dVs;
            dsxpart_dVb = -(dsxpart_dVd + dsxpart_dVg + dsxpart_dVs);

            if (instance.BSIM4rgateMod == 3)
            {   gcgmgmb = (cgdo + cgso + pParam.BSIM4cgbo) * ag0;
                gcgmdb = -cgdo * ag0;
                gcgmsb = -cgso * ag0;
                gcgmbb = -pParam.BSIM4cgbo * ag0;

                gcdgmb = gcgmdb;
                gcsgmb = gcgmsb;
                gcbgmb = gcgmbb;

                gcdgb = gcsgb = gcbgb = 0.0;
                gcggb = gcgdb = gcgsb = gcgbb = 0.0;

                qgmb = pParam.BSIM4cgbo * vgmb;
                qgmid = qgdo + qgso + qgmb;
                qgate = 0.0;
                qbulk = -qgmb;
                qdrn = -qgdo;
                qsrc = -(qgmid + qbulk + qdrn);
            }
            else
            {   gcggb = (cgdo + cgso + pParam.BSIM4cgbo ) * ag0;
                gcgdb = -cgdo * ag0;
                gcgsb = -cgso * ag0;
                gcgbb = -pParam.BSIM4cgbo * ag0;

                gcdgb = gcgdb;
                gcsgb = gcgsb;
                gcbgb = gcgbb;
                gcdgmb = gcsgmb = gcbgmb = 0.0;

                qgb = pParam.BSIM4cgbo * vgb;
                qgate = qgdo + qgso + qgb;
                qbulk = -qgb;
                qdrn = -qgdo;
                qsrc = -(qgate + qbulk + qdrn);
            }

            gcddb = (instance.BSIM4capbd + cgdo) * ag0;
            gcdsb = gcsdb = 0.0;
            gcssb = (instance.BSIM4capbs + cgso) * ag0;

            if (!instance.BSIM4rbodyMod)
            {   gcdbb = -(gcdgb + gcddb + gcdgmb);
                gcsbb = -(gcsgb + gcssb + gcsgmb);
                gcbdb = -instance.BSIM4capbd * ag0;
                gcbsb = -instance.BSIM4capbs * ag0;
                gcdbdb = 0.0; gcsbsb = 0.0;
            }
            else
            {   gcdbb = gcsbb = gcbdb = gcbsb = 0.0;
                gcdbdb = -instance.BSIM4capbd * ag0;
                gcsbsb = -instance.BSIM4capbs * ag0;
            }
            gcbbb = -(gcbdb + gcbgb + gcbsb + gcbgmb);
        }
    }
    else
    {   if (instance.BSIM4trnqsMod == 0)
            {   qsrc = qdrn - qgso;
        if (instance.BSIM4rgateMod == 3)
        {   gcgmgmb = (cgdo + cgso + pParam.BSIM4cgbo) * ag0;
            gcgmdb = -cgdo * ag0;
                gcgmsb = -cgso * ag0;
                gcgmbb = -pParam.BSIM4cgbo * ag0;

                    gcdgmb = gcgmdb;
                    gcsgmb = gcgmsb;
                    gcbgmb = gcgmbb;

                    gcggb = instance.BSIM4cggb * ag0;
                    gcgdb = instance.BSIM4cgsb * ag0;
                    gcgsb = instance.BSIM4cgdb * ag0;
                    gcgbb = -(gcggb + gcgdb + gcgsb);

                    gcdgb = -(instance.BSIM4cggb + instance.BSIM4cbgb
                        + instance.BSIM4cdgb) * ag0;
                    gcsgb = instance.BSIM4cdgb * ag0;
                    gcbgb = instance.BSIM4cbgb * ag0;

                    qgmb = pParam.BSIM4cgbo * vgmb;
                    qgmid = qgdo + qgso + qgmb;
                    qbulk -= qgmb;
                    qdrn = -(qgate + qgmid + qbulk + qsrc);
        }
        else
        {   gcggb = (instance.BSIM4cggb + cgdo + cgso
                        + pParam.BSIM4cgbo ) * ag0;
                    gcgdb = (instance.BSIM4cgsb - cgdo) * ag0;
                    gcgsb = (instance.BSIM4cgdb - cgso) * ag0;
                    gcgbb = -(gcggb + gcgdb + gcgsb);

                    gcdgb = -(instance.BSIM4cggb + instance.BSIM4cbgb
                        + instance.BSIM4cdgb + cgdo) * ag0;
                    gcsgb = (instance.BSIM4cdgb - cgso) * ag0;
                    gcbgb = (instance.BSIM4cbgb - pParam.BSIM4cgbo) * ag0;

                    gcdgmb = gcsgmb = gcbgmb = 0.0;

                    qgb = pParam.BSIM4cgbo * vgb;
                    qgate += qgdo + qgso + qgb;
                    qbulk -= qgb;
                    qdrn = -(qgate + qbulk + qsrc);
        }
                gcddb = (instance.BSIM4capbd + cgdo - (instance.BSIM4cgsb
                    + instance.BSIM4cbsb + instance.BSIM4cdsb)) * ag0;
                gcdsb = -(instance.BSIM4cgdb + instance.BSIM4cbdb
                    + instance.BSIM4cddb) * ag0;

                gcsdb = instance.BSIM4cdsb * ag0;
                gcssb = (instance.BSIM4cddb + instance.BSIM4capbs + cgso) * ag0;

        if (!instance.BSIM4rbodyMod)
        {   gcdbb = -(gcdgb + gcddb + gcdsb + gcdgmb);
                    gcsbb = -(gcsgb + gcsdb + gcssb + gcsgmb);
                    gcbdb = (instance.BSIM4cbsb - instance.BSIM4capbd) * ag0;
                    gcbsb = (instance.BSIM4cbdb - instance.BSIM4capbs) * ag0;
                    gcdbdb = 0.0; gcsbsb = 0.0;
                }
                else
                {   gcdbb = -(gcdgb + gcddb + gcdsb + gcdgmb)
            + instance.BSIM4capbd * ag0;
                    gcsbb = -(instance.BSIM4cddb + instance.BSIM4cdgb
                        + instance.BSIM4cdsb) * ag0;
                    gcbdb = instance.BSIM4cbsb * ag0;
                    gcbsb = instance.BSIM4cbdb * ag0;
                    gcdbdb = -instance.BSIM4capbd * ag0;
            gcsbsb = -instance.BSIM4capbs * ag0;
                }
        gcbbb = -(gcbgb + gcbdb + gcbsb + gcbgmb);

                ggtg = ggtd = ggtb = ggts = 0.0;
        sxpart = 0.4;
                dxpart = 0.6;
        ddxpart_dVd = ddxpart_dVg = ddxpart_dVb = ddxpart_dVs = 0.0;
        dsxpart_dVd = dsxpart_dVg = dsxpart_dVb = dsxpart_dVs = 0.0;
            }
            else
            {   qcheq = instance.BSIM4qchqs;
                CoxWL = model.BSIM4coxe * pParam.BSIM4weffCV * instance.BSIM4nf
                    * pParam.BSIM4leffCV;
                T0 = qdef * ScalingFactor / CoxWL;
                ggtg = instance.BSIM4gtg = T0 * instance.BSIM4gcrgg;
                ggts = instance.BSIM4gts = T0 * instance.BSIM4gcrgd;
                ggtd = instance.BSIM4gtd = T0 * instance.BSIM4gcrgs;
                ggtb = instance.BSIM4gtb = T0 * instance.BSIM4gcrgb;
        gqdef = ScalingFactor * ag0;

                gcqgb = instance.BSIM4cqgb * ag0;
                gcqdb = instance.BSIM4cqsb * ag0;
                gcqsb = instance.BSIM4cqdb * ag0;
                gcqbb = instance.BSIM4cqbb * ag0;

                if (fabs(qcheq) <= 1.0e-5 * CoxWL)
                {   if (model.BSIM4xpart < 0.5)
                    {   sxpart = 0.4;
                    }
                    else if (model.BSIM4xpart > 0.5)
                    {   sxpart = 0.0;
                    }
                    else
                    {   sxpart = 0.5;
                    }
                    dsxpart_dVd = dsxpart_dVg = dsxpart_dVb
                                = dsxpart_dVs = 0.0;
                }
                else
                {   sxpart = qdrn / qcheq;
                    Css = instance.BSIM4cddb;
                    Cds = -(instance.BSIM4cgdb + instance.BSIM4cddb
                        + instance.BSIM4cbdb);
                    dsxpart_dVs = (Css - sxpart * (Css + Cds)) / qcheq;
                    Csg = instance.BSIM4cdgb;
                    Cdg = -(instance.BSIM4cggb + instance.BSIM4cdgb
                        + instance.BSIM4cbgb);
                    dsxpart_dVg = (Csg - sxpart * (Csg + Cdg)) / qcheq;

                    Csd = instance.BSIM4cdsb;
                    Cdd = -(instance.BSIM4cgsb + instance.BSIM4cdsb
                        + instance.BSIM4cbsb);
                    dsxpart_dVd = (Csd - sxpart * (Csd + Cdd)) / qcheq;

                    dsxpart_dVb = -(dsxpart_dVd + dsxpart_dVg + dsxpart_dVs);
                }
                dxpart = 1.0 - sxpart;
                ddxpart_dVd = -dsxpart_dVd;
                ddxpart_dVg = -dsxpart_dVg;
                ddxpart_dVs = -dsxpart_dVs;
                ddxpart_dVb = -(ddxpart_dVd + ddxpart_dVg + ddxpart_dVs);

                if (instance.BSIM4rgateMod == 3)
                {   gcgmgmb = (cgdo + cgso + pParam.BSIM4cgbo) * ag0;
                    gcgmdb = -cgdo * ag0;
                    gcgmsb = -cgso * ag0;
                    gcgmbb = -pParam.BSIM4cgbo * ag0;

                    gcdgmb = gcgmdb;
                    gcsgmb = gcgmsb;
                    gcbgmb = gcgmbb;

                    gcdgb = gcsgb = gcbgb = 0.0;
                    gcggb = gcgdb = gcgsb = gcgbb = 0.0;

                    qgmb = pParam.BSIM4cgbo * vgmb;
                    qgmid = qgdo + qgso + qgmb;
                    qgate = 0.0;
                    qbulk = -qgmb;
                    qdrn = -qgdo;
                    qsrc = -qgso;
                }
                else
                {   gcggb = (cgdo + cgso + pParam.BSIM4cgbo ) * ag0;
                    gcgdb = -cgdo * ag0;
                    gcgsb = -cgso * ag0;
                    gcgbb = -pParam.BSIM4cgbo * ag0;

                    gcdgb = gcgdb;
                    gcsgb = gcgsb;
                    gcbgb = gcgbb;
                    gcdgmb = gcsgmb = gcbgmb = 0.0;

                    qgb = pParam.BSIM4cgbo * vgb;
                    qgate = qgdo + qgso + qgb;
                    qbulk = -qgb;
                    qdrn = -qgdo;
                    qsrc = -qgso;
                }

                gcddb = (instance.BSIM4capbd + cgdo) * ag0;
                gcdsb = gcsdb = 0.0;
                gcssb = (instance.BSIM4capbs + cgso) * ag0;
                if (!instance.BSIM4rbodyMod)
                {   gcdbb = -(gcdgb + gcddb + gcdgmb);
                    gcsbb = -(gcsgb + gcssb + gcsgmb);
                    gcbdb = -instance.BSIM4capbd * ag0;
                    gcbsb = -instance.BSIM4capbs * ag0;
                    gcdbdb = 0.0; gcsbsb = 0.0;
                }
                else
                {   gcdbb = gcsbb = gcbdb = gcbsb = 0.0;
                    gcdbdb = -instance.BSIM4capbd * ag0;
                    gcsbsb = -instance.BSIM4capbs * ag0;
                }
                gcbbb = -(gcbdb + gcbgb + gcbsb + gcbgmb);
            }
        }


    if (instance.BSIM4trnqsMod)
    {  instance.BSIM4states0[BSIM4qcdump] = qdef * ScalingFactor;
        if (CKTmode & MODEINITTRAN)
            instance.BSIM4states1[BSIM4qcdump] =
                            instance.BSIM4states0[BSIM4qcdump];
        //   error = NIintegrate(ckt, &geq, &ceq, 0.0, instance.BSIM4qcdump);
        if (error)
            return(error);
    }

    if (ByPass) goto line860;

    instance.BSIM4states0[BSIM4qg] = qgate;
    instance.BSIM4states0[BSIM4qd] = qdrn
                    - instance.BSIM4states0[BSIM4qbd];
    instance.BSIM4states0[BSIM4qs] = qsrc
                    - instance.BSIM4states0[BSIM4qbs];
    if (instance.BSIM4rgateMod == 3)
        instance.BSIM4states0[BSIM4qgmid] = qgmid;

    if (!instance.BSIM4rbodyMod)
    {   instance.BSIM4states0[BSIM4qb] = qbulk
                        + instance.BSIM4states0[BSIM4qbd]
                        + instance.BSIM4states0[BSIM4qbs];
    }
    else
        instance.BSIM4states0[BSIM4qb] = qbulk;


    /* Store small signal parameters */
    if (CKTmode & MODEINITSMSIG)
    {   goto line1000;
    }

    if (!ChargeComputationNeeded)
        goto line850;

    if (CKTmode & MODEINITTRAN)
    {   instance.BSIM4states1[BSIM4qb] =
            instance.BSIM4states0[BSIM4qb];
        instance.BSIM4states1[BSIM4qg] =
            instance.BSIM4states0[BSIM4qg];
        instance.BSIM4states1[BSIM4qd] =
            instance.BSIM4states0[BSIM4qd];
        if (instance.BSIM4rgateMod == 3)
            instance.BSIM4states1[BSIM4qgmid] =
                instance.BSIM4states0[BSIM4qgmid];
        if (instance.BSIM4rbodyMod)
        {  instance.BSIM4states1[BSIM4qbs] =
                            instance.BSIM4states0[BSIM4qbs];
            instance.BSIM4states1[BSIM4qbd] =
                            instance.BSIM4states0[BSIM4qbd];
        }
    }

    // error = NIintegrate(ckt, &geq, &ceq, 0.0, instance.BSIM4qb);
    // if (error)
    //     return(error);
    // error = NIintegrate(ckt, &geq, &ceq, 0.0, instance.BSIM4qg);
    // if (error)
    //     return(error);
    // error = NIintegrate(ckt, &geq, &ceq, 0.0, instance.BSIM4qd);
    // if (error)
    //     return(error);

    // if (instance.BSIM4rgateMod == 3)
    // {   error = NIintegrate(ckt, &geq, &ceq, 0.0, instance.BSIM4qgmid);
    //     if (error) return(error);
    // }

    // if (instance.BSIM4rbodyMod)
    // {   error = NIintegrate(ckt, &geq, &ceq, 0.0, instance.BSIM4qbs);
    //     if (error)
    //         return(error);
    //     error = NIintegrate(ckt, &geq, &ceq, 0.0, instance.BSIM4qbd);
    //     if (error)
    //         return(error);
    // }

    goto line860;


line850:
    /* Zero gcap and ceqcap if (!ChargeComputationNeeded) */
    ceqqg = ceqqb = ceqqd = 0.0;
    ceqqjd = ceqqjs = 0.0;
    cqcheq = cqdef = 0.0;

    gcdgb = gcddb = gcdsb = gcdbb = 0.0;
    gcsgb = gcsdb = gcssb = gcsbb = 0.0;
    gcggb = gcgdb = gcgsb = gcgbb = 0.0;
    gcbdb = gcbgb = gcbsb = gcbbb = 0.0;

    gcgmgmb = gcgmdb = gcgmsb = gcgmbb = 0.0;
    gcdgmb = gcsgmb = gcbgmb = ceqqgmid = 0.0;
        gcdbdb = gcsbsb = 0.0;

    gqdef = gcqgb = gcqdb = gcqsb = gcqbb = 0.0;
        ggtg = ggtd = ggtb = ggts = 0.0;
        sxpart = (1.0 - (dxpart = (instance.BSIM4mode > 0) ? 0.4 : 0.6));
    ddxpart_dVd = ddxpart_dVg = ddxpart_dVb = ddxpart_dVs = 0.0;
    dsxpart_dVd = dsxpart_dVg = dsxpart_dVb = dsxpart_dVs = 0.0;

    if (instance.BSIM4trnqsMod)
    {   CoxWL = model.BSIM4coxe * pParam.BSIM4weffCV * instance.BSIM4nf
            * pParam.BSIM4leffCV;
        T1 = instance.BSIM4gcrg / CoxWL;
        instance.BSIM4gtau = T1 * ScalingFactor;
    }
    else
        instance.BSIM4gtau = 0.0;

    goto line900;


line860:
    /* Calculate equivalent charge current */

    cqgate = instance.BSIM4states0[BSIM4cqg];
    cqbody = instance.BSIM4states0[BSIM4cqb];
    cqdrn = instance.BSIM4states0[BSIM4cqd];

    ceqqg = cqgate - gcggb * vgb + gcgdb * vbd + gcgsb * vbs;
    ceqqd = cqdrn - gcdgb * vgb - gcdgmb * vgmb + (gcddb + gcdbdb)
        * vbd - gcdbdb * vbd_jct + gcdsb * vbs;
    ceqqb = cqbody - gcbgb * vgb - gcbgmb * vgmb
        + gcbdb * vbd + gcbsb * vbs;


    if (instance.BSIM4rgateMod == 3)
        ceqqgmid = instance.BSIM4states0[BSIM4cqgmid]
                + gcgmdb * vbd + gcgmsb * vbs - gcgmgmb * vgmb;
    else
        ceqqgmid = 0.0;

    if (instance.BSIM4rbodyMod)
    {   ceqqjs = instance.BSIM4states0[BSIM4cqbs] + gcsbsb * vbs_jct;
        ceqqjd = instance.BSIM4states0[BSIM4cqbd] + gcdbdb * vbd_jct;
    }

    if (instance.BSIM4trnqsMod)
    {   T0 = ggtg * vgb - ggtd * vbd - ggts * vbs;
        ceqqg += T0;
        T1 = qdef * instance.BSIM4gtau;
        ceqqd -= dxpart * T0 + T1 * (ddxpart_dVg * vgb - ddxpart_dVd
        * vbd - ddxpart_dVs * vbs);
        cqdef = instance.BSIM4states0[BSIM4cqcdump] - gqdef * qdef;
        cqcheq = instance.BSIM4states0[BSIM4cqcheq]
                - (gcqgb * vgb - gcqdb * vbd - gcqsb * vbs) + T0;
    }

    if (CKTmode & MODEINITTRAN)
    {   instance.BSIM4states1[BSIM4cqb] =
                        instance.BSIM4states0[BSIM4cqb];
        instance.BSIM4states1[BSIM4cqg] =
                        instance.BSIM4states0[BSIM4cqg];
        instance.BSIM4states1[BSIM4cqd] =
                        instance.BSIM4states0[BSIM4cqd];

        if (instance.BSIM4rgateMod == 3)
            instance.BSIM4states1[BSIM4cqgmid] =
                            instance.BSIM4states0[BSIM4cqgmid];

        if (instance.BSIM4rbodyMod)
        {   instance.BSIM4states1[BSIM4cqbs] =
                            instance.BSIM4states0[BSIM4cqbs];
            instance.BSIM4states1[BSIM4cqbd] =
                            instance.BSIM4states0[BSIM4cqbd];
        }
    }


    /*
    *  Load current vector
    */

line900:
    if (instance.BSIM4mode >= 0)
    {   
        Gm = instance.BSIM4gm;
        Gmbs = instance.BSIM4gmbs;
        FwdSum = Gm + Gmbs;
        RevSum = 0.0;

        ceqdrn = model.BSIM4type * (cdrain - instance.BSIM4gds * vds
            - Gm * vgs - Gmbs * vbs);
        ceqbd = model.BSIM4type * (instance.BSIM4csub + instance.BSIM4Igidl
            - (instance.BSIM4gbds + instance.BSIM4ggidld) * vds
            - (instance.BSIM4gbgs + instance.BSIM4ggidlg) * vgs
            - (instance.BSIM4gbbs + instance.BSIM4ggidlb) * vbs);
        ceqbs = model.BSIM4type * (instance.BSIM4Igisl + instance.BSIM4ggisls * vds
            - instance.BSIM4ggislg * vgd - instance.BSIM4ggislb * vbd);

        gbbdp = -(instance.BSIM4gbds);
        gbbsp = instance.BSIM4gbds + instance.BSIM4gbgs + instance.BSIM4gbbs;

        gbdpg = instance.BSIM4gbgs;
        gbdpdp = instance.BSIM4gbds;
        gbdpb = instance.BSIM4gbbs;
        gbdpsp = -(gbdpg + gbdpdp + gbdpb);

        gbspg = 0.0;
        gbspdp = 0.0;
        gbspb = 0.0;
        gbspsp = 0.0;

        if (model.BSIM4igcMod)
        {   
            gIstotg = instance.BSIM4gIgsg + instance.BSIM4gIgcsg;
            gIstotd = instance.BSIM4gIgcsd;
            gIstots = instance.BSIM4gIgss + instance.BSIM4gIgcss;
            gIstotb = instance.BSIM4gIgcsb;
            Istoteq = model.BSIM4type * (instance.BSIM4Igs + instance.BSIM4Igcs
              - gIstotg * vgs - instance.BSIM4gIgcsd * vds
              - instance.BSIM4gIgcsb * vbs);

            gIdtotg = instance.BSIM4gIgdg + instance.BSIM4gIgcdg;
            gIdtotd = instance.BSIM4gIgdd + instance.BSIM4gIgcdd;
            gIdtots = instance.BSIM4gIgcds;
            gIdtotb = instance.BSIM4gIgcdb;
            Idtoteq = model.BSIM4type * (instance.BSIM4Igd + instance.BSIM4Igcd
                - instance.BSIM4gIgdg * vgd - instance.BSIM4gIgcdg * vgs
                - instance.BSIM4gIgcdd * vds - instance.BSIM4gIgcdb * vbs);
        }
        else
        {   
            gIstotg = gIstotd = gIstots = gIstotb = Istoteq = 0.0;
            gIdtotg = gIdtotd = gIdtots = gIdtotb = Idtoteq = 0.0;
        }

        if (model.BSIM4igbMod)
        {   gIbtotg = instance.BSIM4gIgbg;
            gIbtotd = instance.BSIM4gIgbd;
            gIbtots = instance.BSIM4gIgbs;
            gIbtotb = instance.BSIM4gIgbb;
            Ibtoteq = model.BSIM4type * (instance.BSIM4Igb
                    - instance.BSIM4gIgbg * vgs - instance.BSIM4gIgbd * vds
                    - instance.BSIM4gIgbb * vbs);
        }
        else
            gIbtotg = gIbtotd = gIbtots = gIbtotb = Ibtoteq = 0.0;

        if ((model.BSIM4igcMod != 0) || (model.BSIM4igbMod != 0))
        {   gIgtotg = gIstotg + gIdtotg + gIbtotg;
            gIgtotd = gIstotd + gIdtotd + gIbtotd ;
            gIgtots = gIstots + gIdtots + gIbtots;
            gIgtotb = gIstotb + gIdtotb + gIbtotb;
            Igtoteq = Istoteq + Idtoteq + Ibtoteq;
        }
        else
            gIgtotg = gIgtotd = gIgtots = gIgtotb = Igtoteq = 0.0;


        if (instance.BSIM4rgateMod == 2)
            T0 = vges - vgs;
        else if (instance.BSIM4rgateMod == 3)
            T0 = vgms - vgs;
        if (instance.BSIM4rgateMod > 1)
        {   gcrgd = instance.BSIM4gcrgd * T0;
            gcrgg = instance.BSIM4gcrgg * T0;
            gcrgs = instance.BSIM4gcrgs * T0;
            gcrgb = instance.BSIM4gcrgb * T0;
            ceqgcrg = -(gcrgd * vds + gcrgg * vgs
                    + gcrgb * vbs);
            gcrgg -= instance.BSIM4gcrg;
            gcrg = instance.BSIM4gcrg;
        }
        else
        ceqgcrg = gcrg = gcrgd = gcrgg = gcrgs = gcrgb = 0.0;
    }
    else
    {   
        Gm = -instance.BSIM4gm;
        Gmbs = -instance.BSIM4gmbs;
        FwdSum = 0.0;
        RevSum = -(Gm + Gmbs);

        ceqdrn = -model.BSIM4type * (cdrain + instance.BSIM4gds * vds
                + Gm * vgd + Gmbs * vbd);

        ceqbs = model.BSIM4type * (instance.BSIM4csub + instance.BSIM4Igisl
            + (instance.BSIM4gbds + instance.BSIM4ggisls) * vds
            - (instance.BSIM4gbgs + instance.BSIM4ggislg) * vgd
            - (instance.BSIM4gbbs + instance.BSIM4ggislb) * vbd);
        ceqbd = model.BSIM4type * (instance.BSIM4Igidl - instance.BSIM4ggidld * vds
            - instance.BSIM4ggidlg * vgs - instance.BSIM4ggidlb * vbs);

        gbbsp = -(instance.BSIM4gbds);
        gbbdp = instance.BSIM4gbds + instance.BSIM4gbgs + instance.BSIM4gbbs;

        gbdpg = 0.0;
        gbdpsp = 0.0;
        gbdpb = 0.0;
        gbdpdp = 0.0;

        gbspg = instance.BSIM4gbgs;
        gbspsp = instance.BSIM4gbds;
        gbspb = instance.BSIM4gbbs;
        gbspdp = -(gbspg + gbspsp + gbspb);

        if (model.BSIM4igcMod)
        {   
            gIstotg = instance.BSIM4gIgsg + instance.BSIM4gIgcdg;
            gIstotd = instance.BSIM4gIgcds;
            gIstots = instance.BSIM4gIgss + instance.BSIM4gIgcdd;
            gIstotb = instance.BSIM4gIgcdb;
            Istoteq = model.BSIM4type * (instance.BSIM4Igs + instance.BSIM4Igcd
                    - instance.BSIM4gIgsg * vgs - instance.BSIM4gIgcdg * vgd
            + instance.BSIM4gIgcdd * vds - instance.BSIM4gIgcdb * vbd);

            gIdtotg = instance.BSIM4gIgdg + instance.BSIM4gIgcsg;
            gIdtotd = instance.BSIM4gIgdd + instance.BSIM4gIgcss;
            gIdtots = instance.BSIM4gIgcsd;
            gIdtotb = instance.BSIM4gIgcsb;
            Idtoteq = model.BSIM4type * (instance.BSIM4Igd + instance.BSIM4Igcs
                    - (instance.BSIM4gIgdg + instance.BSIM4gIgcsg) * vgd
                    + instance.BSIM4gIgcsd * vds - instance.BSIM4gIgcsb * vbd);
        }
        else
        {   gIstotg = gIstotd = gIstots = gIstotb = Istoteq = 0.0;
            gIdtotg = gIdtotd = gIdtots = gIdtotb = Idtoteq = 0.0;
        }

        if (model.BSIM4igbMod)
        {   gIbtotg = instance.BSIM4gIgbg;
            gIbtotd = instance.BSIM4gIgbs;
            gIbtots = instance.BSIM4gIgbd;
            gIbtotb = instance.BSIM4gIgbb;
            Ibtoteq = model.BSIM4type * (instance.BSIM4Igb
                    - instance.BSIM4gIgbg * vgd + instance.BSIM4gIgbd * vds
                    - instance.BSIM4gIgbb * vbd);
        }
        else
            gIbtotg = gIbtotd = gIbtots = gIbtotb = Ibtoteq = 0.0;

        if ((model.BSIM4igcMod != 0) || (model.BSIM4igbMod != 0))
        {   gIgtotg = gIstotg + gIdtotg + gIbtotg;
            gIgtotd = gIstotd + gIdtotd + gIbtotd ;
            gIgtots = gIstots + gIdtots + gIbtots;
            gIgtotb = gIstotb + gIdtotb + gIbtotb;
            Igtoteq = Istoteq + Idtoteq + Ibtoteq;
        }
        else
            gIgtotg = gIgtotd = gIgtots = gIgtotb = Igtoteq = 0.0;


        if (instance.BSIM4rgateMod == 2)
            T0 = vges - vgs;
        else if (instance.BSIM4rgateMod == 3)
            T0 = vgms - vgs;
        if (instance.BSIM4rgateMod > 1)
        {   gcrgd = instance.BSIM4gcrgs * T0;
            gcrgg = instance.BSIM4gcrgg * T0;
            gcrgs = instance.BSIM4gcrgd * T0;
            gcrgb = instance.BSIM4gcrgb * T0;
            ceqgcrg = -(gcrgg * vgd - gcrgs * vds
                    + gcrgb * vbd);
            gcrgg -= instance.BSIM4gcrg;
            gcrg = instance.BSIM4gcrg;
        }
        else
            ceqgcrg = gcrg = gcrgd = gcrgg = gcrgs = gcrgb = 0.0;
    }

        if (model.BSIM4rdsMod == 1)
        {   ceqgstot = model.BSIM4type * (instance.BSIM4gstotd * vds
                    + instance.BSIM4gstotg * vgs + instance.BSIM4gstotb * vbs);
            /* WDLiu: ceqgstot flowing away from sNodePrime */
            gstot = instance.BSIM4gstot;
            gstotd = instance.BSIM4gstotd;
            gstotg = instance.BSIM4gstotg;
            gstots = instance.BSIM4gstots - gstot;
            gstotb = instance.BSIM4gstotb;

            ceqgdtot = -model.BSIM4type * (instance.BSIM4gdtotd * vds
                    + instance.BSIM4gdtotg * vgs + instance.BSIM4gdtotb * vbs);
            /* WDLiu: ceqgdtot defined as flowing into dNodePrime */
            gdtot = instance.BSIM4gdtot;
            gdtotd = instance.BSIM4gdtotd - gdtot;
            gdtotg = instance.BSIM4gdtotg;
            gdtots = instance.BSIM4gdtots;
            gdtotb = instance.BSIM4gdtotb;
        }
        else
        {   gstot = gstotd = gstotg = gstots = gstotb = ceqgstot = 0.0;
            gdtot = gdtotd = gdtotg = gdtots = gdtotb = ceqgdtot = 0.0;
        }

       if (model.BSIM4type > 0)
           {   ceqjs = (instance.BSIM4cbs - instance.BSIM4gbs * vbs_jct);
               ceqjd = (instance.BSIM4cbd - instance.BSIM4gbd * vbd_jct);
           }
       else
        {   
            ceqjs = -(instance.BSIM4cbs - instance.BSIM4gbs * vbs_jct);
            ceqjd = -(instance.BSIM4cbd - instance.BSIM4gbd * vbd_jct);
            ceqqg = -ceqqg;
            ceqqd = -ceqqd;
            ceqqb = -ceqqb;
            ceqgcrg = -ceqgcrg;

            if (instance.BSIM4trnqsMod)
            {   cqdef = -cqdef;
                cqcheq = -cqcheq;
            }

            if (instance.BSIM4rbodyMod)
            {   ceqqjs = -ceqqjs;
                ceqqjd = -ceqqjd;
            }

           if (instance.BSIM4rgateMod == 3)
           ceqqgmid = -ceqqgmid;
       }


    /*
    *  Loading RHS
    */
    if(!(dNodePrime < 0)){
        RHS(dNodePrime) += (ceqjd - ceqbd + ceqgdtot - ceqdrn - ceqqd + Idtoteq);
    }
    if(!(gNodePrime < 0)){
        RHS(gNodePrime) -= (ceqqg - ceqgcrg + Igtoteq);
    }
   

    // if (instance.BSIM4rgateMod == 2)
    //     RHS(instance.BSIM4gNodeExt - 1) -= ceqgcrg;
    // else if (instance.BSIM4rgateMod == 3)
    //     RHS(instance.BSIM4gNodeMid - 1) -= (ceqqgmid + ceqgcrg);

    if (!instance.BSIM4rbodyMod)
    {   
        if(!(bNodePrime < 0)){
            RHS(bNodePrime) += (ceqbd + ceqbs - ceqjd - ceqjs - ceqqb + Ibtoteq);
        }
        if(!(sNodePrime < 0)){
            RHS(sNodePrime) += (ceqdrn - ceqbs + ceqjs + ceqqg + ceqqb + ceqqd + ceqqgmid - ceqgstot + Istoteq);
        }       
    }
    // else
    // {   
    //     RHS(instance.BSIM4dbNode - 1) -= (ceqjd + ceqqjd);
    //     RHS(instance.BSIM4bNodePrime - 1) += (ceqbd + ceqbs - ceqqb + Ibtoteq);
    //     RHS(instance.BSIM4sbNode - 1) -= (ceqjs + ceqqjs);
    //     RHS(instance.BSIM4sNodePrime - 1) += (ceqdrn - ceqbs + ceqjs + ceqqd
    //         + ceqqg + ceqqb + ceqqjd + ceqqjs + ceqqgmid - ceqgstot + Istoteq);
    // }

    // if (model.BSIM4rdsMod)
    // {   
    //     RHS(instance.BSIM4dNode - 1) -= ceqgdtot;
    //     RHS(instance.BSIM4sNode - 1) += ceqgstot;
    // }

    // if (instance.BSIM4trnqsMod)
    //     RHS(instance.BSIM4qNode - 1) += (cqcheq - cqdef);


    /*
    *  Loading matrix
    */

    if (!instance.BSIM4rbodyMod)
        {   gjbd = instance.BSIM4gbd;
            gjbs = instance.BSIM4gbs;
        }
    else
        gjbd = gjbs = 0.0;

    if (!model.BSIM4rdsMod)
    {   gdpr = instance.BSIM4drainConductance;
        gspr = instance.BSIM4sourceConductance;
    }
    else
        gdpr = gspr = 0.0;

    geltd = instance.BSIM4grgeltd;

    T1 = qdef * instance.BSIM4gtau;

    // // if (instance.BSIM4rgateMod == 1)
    // // {   
    //     //     (*(instance.BSIM4GEgePtr) += geltd);
    //     //     (*(instance.BSIM4GPgePtr) -= geltd);
    //     //     (*(instance.BSIM4GEgpPtr) -= geltd);
    //     //    (*(instance.BSIM4GPgpPtr) += gcggb + geltd - ggtg + gIgtotg);
    //     //    (*(instance.BSIM4GPdpPtr) += gcgdb - ggtd + gIgtotd);
    //     //    (*(instance.BSIM4GPspPtr) += gcgsb - ggts + gIgtots);
    //     //    (*(instance.BSIM4GPbpPtr) += gcgbb - ggtb + gIgtotb);
            
    //     // } /* WDLiu: gcrg already subtracted from all gcrgg below */
    // //     else if (instance.BSIM4rgateMod == 2)
    // //    {   
    //     //     (*(instance.BSIM4GEgePtr) += gcrg);
    //     //     (*(instance.BSIM4GEgpPtr) += gcrgg);
    //     //     (*(instance.BSIM4GEdpPtr) += gcrgd);
    //     //     (*(instance.BSIM4GEspPtr) += gcrgs);
    //     //     (*(instance.BSIM4GEbpPtr) += gcrgb);

    //     //     (*(instance.BSIM4GPgePtr) -= gcrg);
    //     //    (*(instance.BSIM4GPgpPtr) += gcggb  - gcrgg - ggtg + gIgtotg);
    //     //    (*(instance.BSIM4GPdpPtr) += gcgdb - gcrgd - ggtd + gIgtotd);
    //     //    (*(instance.BSIM4GPspPtr) += gcgsb - gcrgs - ggts + gIgtots);
    //     //    (*(instance.BSIM4GPbpPtr) += gcgbb - gcrgb - ggtb + gIgtotb);
    // //    }
    // //    else if (instance.BSIM4rgateMod == 3)
    // //    {   
    //         // (*(instance.BSIM4GEgePtr) += geltd);
    //         //    (*(instance.BSIM4GEgmPtr) -= geltd);
    //         //    (*(instance.BSIM4GMgePtr) -= geltd);
    //         //    (*(instance.BSIM4GMgmPtr) += geltd + gcrg + gcgmgmb);

    //         //    (*(instance.BSIM4GMdpPtr) += gcrgd + gcgmdb);
    //         //    (*(instance.BSIM4GMgpPtr) += gcrgg);
    //         //    (*(instance.BSIM4GMspPtr) += gcrgs + gcgmsb);
    //         //    (*(instance.BSIM4GMbpPtr) += gcrgb + gcgmbb);

    //         //    (*(instance.BSIM4DPgmPtr) += gcdgmb);
    //         //    (*(instance.BSIM4GPgmPtr) -= gcrg);
    //         //    (*(instance.BSIM4SPgmPtr) += gcsgmb);
    //         //    (*(instance.BSIM4BPgmPtr) += gcbgmb);

    //         //    (*(instance.BSIM4GPgpPtr) += gcggb - gcrgg - ggtg + gIgtotg);
    //         //    (*(instance.BSIM4GPdpPtr) += gcgdb - gcrgd - ggtd + gIgtotd);
    //         //    (*(instance.BSIM4GPspPtr) += gcgsb - gcrgs - ggts + gIgtots);
    //         //    (*(instance.BSIM4GPbpPtr) += gcgbb - gcrgb - ggtb + gIgtotb);
    //    }
    // //    else
    // //    {   
                if(!(gNodePrime < 0)){
                    LHS(gNodePrime, gNodePrime) += gcggb - ggtg + gIgtotg; // BSIM4GPgpPtr
                }
                if(!(dNodePrime < 0) && !(gNodePrime < 0)){
                    LHS(gNodePrime, dNodePrime) += gcgdb - ggtd + gIgtotd; // BSIM4GPdpPtr
                }
                if(!(sNodePrime < 0) && !(gNodePrime < 0)){
                    LHS(gNodePrime, sNodePrime) += gcgsb - ggts + gIgtots; // BSIM4GPspPtr
                }
                if(!(bNodePrime < 0) && !(gNodePrime < 0)){
                    LHS(gNodePrime, bNodePrime) += gcgbb - ggtb + gIgtotb; // BSIM4GPbpPtr
                }
            
            
            
            
            
    // //    }

    // //    if (model.BSIM4rdsMod)
    // //    {   (*(instance.BSIM4DgpPtr) += gdtotg);
    // //        (*(instance.BSIM4DspPtr) += gdtots);
    // //         (*(instance.BSIM4DbpPtr) += gdtotb);
    // //         (*(instance.BSIM4SdpPtr) += gstotd);
    // //         (*(instance.BSIM4SgpPtr) += gstotg);
    // //         (*(instance.BSIM4SbpPtr) += gstotb);
    // //    }

    if(!(dNodePrime < 0)){
        // BSIM4DPdpPtr
        LHS(dNodePrime, dNodePrime) += gdpr + instance.BSIM4gds + instance.BSIM4gbd + T1 * ddxpart_dVd
            - gdtotd + RevSum + gcddb + gbdpdp + dxpart * ggtd - gIdtotd; 
    
        if(!(dNode < 0)){
            // BSIM4DPdPtr
            LHS(dNodePrime, dNode) -= gdtotg + gIdtotg; 
            // BSIM4DdpPtr
            LHS(dNode, dNodePrime) -= gdpr - gdtotd;
            // BSIM4DdPtr
            LHS(dNode, dNode) += gdpr + gdtot; 
        }
        if(!(gNodePrime < 0)){
            // BSIM4DPgpPtr
            LHS(dNodePrime, gNodePrime) += Gm + gcdgb - gdtotg + gbdpg - gIdtotg
                + dxpart * ggtg + T1 * ddxpart_dVg; 
        }
        if(!(sNodePrime < 0)){
            // BSIM4DPspPtr
            LHS(dNodePrime, sNodePrime) -= instance.BSIM4gds + gdtots - dxpart * ggts + gIdtots
                - T1 * ddxpart_dVs + FwdSum - gcdsb - gbdpsp; 
            // BSIM4SPdpPtr
            LHS(sNodePrime, dNodePrime) -= instance.BSIM4gds + gstotd + RevSum - gcsdb - gbspdp
                - T1 * dsxpart_dVd - sxpart * ggtd + gIstotd; 
        }
        if(!(bNodePrime < 0)){
            // BSIM4DPbpPtr
            LHS(dNodePrime, bNodePrime) -= gjbd + gdtotb - Gmbs - gcdbb - gbdpb + gIdtotb
                - T1 * ddxpart_dVb - dxpart * ggtb; 
 
        }
    }


    if(!(gNodePrime < 0)){
        if(!(sNodePrime < 0)){
            // BSIM4SPgpPtr
            LHS(sNodePrime, gNodePrime)  += gcsgb - Gm - gstotg + gbspg + sxpart * ggtg
                + T1 * dsxpart_dVg - gIstotg; 
        }
    }

    if(!(sNodePrime < 0)){
        // BSIM4SPspPtr
        LHS(sNodePrime, sNodePrime) += gspr + instance.BSIM4gds + instance.BSIM4gbs + T1 * dsxpart_dVs
            - gstots + FwdSum + gcssb + gbspsp + sxpart * ggts - gIstots; 
        if(!(sNode < 0)){
            // BSIM4SPsPtr
            LHS(sNodePrime, sNode) -= gspr + gstot;
            // BSIM4SspPtr
            LHS(sNode, sNodePrime)  -= gspr - gstots;  
        }
        if(!(bNodePrime < 0)){
            //BSIM4SPbpPtr
            LHS(sNodePrime, bNodePrime) -= gjbs + gstotb + Gmbs - gcsbb - gbspb - sxpart * ggtb
                - T1 * dsxpart_dVb + gIstotb;
        }
    }
        

    if(!(sNode < 0)){
        // BSIM4SsPtr
        LHS(sNode, sNode) += gspr + gstot;
    }

    // BSIM4BPdpPtr
    if(!(bNodePrime < 0) && !(dNodePrime < 0)){
        LHS(bNodePrime, dNodePrime) += gcbdb - gjbd + gbbdp - gIbtotd;
    }
   
    // BSIM4BPgpPtr
    if(!(bNodePrime < 0) && !(gNodePrime < 0)){
        LHS(bNodePrime, gNodePrime) += gcbgb - instance.BSIM4gbgs - gIbtotg; 
    }
   
    // BSIM4BPspPtr
    if(!(bNodePrime < 0) && !(sNodePrime < 0)){
        LHS(bNodePrime, sNodePrime) += gcbsb - gjbs + gbbsp - gIbtots; 
    }

    // BSIM4BPbpPtr
    if(!(bNodePrime < 0)){
        LHS(bNodePrime, bNodePrime) += gjbd + gjbs + gcbbb - instance.BSIM4gbbs
            - gIbtotb; 
    }

    ggidld = instance.BSIM4ggidld;
    ggidlg = instance.BSIM4ggidlg;
    ggidlb = instance.BSIM4ggidlb;
    ggislg = instance.BSIM4ggislg;
    ggisls = instance.BSIM4ggisls;
    ggislb = instance.BSIM4ggislb;
    

    /* stamp gidl */
    // BSIM4DPdpPtr
    if(!(dNodePrime < 0)){
        LHS(dNodePrime, dNodePrime) += ggidld;
    }
    // BSIM4DPgpPtr
    if(!(dNodePrime < 0) && !(gNodePrime < 0)){
        LHS(dNodePrime, gNodePrime) += ggidlg;
    }
    // BSIM4DPspPtr
    if(!(dNodePrime < 0) && !(sNodePrime < 0)){
        LHS(dNodePrime, sNodePrime) -= ggidlg + ggidld + ggidlb;
    }
    // BSIM4DPbpPtr
    if(!(dNodePrime < 0) && !(bNodePrime < 0)){
        LHS(dNodePrime, bNodePrime) += ggidlb;
    }
    // BSIM4BPdpPtr
    if(!(bNodePrime < 0) && !(dNodePrime < 0)){
        LHS(bNodePrime, dNodePrime) -= ggidld;
    }
    // BSIM4BPgpPtr
    if(!(bNodePrime < 0) && !(gNodePrime < 0)){
        LHS(bNodePrime, gNodePrime) -= ggidlg;
    }
    // BSIM4BPspPtr
    if(!(bNodePrime < 0) && !(sNodePrime < 0)){
        LHS(bNodePrime, sNodePrime) += ggidlg + ggidld + ggidlb;
    }
    // BSIM4BPbpPtr
    if(!(bNodePrime < 0)){
        LHS(bNodePrime, bNodePrime) -= ggidlb;
    }
    
    // /* stamp gisl */
    // BSIM4SPdpPtr
    if(!(sNodePrime < 0) && !(dNodePrime < 0)){
        LHS(sNodePrime, dNodePrime) -= ggisls + ggislg + ggislb;
    }
    // BSIM4SPgpPtr
    if(!(sNodePrime < 0) && !(gNodePrime < 0)){
        LHS(sNodePrime, gNodePrime) += ggislg;
    }
    // BSIM4SPspPtr
    if(!(sNodePrime < 0)){
        LHS(sNodePrime, sNodePrime) += ggisls;
    }
    // BSIM4SPbpPtr
    if(!(sNodePrime < 0) && !(bNodePrime < 0)){
        LHS(sNodePrime, bNodePrime) += ggislb;
    }
    // BSIM4BPdpPtr
    if(!(bNodePrime < 0) && !(dNodePrime < 0)){
        LHS(bNodePrime, dNodePrime) += ggislg + ggisls + ggislb;
    }
    // BSIM4BPgpPtr
    if(!(bNodePrime < 0) && !(gNodePrime < 0)){
        LHS(bNodePrime, gNodePrime) -= ggislg;
    }
    // BSIM4BPspPtr
    if(!(bNodePrime < 0) && !(sNodePrime < 0)){
        LHS(bNodePrime, sNodePrime) -= ggisls;
    }
    // BSIM4BPbpPtr
    if(!(bNodePrime < 0)){
        LHS(bNodePrime, bNodePrime) -= ggislb;
    }


   


    // if (instance.BSIM4rbodyMod)
    // {   (*(instance.BSIM4DPdbPtr) += gcdbdb - instance.BSIM4gbd);
    //     (*(instance.BSIM4SPsbPtr) -= instance.BSIM4gbs - gcsbsb);

    //     (*(instance.BSIM4DBdpPtr) += gcdbdb - instance.BSIM4gbd);
    //     (*(instance.BSIM4DBdbPtr) += instance.BSIM4gbd - gcdbdb
    //                             + instance.BSIM4grbpd + instance.BSIM4grbdb);
    //     (*(instance.BSIM4DBbpPtr) -= instance.BSIM4grbpd);
    //     (*(instance.BSIM4DBbPtr) -= instance.BSIM4grbdb);

    //     (*(instance.BSIM4BPdbPtr) -= instance.BSIM4grbpd);
    //     (*(instance.BSIM4BPbPtr) -= instance.BSIM4grbpb);
    //     (*(instance.BSIM4BPsbPtr) -= instance.BSIM4grbps);
    //     (*(instance.BSIM4BPbpPtr) += instance.BSIM4grbpd + instance.BSIM4grbps
    //                             + instance.BSIM4grbpb);
    // /* WDLiu: (gcbbb - instance.BSIM4gbbs) already added to BPbpPtr */

    //     (*(instance.BSIM4SBspPtr) += gcsbsb - instance.BSIM4gbs);
    //     (*(instance.BSIM4SBbpPtr) -= instance.BSIM4grbps);
    //     (*(instance.BSIM4SBbPtr) -= instance.BSIM4grbsb);
    //     (*(instance.BSIM4SBsbPtr) += instance.BSIM4gbs - gcsbsb
    //                             + instance.BSIM4grbps + instance.BSIM4grbsb);

    //     (*(instance.BSIM4BdbPtr) -= instance.BSIM4grbdb);
    //     (*(instance.BSIM4BbpPtr) -= instance.BSIM4grbpb);
    //     (*(instance.BSIM4BsbPtr) -= instance.BSIM4grbsb);
    //     (*(instance.BSIM4BbPtr) += instance.BSIM4grbsb + instance.BSIM4grbdb
    //                             + instance.BSIM4grbpb);
    // }

    // if (instance.BSIM4trnqsMod)
    // {   (*(instance.BSIM4QqPtr) += gqdef + instance.BSIM4gtau);
    //     (*(instance.BSIM4QgpPtr) += ggtg - gcqgb);
    //     (*(instance.BSIM4QdpPtr) += ggtd - gcqdb);
    //     (*(instance.BSIM4QspPtr) += ggts - gcqsb);
    //     (*(instance.BSIM4QbpPtr) += ggtb - gcqbb);

    //     (*(instance.BSIM4DPqPtr) += dxpart * instance.BSIM4gtau);
    //     (*(instance.BSIM4SPqPtr) += sxpart * instance.BSIM4gtau);
    //     (*(instance.BSIM4GPqPtr) -= instance.BSIM4gtau);
    // }

line1000:  ;


return(true);
}

void updateState1(BSIM4V82 &inst){
    inst.BSIM4states1 = inst.BSIM4states0;
}


} // namespace bsim4


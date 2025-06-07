/*
    b4acld.h - BSIM4v4.8.2
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
#include <armadillo>
#include <cmath>
#include "bsim4v82.hpp"
#include "CKT.hpp"

namespace bsim4{

// For safe modification of MNA matrix
/*  
    These checks are critical to prevent out-of-bounds access, 
    which could lead to runtime errors or undefined behavior, especially
    since negative indices are invalid for Armadillo matrices.
*/
inline void acStamp(arma::cx_dmat& mat, int row, int col, double real_val, double imag_val,
                    const std::array<bool, 12> &BSIM4nodeValid, BSIM4V82::NodeType row_node, BSIM4V82::NodeType col_node){
    if (BSIM4nodeValid[row_node] && BSIM4nodeValid[col_node]) {
        mat(row, col).real(mat(row, col).real() + real_val);
        mat(row, col).imag(mat(row, col).imag() + imag_val);
    }
}

int
BSIM4acLoad(const CKTcircuit &ckt, const BSIM4model &model, BSIM4V82 &instance, const SPICECompatible &spice, double omega,
            arma::cx_dmat &LHS, arma::cx_dmat &RHS)
{
double gjbd, gjbs, geltd, gcrg, gcrgg, gcrgd, gcrgs, gcrgb;
double xcbgb, xcbdb, xcbsb, xcbbb;
double xcggbr, xcgdbr, xcgsbr, xcgbbr, xcggbi, xcgdbi, xcgsbi, xcgbbi;
double Cggr, Cgdr, Cgsr, Cgbr, Cggi, Cgdi, Cgsi, Cgbi;
double xcddbr, xcdgbr, xcdsbr, xcdbbr, xcsdbr, xcsgbr, xcssbr, xcsbbr;
double xcddbi, xcdgbi, xcdsbi, xcdbbi, xcsdbi, xcsgbi, xcssbi, xcsbbi;
double xcdbdb, xcsbsb=0.0, xcgmgmb=0.0, xcgmdb=0.0, xcgmsb=0.0, xcdgmb, xcsgmb;
double xcgmbb=0.0, xcbgmb;
double capbd, capbs; // omega deleted here
double gstot, gstotd, gstotg, gstots, gstotb, gspr;
double gdtot, gdtotd, gdtotg, gdtots, gdtotb, gdpr;
double gIstotg, gIstotd, gIstots, gIstotb;
double gIdtotg, gIdtotd, gIdtots, gIdtotb;
double gIbtotg, gIbtotd, gIbtots, gIbtotb;
double gIgtotg, gIgtotd, gIgtots, gIgtotb;
double cgso, cgdo, cgbo;
double gbspsp, gbbdp, gbbsp, gbspg, gbspb;
double gbspdp, gbdpdp, gbdpg, gbdpb, gbdpsp;
double T0=0.0, T1, T2, T3;
double Csg, Csd, Css;
double Cdgr, Cddr, Cdsr, Cdbr, Csgr, Csdr, Cssr, Csbr;
double Cdgi, Cddi, Cdsi, Cdbi, Csgi, Csdi, Cssi, Csbi;
double gmr, gmi, gmbsr, gmbsi, gdsr, gdsi;
double FwdSumr, RevSumr, Gmr, Gmbsr;
double FwdSumi, RevSumi, Gmi, Gmbsi;
struct bsim4SizeDependParam pParam;
double ggidld, ggidlg, ggidlb, ggislg, ggislb, ggisls;

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

// double m;

    pParam = *instance.pParam;
    capbd = instance.BSIM4capbd;
    capbs = instance.BSIM4capbs;
    cgso = instance.BSIM4cgso;
    cgdo = instance.BSIM4cgdo;
    cgbo = pParam.BSIM4cgbo;

    Csd = -(instance.BSIM4cddb + instance.BSIM4cgdb + instance.BSIM4cbdb);
    Csg = -(instance.BSIM4cdgb + instance.BSIM4cggb + instance.BSIM4cbgb);
    Css = -(instance.BSIM4cdsb + instance.BSIM4cgsb + instance.BSIM4cbsb);

    if (instance.BSIM4acnqsMod)
    {   T0 = omega * instance.BSIM4taunet;
        T1 = T0 * T0;
        T2 = 1.0 / (1.0 + T1);
        T3 = T0 * T2;

        gmr = instance.BSIM4gm * T2;
        gmbsr = instance.BSIM4gmbs * T2;
        gdsr = instance.BSIM4gds * T2;

        gmi = -instance.BSIM4gm * T3;
        gmbsi = -instance.BSIM4gmbs * T3;
        gdsi = -instance.BSIM4gds * T3;

        Cddr = instance.BSIM4cddb * T2;
        Cdgr = instance.BSIM4cdgb * T2;
        Cdsr = instance.BSIM4cdsb * T2;
        Cdbr = -(Cddr + Cdgr + Cdsr);

        /* WDLiu: Cxyi mulitplied by jomega below, and actually to be of conductance */
        Cddi = instance.BSIM4cddb * T3 * omega;
        Cdgi = instance.BSIM4cdgb * T3 * omega;
        Cdsi = instance.BSIM4cdsb * T3 * omega;
        Cdbi = -(Cddi + Cdgi + Cdsi);

        Csdr = Csd * T2;
        Csgr = Csg * T2;
        Cssr = Css * T2;
        Csbr = -(Csdr + Csgr + Cssr);

        Csdi = Csd * T3 * omega;
        Csgi = Csg * T3 * omega;
        Cssi = Css * T3 * omega;
        Csbi = -(Csdi + Csgi + Cssi);

        Cgdr = -(Cddr + Csdr + instance.BSIM4cbdb);
        Cggr = -(Cdgr + Csgr + instance.BSIM4cbgb);
        Cgsr = -(Cdsr + Cssr + instance.BSIM4cbsb);
        Cgbr = -(Cgdr + Cggr + Cgsr);

        Cgdi = -(Cddi + Csdi);
        Cggi = -(Cdgi + Csgi);
        Cgsi = -(Cdsi + Cssi);
        Cgbi = -(Cgdi + Cggi + Cgsi);
    }
    else /* QS */
    {   gmr = instance.BSIM4gm;
        gmbsr = instance.BSIM4gmbs;
        gdsr = instance.BSIM4gds;
        gmi = gmbsi = gdsi = 0.0;

        Cddr = instance.BSIM4cddb;
        Cdgr = instance.BSIM4cdgb;
        Cdsr = instance.BSIM4cdsb;
        Cdbr = -(Cddr + Cdgr + Cdsr);
        Cddi = Cdgi = Cdsi = Cdbi = 0.0;

        Csdr = Csd;
        Csgr = Csg;
        Cssr = Css;
        Csbr = -(Csdr + Csgr + Cssr);
        Csdi = Csgi = Cssi = Csbi = 0.0;

        Cgdr = instance.BSIM4cgdb;
        Cggr = instance.BSIM4cggb;
        Cgsr = instance.BSIM4cgsb;
        Cgbr = -(Cgdr + Cggr + Cgsr);
        Cgdi = Cggi = Cgsi = Cgbi = 0.0;
    }


    if (instance.BSIM4mode >= 0) 
    {   Gmr = gmr;
        Gmbsr = gmbsr;
        FwdSumr = Gmr + Gmbsr;
        RevSumr = 0.0;
        Gmi = gmi;
        Gmbsi = gmbsi;
        FwdSumi = Gmi + Gmbsi;
        RevSumi = 0.0;

        gbbdp = -(instance.BSIM4gbds);
        gbbsp = instance.BSIM4gbds + instance.BSIM4gbgs + instance.BSIM4gbbs;
        gbdpg = instance.BSIM4gbgs;
        gbdpdp = instance.BSIM4gbds;
        gbdpb = instance.BSIM4gbbs;
        gbdpsp = -(gbdpg + gbdpdp + gbdpb);

        gbspdp = 0.0;
        gbspg = 0.0;
        gbspb = 0.0;
        gbspsp = 0.0;

        if (model.BSIM4igcMod)
        {   gIstotg = instance.BSIM4gIgsg + instance.BSIM4gIgcsg;
            gIstotd = instance.BSIM4gIgcsd;
            gIstots = instance.BSIM4gIgss + instance.BSIM4gIgcss;
            gIstotb = instance.BSIM4gIgcsb;

            gIdtotg = instance.BSIM4gIgdg + instance.BSIM4gIgcdg;
            gIdtotd = instance.BSIM4gIgdd + instance.BSIM4gIgcdd;
            gIdtots = instance.BSIM4gIgcds;
            gIdtotb = instance.BSIM4gIgcdb;
        }
        else
        {   gIstotg = gIstotd = gIstots = gIstotb = 0.0;
            gIdtotg = gIdtotd = gIdtots = gIdtotb = 0.0;
        }

        if (model.BSIM4igbMod)
        {   gIbtotg = instance.BSIM4gIgbg;
            gIbtotd = instance.BSIM4gIgbd;
            gIbtots = instance.BSIM4gIgbs;
            gIbtotb = instance.BSIM4gIgbb;
        }
        else
            gIbtotg = gIbtotd = gIbtots = gIbtotb = 0.0;

        if ((model.BSIM4igcMod != 0) || (model.BSIM4igbMod != 0))
        {   gIgtotg = gIstotg + gIdtotg + gIbtotg;
            gIgtotd = gIstotd + gIdtotd + gIbtotd ;
            gIgtots = gIstots + gIdtots + gIbtots;
            gIgtotb = gIstotb + gIdtotb + gIbtotb;
        }
        else
            gIgtotg = gIgtotd = gIgtots = gIgtotb = 0.0;

        if (instance.BSIM4rgateMod == 2)
            T0 = (instance.BSIM4states0[BSIM4vges])
                - (instance.BSIM4states0[BSIM4vgs]);
        else if (instance.BSIM4rgateMod == 3)
            T0 = (instance.BSIM4states0[BSIM4vgms])
                - (instance.BSIM4states0[BSIM4vgs]);
        if (instance.BSIM4rgateMod > 1)
        {   gcrgd = instance.BSIM4gcrgd * T0;
            gcrgg = instance.BSIM4gcrgg * T0;
            gcrgs = instance.BSIM4gcrgs * T0;
            gcrgb = instance.BSIM4gcrgb * T0;
            gcrgg -= instance.BSIM4gcrg;
            gcrg = instance.BSIM4gcrg;
        }
        else
            gcrg = gcrgd = gcrgg = gcrgs = gcrgb = 0.0;

        if (instance.BSIM4rgateMod == 3)
        {   xcgmgmb = (cgdo + cgso + pParam.BSIM4cgbo) * omega;
            xcgmdb = -cgdo * omega;
            xcgmsb = -cgso * omega;
            xcgmbb = -pParam.BSIM4cgbo * omega;

            xcdgmb = xcgmdb;
            xcsgmb = xcgmsb;
            xcbgmb = xcgmbb;

            xcggbr = Cggr * omega;
            xcgdbr = Cgdr * omega;
            xcgsbr = Cgsr * omega;
            xcgbbr = -(xcggbr + xcgdbr + xcgsbr);

            xcdgbr = Cdgr * omega;
            xcsgbr = Csgr * omega;
            xcbgb = instance.BSIM4cbgb * omega;
        }
        else
        {   xcggbr = (Cggr + cgdo + cgso + pParam.BSIM4cgbo ) * omega;
            xcgdbr = (Cgdr - cgdo) * omega;
            xcgsbr = (Cgsr - cgso) * omega;
            xcgbbr = -(xcggbr + xcgdbr + xcgsbr);

            xcdgbr = (Cdgr - cgdo) * omega;
            xcsgbr = (Csgr - cgso) * omega;
            xcbgb = (instance.BSIM4cbgb - pParam.BSIM4cgbo) * omega;

            xcdgmb = xcsgmb = xcbgmb = 0.0;
        }
        xcddbr = (Cddr + instance.BSIM4capbd + cgdo) * omega;
        xcdsbr = Cdsr * omega;
        xcsdbr = Csdr * omega;
        xcssbr = (instance.BSIM4capbs + cgso + Cssr) * omega;

        if (!instance.BSIM4rbodyMod)
        {   xcdbbr = -(xcdgbr + xcddbr + xcdsbr + xcdgmb);
            xcsbbr = -(xcsgbr + xcsdbr + xcssbr + xcsgmb);

            xcbdb = (instance.BSIM4cbdb - instance.BSIM4capbd) * omega;
            xcbsb = (instance.BSIM4cbsb - instance.BSIM4capbs) * omega;
            xcdbdb = 0.0;
        }
        else
        {   xcdbbr = Cdbr * omega;
            xcsbbr = -(xcsgbr + xcsdbr + xcssbr + xcsgmb)
                    + instance.BSIM4capbs * omega;

            xcbdb = instance.BSIM4cbdb * omega;
            xcbsb = instance.BSIM4cbsb * omega;

            xcdbdb = -instance.BSIM4capbd * omega;
            xcsbsb = -instance.BSIM4capbs * omega;
        }
        xcbbb = -(xcbdb + xcbgb + xcbsb + xcbgmb);

        xcdgbi = Cdgi;
        xcsgbi = Csgi;
        xcddbi = Cddi;
        xcdsbi = Cdsi;
        xcsdbi = Csdi;
        xcssbi = Cssi;
        xcdbbi = Cdbi;
        xcsbbi = Csbi;
        xcggbi = Cggi;
        xcgdbi = Cgdi;
        xcgsbi = Cgsi;
        xcgbbi = Cgbi;
    } 
    else /* Reverse mode */
    {   Gmr = -gmr;
        Gmbsr = -gmbsr;
        FwdSumr = 0.0;
        RevSumr = -(Gmr + Gmbsr);
        Gmi = -gmi;
        Gmbsi = -gmbsi;
        FwdSumi = 0.0;
        RevSumi = -(Gmi + Gmbsi);

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
        {   gIstotg = instance.BSIM4gIgsg + instance.BSIM4gIgcdg;
            gIstotd = instance.BSIM4gIgcds;
            gIstots = instance.BSIM4gIgss + instance.BSIM4gIgcdd;
            gIstotb = instance.BSIM4gIgcdb;

            gIdtotg = instance.BSIM4gIgdg + instance.BSIM4gIgcsg;
            gIdtotd = instance.BSIM4gIgdd + instance.BSIM4gIgcss;
            gIdtots = instance.BSIM4gIgcsd;
            gIdtotb = instance.BSIM4gIgcsb;
        }
        else
        {   gIstotg = gIstotd = gIstots = gIstotb = 0.0;
            gIdtotg = gIdtotd = gIdtots = gIdtotb  = 0.0;
        }

        if (model.BSIM4igbMod)
        {   gIbtotg = instance.BSIM4gIgbg;
            gIbtotd = instance.BSIM4gIgbs;
            gIbtots = instance.BSIM4gIgbd;
            gIbtotb = instance.BSIM4gIgbb;
        }
        else
            gIbtotg = gIbtotd = gIbtots = gIbtotb = 0.0;

        if ((model.BSIM4igcMod != 0) || (model.BSIM4igbMod != 0))
        {   gIgtotg = gIstotg + gIdtotg + gIbtotg;
            gIgtotd = gIstotd + gIdtotd + gIbtotd ;
            gIgtots = gIstots + gIdtots + gIbtots;
            gIgtotb = gIstotb + gIdtotb + gIbtotb;
        }
        else
            gIgtotg = gIgtotd = gIgtots = gIgtotb = 0.0;

        if (instance.BSIM4rgateMod == 2)
            T0 = (instance.BSIM4states0[BSIM4vges])
                - (instance.BSIM4states0[BSIM4vgs]);
        else if (instance.BSIM4rgateMod == 3)
            T0 = (instance.BSIM4states0[BSIM4vgms])
                - (instance.BSIM4states0[BSIM4vgs]);
        if (instance.BSIM4rgateMod > 1)
        {   gcrgd = instance.BSIM4gcrgs * T0;
            gcrgg = instance.BSIM4gcrgg * T0;
            gcrgs = instance.BSIM4gcrgd * T0;
            gcrgb = instance.BSIM4gcrgb * T0;
            gcrgg -= instance.BSIM4gcrg;
            gcrg = instance.BSIM4gcrg;
        }
        else
            gcrg = gcrgd = gcrgg = gcrgs = gcrgb = 0.0;

        if (instance.BSIM4rgateMod == 3)
        {   xcgmgmb = (cgdo + cgso + pParam.BSIM4cgbo) * omega;
            xcgmdb = -cgdo * omega;
            xcgmsb = -cgso * omega;
            xcgmbb = -pParam.BSIM4cgbo * omega;

            xcdgmb = xcgmdb;
            xcsgmb = xcgmsb;
            xcbgmb = xcgmbb;

            xcggbr = Cggr * omega;
            xcgdbr = Cgsr * omega;
            xcgsbr = Cgdr * omega;
            xcgbbr = -(xcggbr + xcgdbr + xcgsbr);

            xcdgbr = Csgr * omega;
            xcsgbr = Cdgr * omega;
            xcbgb = instance.BSIM4cbgb * omega;
        }
        else
        {   xcggbr = (Cggr + cgdo + cgso + pParam.BSIM4cgbo ) * omega;
            xcgdbr = (Cgsr - cgdo) * omega;
            xcgsbr = (Cgdr - cgso) * omega;
            xcgbbr = -(xcggbr + xcgdbr + xcgsbr);

            xcdgbr = (Csgr - cgdo) * omega;
            xcsgbr = (Cdgr - cgso) * omega;
            xcbgb = (instance.BSIM4cbgb - pParam.BSIM4cgbo) * omega;

            xcdgmb = xcsgmb = xcbgmb = 0.0;
        }
        xcddbr = (instance.BSIM4capbd + cgdo + Cssr) * omega;
        xcdsbr = Csdr * omega;
        xcsdbr = Cdsr * omega;
        xcssbr = (Cddr + instance.BSIM4capbs + cgso) * omega;

        if (!instance.BSIM4rbodyMod)
        {   xcdbbr = -(xcdgbr + xcddbr + xcdsbr + xcdgmb);
            xcsbbr = -(xcsgbr + xcsdbr + xcssbr + xcsgmb);

            xcbdb = (instance.BSIM4cbsb - instance.BSIM4capbd) * omega;
            xcbsb = (instance.BSIM4cbdb - instance.BSIM4capbs) * omega;
            xcdbdb = 0.0;
        }
        else
        {   xcdbbr = -(xcdgbr + xcddbr + xcdsbr + xcdgmb)
                    + instance.BSIM4capbd * omega;
            xcsbbr = Cdbr * omega;

            xcbdb = instance.BSIM4cbsb * omega;
            xcbsb = instance.BSIM4cbdb * omega;
            xcdbdb = -instance.BSIM4capbd * omega;
            xcsbsb = -instance.BSIM4capbs * omega;
        }
        xcbbb = -(xcbgb + xcbdb + xcbsb + xcbgmb);

        xcdgbi = Csgi;
        xcsgbi = Cdgi;
        xcddbi = Cssi;
        xcdsbi = Csdi;
        xcsdbi = Cdsi;
        xcssbi = Cddi;
        xcdbbi = Csbi;
        xcsbbi = Cdbi;
        xcggbi = Cggi;
        xcgdbi = Cgsi;
        xcgsbi = Cgdi;
        xcgbbi = Cgbi;
    }

    if (model.BSIM4rdsMod == 1)
    {   gstot = instance.BSIM4gstot;
        gstotd = instance.BSIM4gstotd;
        gstotg = instance.BSIM4gstotg;
        gstots = instance.BSIM4gstots - gstot;
        gstotb = instance.BSIM4gstotb;

        gdtot = instance.BSIM4gdtot;
        gdtotd = instance.BSIM4gdtotd - gdtot;
        gdtotg = instance.BSIM4gdtotg;
        gdtots = instance.BSIM4gdtots;
        gdtotb = instance.BSIM4gdtotb;
    }
    else
    {   gstot = gstotd = gstotg = gstots = gstotb = 0.0;
        gdtot = gdtotd = gdtotg = gdtots = gdtotb = 0.0;
    }


    /*
    * Loading AC matrix
    */
    // m = instance.BSIM4m;

    // MNA Matrix Stamping
    // Real part, Different nodes are distinguished by uppercase and lowercase letters
    double  BSIM4DPbpReal, BSIM4GPbpReal, BSIM4SPbpReal, BSIM4BPdpReal, BSIM4BPgpReal, BSIM4BPspReal, BSIM4BPbpReal, BSIM4DdReal, BSIM4GPgpReal, BSIM4SsReal, 
        BSIM4DPdpReal, BSIM4SPspReal, BSIM4DdpReal, BSIM4GPdpReal, BSIM4GPspReal, BSIM4SspReal, BSIM4DPspReal, BSIM4DPdReal, BSIM4DPgpReal, BSIM4SPgpReal,
        BSIM4SPsReal, BSIM4SPdpReal, BSIM4QqReal, BSIM4QbpReal, BSIM4QdpReal, BSIM4QspReal, BSIM4QgpReal, BSIM4DPqReal, BSIM4SPqReal, BSIM4GPqReal,
        BSIM4GEgeReal, BSIM4GEgpReal, BSIM4GPgeReal, BSIM4GEdpReal, BSIM4GEspReal, BSIM4GEbpReal, BSIM4GMdpReal, BSIM4GMgpReal, BSIM4GMgmReal, BSIM4GMgeReal,
        BSIM4GMspReal, BSIM4GMbpReal, BSIM4DPgmReal, BSIM4GPgmReal, BSIM4GEgmReal, BSIM4SPgmReal, BSIM4BPgmReal, BSIM4DPdbReal, BSIM4SPsbReal, BSIM4DBdpReal, 
        BSIM4DBdbReal, BSIM4DBbpReal, BSIM4DBbReal, BSIM4BPdbReal, BSIM4BPbReal, BSIM4BPsbReal, BSIM4SBspReal, BSIM4SBbpReal, BSIM4SBbReal, BSIM4SBsbReal, 
        BSIM4BdbReal, BSIM4BbpReal, BSIM4BsbReal, BSIM4BbReal, BSIM4DgpReal, BSIM4DspReal, BSIM4DbpReal, BSIM4SdpReal, BSIM4SgpReal, BSIM4SbpReal;
    // Imaginary part
    double  BSIM4DPbpImag, BSIM4GPbpImag, BSIM4SPbpImag, BSIM4BPdpImag, BSIM4BPgpImag, BSIM4BPspImag, BSIM4BPbpImag, BSIM4DdImag, BSIM4GPgpImag, BSIM4SsImag, 
        BSIM4DPdpImag, BSIM4SPspImag, BSIM4DdpImag, BSIM4GPdpImag, BSIM4GPspImag, BSIM4SspImag, BSIM4DPspImag, BSIM4DPdImag, BSIM4DPgpImag, BSIM4SPgpImag,
        BSIM4SPsImag, BSIM4SPdpImag, BSIM4QqImag, BSIM4QbpImag, BSIM4QdpImag, BSIM4QspImag, BSIM4QgpImag, BSIM4DPqImag, BSIM4SPqImag, BSIM4GPqImag,
        BSIM4GEgeImag, BSIM4GEgpImag, BSIM4GPgeImag, BSIM4GEdpImag, BSIM4GEspImag, BSIM4GEbpImag, BSIM4GMdpImag, BSIM4GMgpImag, BSIM4GMgmImag, BSIM4GMgeImag,
        BSIM4GMspImag, BSIM4GMbpImag, BSIM4DPgmImag, BSIM4GPgmImag, BSIM4GEgmImag, BSIM4SPgmImag, BSIM4BPgmImag, BSIM4DPdbImag, BSIM4SPsbImag, BSIM4DBdpImag, 
        BSIM4DBdbImag, BSIM4DBbpImag, BSIM4DBbImag, BSIM4BPdbImag, BSIM4BPbImag, BSIM4BPsbImag, BSIM4SBspImag, BSIM4SBbpImag, BSIM4SBbImag, BSIM4SBsbImag, 
        BSIM4BdbImag, BSIM4BbpImag, BSIM4BsbImag, BSIM4BbImag, BSIM4DgpImag, BSIM4DspImag, BSIM4DbpImag, BSIM4SdpImag, BSIM4SgpImag, BSIM4SbpImag;

    if (!model.BSIM4rdsMod)
    {   gdpr = instance.BSIM4drainConductance;
        gspr = instance.BSIM4sourceConductance;
    }
    else
        gdpr = gspr = 0.0;

    if (!instance.BSIM4rbodyMod)
    {   gjbd = instance.BSIM4gbd;
        gjbs = instance.BSIM4gbs;
    }
    else
        gjbd = gjbs = 0.0;

    geltd = instance.BSIM4grgeltd;

    if (instance.BSIM4rgateMod == 1)
    {   (BSIM4GEgeReal) +=  geltd;
        (BSIM4GPgeReal) -=  geltd;
        (BSIM4GEgpReal) -=  geltd;

        (BSIM4GPgpImag) +=  xcggbr;    // BSIM4GPgpPtr + 1
        (BSIM4GPgpReal) +=  (geltd + xcggbi + gIgtotg);
        (BSIM4GPdpImag) +=  xcgdbr;    // BSIM4GPdpPtr + 1
        (BSIM4GPdpReal) +=  (xcgdbi + gIgtotd);
        (BSIM4GPspImag) +=  xcgsbr;    // BSIM4GPspPtr + 1
        (BSIM4GPspReal) +=  (xcgsbi + gIgtots);
        (BSIM4GPbpImag) +=  xcgbbr;    // BSIM4GPbpPtr + 1
        (BSIM4GPbpReal) +=  (xcgbbi + gIgtotb);
    } /* WDLiu: gcrg already subtracted from all gcrgg below */
    else if (instance.BSIM4rgateMod == 2)
    {   (BSIM4GEgeReal) +=  gcrg;
        (BSIM4GEgpReal) +=  gcrgg;
        (BSIM4GEdpReal) +=  gcrgd;
        (BSIM4GEspReal) +=  gcrgs;
        (BSIM4GEbpReal) +=  gcrgb;

        (BSIM4GPgeReal) -=  gcrg;
        (BSIM4GPgpImag) +=  xcggbr;    // BSIM4GPgpPtr + 1
        (BSIM4GPgpReal) -=  (gcrgg - xcggbi - gIgtotg);
        (BSIM4GPdpImag) +=  xcgdbr;    // BSIM4GPdpPtr + 1
        (BSIM4GPdpReal) -=  (gcrgd - xcgdbi - gIgtotd);
        (BSIM4GPspImag) +=  xcgsbr;    // BSIM4GPspPtr + 1
        (BSIM4GPspReal) -=  (gcrgs - xcgsbi - gIgtots);
        (BSIM4GPbpImag) +=  xcgbbr;    // BSIM4GPbpPtr + 1
        (BSIM4GPbpReal) -=  (gcrgb - xcgbbi - gIgtotb);
    }
    else if (instance.BSIM4rgateMod == 3)
    {   (BSIM4GEgeReal) +=  geltd;
        (BSIM4GEgmReal) -=  geltd;
        (BSIM4GMgeReal) -=  geltd;
        (BSIM4GMgmReal) +=  (geltd + gcrg);
        (BSIM4GMgmImag) +=  xcgmgmb;   // BSIM4GMgmPtr + 1

        (BSIM4GMdpReal) +=  gcrgd;
        (BSIM4GMdpImag) +=  xcgmdb;    // BSIM4GMdpPtr + 1
        (BSIM4GMgpReal) +=  gcrgg;
        (BSIM4GMspReal) +=  gcrgs;
        (BSIM4GMspImag) +=  xcgmsb;    // BSIM4GMspPtr + 1
        (BSIM4GMbpReal) +=  gcrgb;
        (BSIM4GMbpImag) +=  xcgmbb;    // BSIM4GMbpPtr + 1

        (BSIM4DPgmImag) +=  xcdgmb;    // BSIM4DPgmPtr + 1
        (BSIM4GPgmReal) -=  gcrg;
        (BSIM4SPgmImag) +=  xcsgmb;    // BSIM4SPgmPtr + 1
        (BSIM4BPgmImag) +=  xcbgmb;    // BSIM4BPgmPtr + 1

        (BSIM4GPgpReal) -=  (gcrgg - xcggbi - gIgtotg);
        (BSIM4GPgpImag) +=  xcggbr;    // BSIM4GPgpPtr + 1
        (BSIM4GPdpReal) -=  (gcrgd - xcgdbi - gIgtotd);
        (BSIM4GPdpImag) +=  xcgdbr;    // BSIM4GPdpPtr + 1
        (BSIM4GPspReal) -=  (gcrgs - xcgsbi - gIgtots);
        (BSIM4GPspImag) +=  xcgsbr;    // BSIM4GPspPtr + 1
        (BSIM4GPbpReal) -=  (gcrgb - xcgbbi - gIgtotb);
        (BSIM4GPbpImag) +=  xcgbbr;    // BSIM4GPbpPtr + 1
    }
    else
    {   (BSIM4GPgpImag) +=  xcggbr;    // BSIM4GPgpPtr + 1
        (BSIM4GPgpReal) +=  (xcggbi + gIgtotg);
        (BSIM4GPdpImag) +=  xcgdbr;    // BSIM4GPdpPtr + 1
        (BSIM4GPdpReal) +=  (xcgdbi + gIgtotd);
        (BSIM4GPspImag) +=  xcgsbr;    // BSIM4GPspPtr + 1
        (BSIM4GPspReal) +=  (xcgsbi + gIgtots);
        (BSIM4GPbpImag) +=  xcgbbr;    // BSIM4GPbpPtr + 1
        (BSIM4GPbpReal) +=  (xcgbbi + gIgtotb);
    }

    if (model.BSIM4rdsMod)
    {   ((BSIM4DgpReal) +=  gdtotg);
        ((BSIM4DspReal) +=  gdtots);
        ((BSIM4DbpReal) +=  gdtotb);
        ((BSIM4SdpReal) +=  gstotd);
        ((BSIM4SgpReal) +=  gstotg);
        ((BSIM4SbpReal) +=  gstotb);
    }

    (BSIM4DPdpImag) +=  (xcddbr + gdsi + RevSumi);     // BSIM4DPdpPtr + 1
    (BSIM4DPdpReal) +=  (gdpr + xcddbi + gdsr + instance.BSIM4gbd 
                            - gdtotd + RevSumr + gbdpdp - gIdtotd);
    (BSIM4DPdReal) -=  (gdpr + gdtot);
    (BSIM4DPgpImag) +=  (xcdgbr + Gmi);    // BSIM4DPgpPtr + 1
    (BSIM4DPgpReal) +=  (Gmr + xcdgbi - gdtotg + gbdpg - gIdtotg);
    (BSIM4DPspImag) +=  (xcdsbr - gdsi - FwdSumi);   // BSIM4DPspPtr + 1
    (BSIM4DPspReal) -=  (gdsr - xcdsbi + FwdSumr + gdtots - gbdpsp + gIdtots);
    (BSIM4DPbpImag) +=  (xcdbbr + Gmbsi);      // BSIM4DPbpPtr + 1
    (BSIM4DPbpReal) -=  (gjbd + gdtotb - xcdbbi - Gmbsr - gbdpb + gIdtotb);

    (BSIM4DdpReal) -=  (gdpr - gdtotd);
    (BSIM4DdReal) +=  (gdpr + gdtot);

    (BSIM4SPdpImag) +=  (xcsdbr - gdsi - RevSumi);   // BSIM4SPdpPtr + 1
    (BSIM4SPdpReal) -=  (gdsr - xcsdbi + gstotd + RevSumr - gbspdp + gIstotd);
    (BSIM4SPgpImag) +=  (xcsgbr - Gmi);    // BSIM4SPgpPtr + 1
    (BSIM4SPgpReal) -=  (Gmr - xcsgbi + gstotg - gbspg + gIstotg);
    (BSIM4SPspImag) +=  (xcssbr + gdsi + FwdSumi);     // BSIM4SPspPtr + 1
    (BSIM4SPspReal) +=  (gspr + xcssbi + gdsr + instance.BSIM4gbs
                            - gstots + FwdSumr + gbspsp - gIstots);
    (BSIM4SPsReal) -=  (gspr + gstot);
    (BSIM4SPbpImag) +=  (xcsbbr - Gmbsi);   // BSIM4SPbpPtr + 1
    (BSIM4SPbpReal) -=  (gjbs + gstotb - xcsbbi + Gmbsr - gbspb + gIstotb);

    (BSIM4SspReal) -=  (gspr - gstots);
    (BSIM4SsReal) +=  (gspr + gstot);

    (BSIM4BPdpImag) +=  xcbdb;     // BSIM4BPdpPtr + 1
    (BSIM4BPdpReal) -=  (gjbd - gbbdp + gIbtotd);
    (BSIM4BPgpImag) +=  xcbgb;     // BSIM4BPgpPtr + 1
    (BSIM4BPgpReal) -=  (instance.BSIM4gbgs + gIbtotg);
    (BSIM4BPspImag) +=  xcbsb;     // BSIM4BPspPtr + 1
    (BSIM4BPspReal) -=  (gjbs - gbbsp + gIbtots);
    (BSIM4BPbpImag) +=  xcbbb;     // BSIM4BPbpPtr + 1
    (BSIM4BPbpReal) +=  (gjbd + gjbs - instance.BSIM4gbbs
                            - gIbtotb);
    ggidld = instance.BSIM4ggidld;
    ggidlg = instance.BSIM4ggidlg;
    ggidlb = instance.BSIM4ggidlb;
    ggislg = instance.BSIM4ggislg;
    ggisls = instance.BSIM4ggisls;
    ggislb = instance.BSIM4ggislb;

    /* stamp gidl */
    (BSIM4DPdpReal) +=  ggidld;
    (BSIM4DPgpReal) +=  ggidlg;
    (BSIM4DPspReal) -=  ((ggidlg + ggidld) + ggidlb);
    (BSIM4DPbpReal) +=  ggidlb;
    (BSIM4BPdpReal) -=  ggidld;
    (BSIM4BPgpReal) -=  ggidlg;
    (BSIM4BPspReal) +=  ((ggidlg + ggidld) + ggidlb);
    (BSIM4BPbpReal) -=  ggidlb;
    /* stamp gisl */
    (BSIM4SPdpReal) -=  ((ggisls + ggislg) + ggislb);
    (BSIM4SPgpReal) +=  ggislg;
    (BSIM4SPspReal) +=  ggisls;
    (BSIM4SPbpReal) +=  ggislb;
    (BSIM4BPdpReal) +=  ((ggislg + ggisls) + ggislb);
    (BSIM4BPgpReal) -=  ggislg;
    (BSIM4BPspReal) -=  ggisls;
    (BSIM4BPbpReal) -=  ggislb;

    if (instance.BSIM4rbodyMod)
    {   (BSIM4DPdbImag) +=  xcdbdb;      // BSIM4DPdbPtr + 1
        (BSIM4DPdbReal) -=  instance.BSIM4gbd;
        (BSIM4SPsbImag) +=  xcsbsb;      // BSIM4SPsbPtr + 1
        (BSIM4SPsbReal) -=  instance.BSIM4gbs;

        (BSIM4DBdpImag) +=  xcdbdb;      // BSIM4DBdpPtr + 1
        (BSIM4DBdpReal) -=  instance.BSIM4gbd;
        (BSIM4DBdbImag) -=  xcdbdb;      // BSIM4DBdbPtr + 1
        (BSIM4DBdbReal) +=  (instance.BSIM4gbd + instance.BSIM4grbpd 
                                + instance.BSIM4grbdb);
        (BSIM4DBbpReal) -=  instance.BSIM4grbpd;
        (BSIM4DBbReal) -=  instance.BSIM4grbdb;

        (BSIM4BPdbReal) -=  instance.BSIM4grbpd;
        (BSIM4BPbReal) -=  instance.BSIM4grbpb;
        (BSIM4BPsbReal) -=  instance.BSIM4grbps;
        (BSIM4BPbpReal) +=  (instance.BSIM4grbpd + instance.BSIM4grbps 
                                + instance.BSIM4grbpb);
        /* WDLiu: (-instance.BSIM4gbbs) already added to BPbpPtr */

        (BSIM4SBspImag) +=  xcsbsb;      // BSIM4SBspPtr + 1
        (BSIM4SBspReal) -=  instance.BSIM4gbs;
        (BSIM4SBbpReal) -=  instance.BSIM4grbps;
        (BSIM4SBbReal) -=  instance.BSIM4grbsb;
        (BSIM4SBsbImag) -=  xcsbsb;      // BSIM4SBsbPtr + 1
        (BSIM4SBsbReal) +=  (instance.BSIM4gbs
                                + instance.BSIM4grbps + instance.BSIM4grbsb);

        (BSIM4BdbReal) -=  instance.BSIM4grbdb;
        (BSIM4BbpReal) -=  instance.BSIM4grbpb;
        (BSIM4BsbReal) -=  instance.BSIM4grbsb;
        (BSIM4BbReal) +=  (instance.BSIM4grbsb + instance.BSIM4grbdb
                            + instance.BSIM4grbpb);
    }

    /*
    * WDLiu: The internal charge node generated for transient NQS is not needed for
    *        AC NQS. The following is not doing a real job, but we have to keep it;
    *        otherwise a singular AC NQS matrix may occur if the transient NQS is on.
    *        The charge node is isolated from the instance.
    */
    if (instance.BSIM4trnqsMod)
    {   ((BSIM4QqReal) +=  1.0);
        ((BSIM4QgpReal) += 0.0);
        ((BSIM4QdpReal) += 0.0);
        ((BSIM4QspReal) += 0.0);
        ((BSIM4QbpReal) += 0.0);

        ((BSIM4DPqReal) += 0.0);
        ((BSIM4SPqReal) += 0.0);
        ((BSIM4GPqReal) += 0.0);
    }

    // Now stamping the real and imaginary parts to the MNA matrix
    const std::array<bool,12> &nodeValid = instance.BSIM4nodeValid;

    // BSIM4DPbpPtr
    acStamp(LHS, dNodePrime, bNodePrime, 
        BSIM4DPbpReal, BSIM4DPbpImag, nodeValid, 
        BSIM4V82::D_NODE_PRIME, BSIM4V82::B_NODE_PRIME);
    // BSIM4GPbpPtr
    acStamp(LHS, gNodePrime, bNodePrime, 
        BSIM4GPbpReal, BSIM4GPbpImag, nodeValid, 
        BSIM4V82::G_NODE_PRIME, BSIM4V82::B_NODE_PRIME);
    // BSIM4SPbpPtr
    acStamp(LHS, sNodePrime, bNodePrime, 
        BSIM4SPbpReal, BSIM4SPbpImag, nodeValid, 
        BSIM4V82::S_NODE_PRIME, BSIM4V82::B_NODE_PRIME);


    // BSIM4BPdpPtr
    acStamp(LHS, dNodePrime, bNodePrime, 
        BSIM4BPdpReal, BSIM4BPdpImag, nodeValid, 
        BSIM4V82::D_NODE_PRIME, BSIM4V82::B_NODE_PRIME);
    // BSIM4BPgpPtr
    acStamp(LHS, bNodePrime, gNodePrime, 
        BSIM4BPgpReal, BSIM4BPgpImag, nodeValid, 
        BSIM4V82::B_NODE_PRIME, BSIM4V82::G_NODE_PRIME);
    // BSIM4BPspPtr
    acStamp(LHS, bNodePrime, sNodePrime, 
        BSIM4BPspReal, BSIM4BPspImag, nodeValid, 
        BSIM4V82::B_NODE_PRIME, BSIM4V82::S_NODE_PRIME);
    // BSIM4BPbpPtr
    acStamp(LHS, bNodePrime, bNodePrime, 
        BSIM4BPbpReal, BSIM4BPbpImag, nodeValid, 
        BSIM4V82::B_NODE_PRIME, BSIM4V82::B_NODE_PRIME);


    // BSIM4DdPtr
    acStamp(LHS, dNode, dNode, 
        BSIM4DdReal, BSIM4DdImag, nodeValid, 
        BSIM4V82::D_NODE, BSIM4V82::D_NODE);
    // BSIM4GPgpPtr
    acStamp(LHS, gNodePrime, gNodePrime, 
        BSIM4GPgpReal, BSIM4GPgpImag, nodeValid, 
        BSIM4V82::G_NODE_PRIME, BSIM4V82::G_NODE_PRIME);
    // BSIM4SsPtr
    acStamp(LHS, sNode, sNode, 
        BSIM4SsReal, BSIM4SsImag, nodeValid, 
        BSIM4V82::S_NODE, BSIM4V82::S_NODE);
    // BSIM4DPdpPtr
    acStamp(LHS, dNodePrime, dNodePrime, 
        BSIM4DPdpReal, BSIM4DPdpImag, nodeValid, 
        BSIM4V82::D_NODE_PRIME, BSIM4V82::D_NODE_PRIME);
    // BSIM4SPspPtr
    acStamp(LHS, sNodePrime, sNodePrime, 
        BSIM4SPspReal, BSIM4SPspImag, nodeValid, 
        BSIM4V82::S_NODE_PRIME, BSIM4V82::S_NODE_PRIME);
    // BSIM4DdpPtr
    acStamp(LHS, dNode, dNodePrime, 
        BSIM4DdpReal, BSIM4DdpImag, nodeValid, 
        BSIM4V82::D_NODE, BSIM4V82::D_NODE_PRIME);
    // BSIM4GPdpPtr
    acStamp(LHS, gNodePrime, dNodePrime, 
        BSIM4GPdpReal, BSIM4GPdpImag, nodeValid, 
        BSIM4V82::G_NODE_PRIME, BSIM4V82::D_NODE_PRIME);
    // BSIM4GPspPtr
    acStamp(LHS, gNodePrime, sNodePrime, 
        BSIM4GPspReal, BSIM4GPspImag, nodeValid, 
        BSIM4V82::G_NODE_PRIME, BSIM4V82::S_NODE_PRIME);
    // BSIM4SspPtr
    acStamp(LHS, sNode, sNodePrime, 
        BSIM4SspReal, BSIM4SspImag, nodeValid, 
        BSIM4V82::S_NODE, BSIM4V82::S_NODE_PRIME);
    // BSIM4DPspPtr
    acStamp(LHS, dNodePrime, sNodePrime, 
        BSIM4DPspReal, BSIM4DPspImag, nodeValid, 
        BSIM4V82::D_NODE_PRIME, BSIM4V82::S_NODE_PRIME);
    // BSIM4DPdPtr
    acStamp(LHS, dNodePrime, dNode, 
        BSIM4DPdReal, BSIM4DPdImag, nodeValid, 
        BSIM4V82::D_NODE_PRIME, BSIM4V82::D_NODE);
    // BSIM4DPgpPtr
    acStamp(LHS, dNodePrime, gNodePrime, 
        BSIM4DPgpReal, BSIM4DPgpImag, nodeValid, 
        BSIM4V82::D_NODE_PRIME, BSIM4V82::G_NODE_PRIME);
    // BSIM4SPgpPtr
    acStamp(LHS, sNodePrime, gNodePrime, 
        BSIM4SPgpReal, BSIM4SPgpImag, nodeValid, 
        BSIM4V82::S_NODE_PRIME, BSIM4V82::G_NODE_PRIME);
    // BSIM4SPsPtr
    acStamp(LHS, sNodePrime, sNode, 
        BSIM4SPsReal, BSIM4SPsImag, nodeValid, 
        BSIM4V82::S_NODE_PRIME, BSIM4V82::S_NODE);
    // BSIM4SPdpPtr
    acStamp(LHS, sNodePrime, dNodePrime, 
        BSIM4SPdpReal, BSIM4SPdpImag, nodeValid, 
        BSIM4V82::S_NODE_PRIME, BSIM4V82::D_NODE_PRIME);


    // BSIM4QqPtr
    acStamp(LHS, qNode, qNode, 
        BSIM4QqReal, BSIM4QqImag, nodeValid, 
        BSIM4V82::Q_NODE, BSIM4V82::Q_NODE);
    // BSIM4QbpPtr
    acStamp(LHS, qNode, bNodePrime, 
        BSIM4QbpReal, BSIM4QbpImag, nodeValid, 
        BSIM4V82::Q_NODE, BSIM4V82::B_NODE_PRIME);
    // BSIM4QdpPtr
    acStamp(LHS, qNode, dNodePrime, 
        BSIM4QdpReal, BSIM4QdpImag, nodeValid, 
        BSIM4V82::Q_NODE, BSIM4V82::D_NODE_PRIME);
    // BSIM4QspPtr
    acStamp(LHS, qNode, sNodePrime, 
        BSIM4QspReal, BSIM4QspImag, nodeValid, 
        BSIM4V82::Q_NODE, BSIM4V82::S_NODE_PRIME);
    // BSIM4QgpPtr
    acStamp(LHS, qNode, gNodePrime, 
        BSIM4QgpReal, BSIM4QgpImag, nodeValid, 
        BSIM4V82::Q_NODE, BSIM4V82::G_NODE_PRIME);
    // BSIM4DPqPtr
    acStamp(LHS, dNodePrime, qNode, 
        BSIM4DPqReal, BSIM4DPqImag, nodeValid, 
        BSIM4V82::D_NODE_PRIME, BSIM4V82::Q_NODE);
    // BSIM4SPqPtr
    acStamp(LHS, sNodePrime, qNode, 
        BSIM4SPqReal, BSIM4SPqImag, nodeValid, 
        BSIM4V82::S_NODE_PRIME, BSIM4V82::Q_NODE);
    // BSIM4GPqPtr
    acStamp(LHS, gNodePrime, qNode, 
        BSIM4GPqReal, BSIM4GPqImag, nodeValid, 
        BSIM4V82::G_NODE_PRIME, BSIM4V82::Q_NODE);
        
    // if (inst.BSIM4rgateMod != 0)
        // BSIM4GEgePtr
        // BSIM4GEgpPtr
        // BSIM4GPgePtr
        // BSIM4GEdpPtr
        // BSIM4GEspPtr
        // BSIM4GEbpPtr
        // BSIM4GMdpPtr
        // BSIM4GMgpPtr
        // BSIM4GMgmPtr
        // BSIM4GMgePtr
        // BSIM4GMspPtr
        // BSIM4GMbpPtr
        // BSIM4DPgmPtr
        // BSIM4GPgmPtr
        // BSIM4GEgmPtr
        // BSIM4SPgmPtr
        // BSIM4BPgmPtr
    
    // if ((inst.BSIM4rbodyMod ==1) || (inst.BSIM4rbodyMod ==2))
        // BSIM4DPdbPtr
        // BSIM4SPsbPtr
        // BSIM4DBdpPtr
        // BSIM4DBdbPtr
        // BSIM4DBbpPtr
        // BSIM4DBbPtr
        // BSIM4BPdbPtr
        // BSIM4BPbPtr
        // BSIM4BPsbPtr
        // BSIM4SBspPtr
        // BSIM4SBbpPtr
        // BSIM4SBbPtr
        // BSIM4SBsbPtr
        // BSIM4BdbPtr
        // BSIM4BbpPtr
        // BSIM4BsbPtr
        // BSIM4BbPtr
    
    // if (model.BSIM4rdsMod)
        // BSIM4DgpPtr
        // BSIM4DspPtr
        // BSIM4DbpPtr
        // BSIM4SdpPtr
        // BSIM4SgpPtr
        // BSIM4SbpPtr

    return 0;
}

} // namespace bsim4
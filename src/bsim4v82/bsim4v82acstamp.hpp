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

    // if (instance.BSIM4rgateMod == 1)
    // {   *(instance.BSIM4GEgePtr) +=  geltd;
    //     *(instance.BSIM4GPgePtr) -=  geltd;
    //     *(instance.BSIM4GEgpPtr) -=  geltd;

    //     *(instance.BSIM4GPgpPtr +1) +=  xcggbr;
    //     *(instance.BSIM4GPgpPtr) +=  (geltd + xcggbi + gIgtotg);
    //     *(instance.BSIM4GPdpPtr +1) +=  xcgdbr;
    //     *(instance.BSIM4GPdpPtr) +=  (xcgdbi + gIgtotd);
    //     *(instance.BSIM4GPspPtr +1) +=  xcgsbr;
    //     *(instance.BSIM4GPspPtr) +=  (xcgsbi + gIgtots);
    //     *(instance.BSIM4GPbpPtr +1) +=  xcgbbr;
    //     *(instance.BSIM4GPbpPtr) +=  (xcgbbi + gIgtotb);
    // } /* WDLiu: gcrg already subtracted from all gcrgg below */
    // else if (instance.BSIM4rgateMod == 2)
    // {   *(instance.BSIM4GEgePtr) +=  gcrg;
    //     *(instance.BSIM4GEgpPtr) +=  gcrgg;
    //     *(instance.BSIM4GEdpPtr) +=  gcrgd;
    //     *(instance.BSIM4GEspPtr) +=  gcrgs;
    //     *(instance.BSIM4GEbpPtr) +=  gcrgb;

    //     *(instance.BSIM4GPgePtr) -=  gcrg;
    //     *(instance.BSIM4GPgpPtr +1) +=  xcggbr;
    //     *(instance.BSIM4GPgpPtr) -=  (gcrgg - xcggbi - gIgtotg);
    //     *(instance.BSIM4GPdpPtr +1) +=  xcgdbr;
    //     *(instance.BSIM4GPdpPtr) -=  (gcrgd - xcgdbi - gIgtotd);
    //     *(instance.BSIM4GPspPtr +1) +=  xcgsbr;
    //     *(instance.BSIM4GPspPtr) -=  (gcrgs - xcgsbi - gIgtots);
    //     *(instance.BSIM4GPbpPtr +1) +=  xcgbbr;
    //     *(instance.BSIM4GPbpPtr) -=  (gcrgb - xcgbbi - gIgtotb);
    // }
    // else if (instance.BSIM4rgateMod == 3)
    // {   *(instance.BSIM4GEgePtr) +=  geltd;
    //     *(instance.BSIM4GEgmPtr) -=  geltd;
    //     *(instance.BSIM4GMgePtr) -=  geltd;
    //     *(instance.BSIM4GMgmPtr) +=  (geltd + gcrg);
    //     *(instance.BSIM4GMgmPtr +1) +=  xcgmgmb;

    //     *(instance.BSIM4GMdpPtr) +=  gcrgd;
    //     *(instance.BSIM4GMdpPtr +1) +=  xcgmdb;
    //     *(instance.BSIM4GMgpPtr) +=  gcrgg;
    //     *(instance.BSIM4GMspPtr) +=  gcrgs;
    //     *(instance.BSIM4GMspPtr +1) +=  xcgmsb;
    //     *(instance.BSIM4GMbpPtr) +=  gcrgb;
    //     *(instance.BSIM4GMbpPtr +1) +=  xcgmbb;

    //     *(instance.BSIM4DPgmPtr +1) +=  xcdgmb;
    //     *(instance.BSIM4GPgmPtr) -=  gcrg;
    //     *(instance.BSIM4SPgmPtr +1) +=  xcsgmb;
    //     *(instance.BSIM4BPgmPtr +1) +=  xcbgmb;

    //     *(instance.BSIM4GPgpPtr) -=  (gcrgg - xcggbi - gIgtotg);
    //     *(instance.BSIM4GPgpPtr +1) +=  xcggbr;
    //     *(instance.BSIM4GPdpPtr) -=  (gcrgd - xcgdbi - gIgtotd);
    //     *(instance.BSIM4GPdpPtr +1) +=  xcgdbr;
    //     *(instance.BSIM4GPspPtr) -=  (gcrgs - xcgsbi - gIgtots);
    //     *(instance.BSIM4GPspPtr +1) +=  xcgsbr;
    //     *(instance.BSIM4GPbpPtr) -=  (gcrgb - xcgbbi - gIgtotb);
    //     *(instance.BSIM4GPbpPtr +1) +=  xcgbbr;
    // }
    // else
    // {   
        // BSIM4GPgpPtr
        if(!(gNodePrime < 0)){
            arma::cx_double &BSIM4GPgpPtr = LHS(gNodePrime, gNodePrime);
            BSIM4GPgpPtr.imag(BSIM4GPgpPtr.imag() + xcggbr); // BSIM4GPgpPtr +1
            BSIM4GPgpPtr.real(BSIM4GPgpPtr.real() + (xcggbi + gIgtotg)); // BSIM4GPgpPtr
        }
        // BSIM4GPdpPtr
        if(!(dNodePrime < 0) && !(gNodePrime < 0)){
            arma::cx_double &BSIM4GPdpPtr = LHS(gNodePrime, dNodePrime);
            BSIM4GPdpPtr.imag(BSIM4GPdpPtr.imag() + xcgdbr); // BSIM4GPdpPtr +1
            BSIM4GPdpPtr.real(BSIM4GPdpPtr.real() + (xcgdbi + gIgtotd)); // BSIM4GPdpPtr
        }
        // BSIM4GPspPtr
        if(!(sNodePrime < 0) && !(gNodePrime < 0)){
            arma::cx_double &BSIM4GPspPtr =  LHS(gNodePrime, sNodePrime);
            BSIM4GPspPtr.imag(BSIM4GPspPtr.imag() + xcgsbr); // BSIM4GPspPtr +1
            BSIM4GPspPtr.real(BSIM4GPspPtr.real() + (xcgsbi + gIgtots)); // BSIM4GPspPtr
        }
        // BSIM4GPbpPtr
        if(!(bNodePrime < 0) && !(gNodePrime < 0)){
            arma::cx_double &BSIM4GPbpPtr = LHS(gNodePrime, bNodePrime);
            BSIM4GPbpPtr.imag(BSIM4GPbpPtr.imag() + xcgbbr); // BSIM4GPbpPtr +1
            BSIM4GPbpPtr.real(BSIM4GPbpPtr.real() + (xcgbbi + gIgtotb)); // BSIM4GPbpPtr
        }

            
    // }

    // if (model.BSIM4rdsMod)
    // {   (*(instance.BSIM4DgpPtr) +=  gdtotg);
    //     (*(instance.BSIM4DspPtr) +=  gdtots);
    //     (*(instance.BSIM4DbpPtr) +=  gdtotb);
    //     (*(instance.BSIM4SdpPtr) +=  gstotd);
    //     (*(instance.BSIM4SgpPtr) +=  gstotg);
    //     (*(instance.BSIM4SbpPtr) +=  gstotb);
    // }

    
    if(!(dNodePrime < 0)){
        // BSIM4DPdpPtr
        arma::cx_double &BSIM4DPdpPtr = LHS(dNodePrime, dNodePrime);
        BSIM4DPdpPtr.imag(BSIM4DPdpPtr.imag() + xcddbr + gdsi + RevSumi); // BSIM4DPdpPtr +1
        BSIM4DPdpPtr.real(BSIM4DPdpPtr.real() + (gdpr + xcddbi + gdsr + instance.BSIM4gbd 
                            - gdtotd + RevSumr + gbdpdp - gIdtotd)); // BSIM4DPdpPtr
                
        if(!(dNode < 0)){
            // BSIM4DPdPtr
            arma::cx_double &BSIM4DPdPtr =  LHS(dNodePrime, dNode);
            BSIM4DPdPtr.real(BSIM4DPdPtr.real() - (gdpr + gdtot)); 
            // BSIM4DdpPtr
            arma::cx_double &BSIM4DdpPtr = LHS(dNode, dNodePrime);
            BSIM4DdpPtr.real(BSIM4DdpPtr.real() - (gdpr - gdtotd)); 
            // BSIM4DdPtr
            arma::cx_double &BSIM4DdPtr = LHS(dNode, dNode);
            BSIM4DdPtr.real(BSIM4DdPtr.real() + (gdpr + gdtot));
        }

        // BSIM4DPgpPtr
        if(!(gNodePrime < 0)){
            arma::cx_double &BSIM4DPgpPtr = LHS(dNodePrime, gNodePrime);
            BSIM4DPgpPtr.imag(BSIM4DPgpPtr.imag() + (xcdgbr + Gmi)); // BSIM4DPgpPtr +1
            BSIM4DPgpPtr.real(BSIM4DPgpPtr.real() + (Gmr + xcdgbi - gdtotg + gbdpg - gIdtotg)); // BSIM4DPgpPtr
        }

        if(!(sNodePrime < 0)){
            // BSIM4DPspPtr
            arma::cx_double &BSIM4DPspPtr = LHS(dNodePrime, sNodePrime);
            BSIM4DPspPtr.imag(BSIM4DPspPtr.imag() + (xcdsbr - gdsi - FwdSumi)); // BSIM4DPspPtr +1
            BSIM4DPspPtr.real(BSIM4DPspPtr.real() - (gdsr - xcdsbi + FwdSumr + gdtots - gbdpsp + gIdtots)); // BSIM4DPspPtr
            // BSIM4SPdpPtr
            arma::cx_double &BSIM4SPdpPtr = LHS(sNodePrime, dNodePrime);
            BSIM4SPdpPtr.imag(BSIM4SPdpPtr.imag() + (xcsdbr - gdsi - RevSumi)); // BSIM4SPdpPtr +1
            BSIM4SPdpPtr.real(BSIM4SPdpPtr.real() - (gdsr - xcsdbi + gstotd + RevSumr - gbspdp + gIstotd)); // BSIM4SPdpPtr
        }

        if(!(bNodePrime < 0)){
            // BSIM4DPbpPtr
            arma::cx_double &BSIM4DPbpPtr = LHS(dNodePrime, bNodePrime);
            BSIM4DPbpPtr.imag(BSIM4DPbpPtr.imag() + (xcdbbr + Gmbsi)); // BSIM4DPbpPtr +1
            BSIM4DPbpPtr.real(BSIM4DPbpPtr.real() - (gjbd + gdtotb - xcdbbi - Gmbsr - gbdpb + gIdtotb)); // BSIM4DPbpPtr
            // BSIM4BPdpPtr
            arma::cx_double &BSIM4BPdpPtr = LHS(bNodePrime, dNodePrime);
            BSIM4BPdpPtr.imag(BSIM4BPdpPtr.imag() + xcbdb); // BSIM4BPdpPtr +1
            BSIM4BPdpPtr.real(BSIM4BPdpPtr.real() - (gjbd - gbbdp + gIbtotd)); // BSIM4BPdpPtr
        }
    }

    if(!(sNodePrime < 0)){
        // BSIM4SPspPtr
        arma::cx_double &BSIM4SPspPtr = LHS(sNodePrime, sNodePrime);
        BSIM4SPspPtr.imag(BSIM4SPspPtr.imag() + xcssbr + gdsi + FwdSumi); // BSIM4SPspPtr +1
        BSIM4SPspPtr.real(BSIM4SPspPtr.real() + (gspr + xcssbi + gdsr + instance.BSIM4gbs
                            - gstots + FwdSumr + gbspsp - gIstots)); // BSIM4SPspPtr

        if(!(gNodePrime < 0)){
            // BSIM4SPgpPtr
            arma::cx_double &BSIM4SPgpPtr =  LHS(sNodePrime, gNodePrime);
            BSIM4SPgpPtr.imag(BSIM4SPgpPtr.imag() + (xcsgbr - Gmi)); // BSIM4SPgpPtr +1
            BSIM4SPgpPtr.real(BSIM4SPgpPtr.real() - (Gmr - xcsgbi + gstotg - gbspg + gIstotg)); // BSIM4SPgpPtr
        }

        if(!(sNode < 0)){
            // BSIM4SPsPtr
            arma::cx_double &BSIM4SPsPtr =  LHS(sNodePrime, sNode);
            BSIM4SPsPtr.real(BSIM4SPsPtr.real() - (gspr + gstot)); 
            // BSIM4SspPtr
            arma::cx_double &BSIM4SspPtr =  LHS(sNode, sNodePrime);
            BSIM4SspPtr.real(BSIM4SspPtr.real() - (gspr - gstots)); 
        }

        if(!(bNodePrime < 0)){
            // BSIM4SPbpPtr
            arma::cx_double &BSIM4SPbpPtr = LHS(sNodePrime, bNodePrime);
            BSIM4SPbpPtr.imag(BSIM4SPbpPtr.imag() + (xcsbbr - Gmbsi)); // BSIM4SPbpPtr +1
            BSIM4SPbpPtr.real(BSIM4SPbpPtr.real() - (gjbs + gstotb - xcsbbi + Gmbsr - gbspb + gIstotb)); // BSIM4SPbpPtr
            // BSIM4BPspPtr
            arma::cx_double &BSIM4BPspPtr = LHS(bNodePrime, sNodePrime);
            BSIM4BPspPtr.imag(BSIM4BPspPtr.imag() + xcbsb); // BSIM4BPspPtr +1
            BSIM4BPspPtr.real(BSIM4BPspPtr.real() - (gjbs - gbbsp + gIbtots)); // BSIM4BPspPtr
        }
    }

    if(!(sNode < 0)){
        // BSIM4SsPtr
        arma::cx_double &BSIM4SsPtr = LHS(sNode, sNode);
        BSIM4SsPtr.real(BSIM4SsPtr.real() + (gspr + gstot));
    }

    if(!(bNodePrime < 0)){
        // BSIM4BPbpPtr
        arma::cx_double &BSIM4BPbpPtr =  LHS(bNodePrime, bNodePrime);
        BSIM4BPbpPtr.imag(BSIM4BPbpPtr.imag() + xcbbb); // BSIM4BPbpPtr +1
        BSIM4BPbpPtr.real(BSIM4BPbpPtr.real() + (gjbd + gjbs - instance.BSIM4gbbs - gIbtotb)); // BSIM4BPbpPtr

        if(!(gNodePrime < 0)){
            //BSIM4BPgpPtr
            arma::cx_double &BSIM4BPgpPtr =  LHS(bNodePrime, gNodePrime);
            BSIM4BPgpPtr.imag(BSIM4BPgpPtr.imag() + xcbgb); // BSIM4BPgpPtr +1
            BSIM4BPgpPtr.real(BSIM4BPgpPtr.real() - (instance.BSIM4gbgs + gIbtotg)); // BSIM4BPgpPtr
        }
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
        LHS(dNodePrime, dNodePrime).real(LHS(dNodePrime, dNodePrime).real() + ggidld);
    }
    // BSIM4DPgpPtr
    if(!(dNodePrime < 0) && !(gNodePrime < 0)){
        LHS(dNodePrime, gNodePrime).real(LHS(dNodePrime, gNodePrime).real() + ggidlg);
    }
    // BSIM4DPspPtr
    if(!(dNodePrime < 0) && !(sNodePrime < 0)){
        LHS(dNodePrime, sNodePrime).real(LHS(dNodePrime, sNodePrime).real() - ((ggidlg + ggidld) + ggidlb));
    }
    // BSIM4DPbpPtr
    if(!(dNodePrime < 0) && !(bNodePrime < 0)){
        LHS(dNodePrime, bNodePrime).real(LHS(dNodePrime, bNodePrime).real() + ggidlb);
    }
    // BSIM4BPdpPtr
    if(!(bNodePrime < 0) && !(dNodePrime < 0)){
        LHS(bNodePrime, dNodePrime).real(LHS(bNodePrime, dNodePrime).real() - ggidld);
    }
    // BSIM4BPgpPtr
    if(!(bNodePrime < 0) && !(gNodePrime < 0)){
        LHS(bNodePrime, gNodePrime).real(LHS(bNodePrime, gNodePrime).real() - ggidlg);
    }
    // BSIM4BPspPtr
    if(!(bNodePrime < 0) && !(sNodePrime < 0)){
        LHS(bNodePrime, sNodePrime).real(LHS(bNodePrime, sNodePrime).real() + (ggidlg + ggidld) + ggidlb);
    }
    // BSIM4BPbpPtr
    if(!(bNodePrime < 0)){
        LHS(bNodePrime, bNodePrime).real(LHS(bNodePrime, bNodePrime).real() - ggidlb);
    }
    
    /* stamp gisl */
    // BSIM4SPdpPtr
    if(!(sNodePrime < 0) && !(dNodePrime < 0)){
        LHS(sNodePrime, dNodePrime).real(LHS(sNodePrime, dNodePrime).real() - ((ggisls + ggislg) + ggislb));
    }
    // BSIM4SPgpPtr
    if(!(sNodePrime < 0) && !(gNodePrime < 0)){
        LHS(sNodePrime, gNodePrime).real(LHS(sNodePrime, gNodePrime).real() + ggislg);
    }
    // BSIM4SPspPtr
    if(!(sNodePrime < 0)){
        LHS(sNodePrime, sNodePrime).real(LHS(sNodePrime, sNodePrime).real() + ggisls);
    }
    // BSIM4SPbpPtr
    if(!(sNodePrime < 0) && !(bNodePrime < 0)){
        LHS(sNodePrime, bNodePrime).real(LHS(sNodePrime, bNodePrime).real() + ggislb);
    }
    // BSIM4BPdpPtr
    if(!(bNodePrime < 0) && !(dNodePrime < 0)){
        LHS(bNodePrime, dNodePrime).real(LHS(bNodePrime, dNodePrime).real() + ((ggislg + ggisls) + ggislb));
    }
    // BSIM4BPgpPtr
    if(!(bNodePrime < 0) && !(gNodePrime < 0)){
        LHS(bNodePrime, gNodePrime).real(LHS(bNodePrime, gNodePrime).real() - ggislg);
    }
    // BSIM4BPspPtr
    if(!(bNodePrime < 0) && !(sNodePrime < 0)){
        LHS(bNodePrime, sNodePrime).real(LHS(bNodePrime, sNodePrime).real() - ggisls);
    }
    // BSIM4BPbpPtr
    if(!(bNodePrime < 0)){
        LHS(bNodePrime, bNodePrime).real(LHS(bNodePrime, bNodePrime).real() - ggislb);
    }

    // if (instance.BSIM4rbodyMod)
    // {   (*(instance.BSIM4DPdbPtr +1) +=  xcdbdb);
    //     (*(instance.BSIM4DPdbPtr) -=  instance.BSIM4gbd);
    //     (*(instance.BSIM4SPsbPtr +1) +=  xcsbsb);
    //     (*(instance.BSIM4SPsbPtr) -=  instance.BSIM4gbs);

    //     (*(instance.BSIM4DBdpPtr +1) +=  xcdbdb);
    //     (*(instance.BSIM4DBdpPtr) -=  instance.BSIM4gbd);
    //     (*(instance.BSIM4DBdbPtr +1) -=  xcdbdb);
    //     (*(instance.BSIM4DBdbPtr) +=  (instance.BSIM4gbd + instance.BSIM4grbpd 
    //                             + instance.BSIM4grbdb));
    //     (*(instance.BSIM4DBbpPtr) -=  instance.BSIM4grbpd);
    //     (*(instance.BSIM4DBbPtr) -=  instance.BSIM4grbdb);

    //     (*(instance.BSIM4BPdbPtr) -=  instance.BSIM4grbpd);
    //     (*(instance.BSIM4BPbPtr) -=  instance.BSIM4grbpb);
    //     (*(instance.BSIM4BPsbPtr) -=  instance.BSIM4grbps);
    //     (*(instance.BSIM4BPbpPtr) +=  (instance.BSIM4grbpd + instance.BSIM4grbps 
    //                             + instance.BSIM4grbpb));
    //     /* WDLiu: (-instance.BSIM4gbbs) already added to BPbpPtr */

    //     (*(instance.BSIM4SBspPtr +1) +=  xcsbsb);
    //     (*(instance.BSIM4SBspPtr) -=  instance.BSIM4gbs);
    //     (*(instance.BSIM4SBbpPtr) -=  instance.BSIM4grbps);
    //     (*(instance.BSIM4SBbPtr) -=  instance.BSIM4grbsb);
    //     (*(instance.BSIM4SBsbPtr +1) -=  xcsbsb);
    //     (*(instance.BSIM4SBsbPtr) +=  (instance.BSIM4gbs
    //                             + instance.BSIM4grbps + instance.BSIM4grbsb));

    //     (*(instance.BSIM4BdbPtr) -=  instance.BSIM4grbdb);
    //     (*(instance.BSIM4BbpPtr) -=  instance.BSIM4grbpb);
    //     (*(instance.BSIM4BsbPtr) -=  instance.BSIM4grbsb);
    //     (*(instance.BSIM4BbPtr) +=  (instance.BSIM4grbsb + instance.BSIM4grbdb
    //                         + instance.BSIM4grbpb));
    // }


    /*
    * WDLiu: The internal charge node generated for transient NQS is not needed for
    *        AC NQS. The following is not doing a real job, but we have to keep it;
    *        otherwise a singular AC NQS matrix may occur if the transient NQS is on.
    *        The charge node is isolated from the instance.
    */
    // if (instance.BSIM4trnqsMod)
    // {   (*(instance.BSIM4QqPtr) +=  1.0);
    //     (*(instance.BSIM4QgpPtr) += 0.0);
    //     (*(instance.BSIM4QdpPtr) += 0.0);
    //     (*(instance.BSIM4QspPtr) += 0.0);
    //     (*(instance.BSIM4QbpPtr) += 0.0);

    //     (*(instance.BSIM4DPqPtr) += 0.0);
    //     (*(instance.BSIM4SPqPtr) += 0.0);
    //     (*(instance.BSIM4GPqPtr) += 0.0);
    // }
    
    return 0;
}

} // namespace bsim4
#pragma once
#include "bsim4v82.hpp"
#include "VariantValue.hpp"

namespace bsim4{
VariantValue BSIM4ask(const BSIM4V82 &inst, int which)
{
// BSIM4instance *here = (BSIM4instance*)inst;

    switch(which) 
    {   case BSIM4_L:
            return VariantValue(inst.BSIM4l);
        case BSIM4_W:
            return VariantValue(inst.BSIM4w);
        // case BSIM4_M:
        //     return VariantValue(inst.BSIM4m);
        case BSIM4_NF:
            return VariantValue(inst.BSIM4nf);
        case BSIM4_MIN:
            return VariantValue(inst.BSIM4min);
        case BSIM4_AS:
            return VariantValue(inst.BSIM4sourceArea);
        case BSIM4_AD:
            return VariantValue(inst.BSIM4drainArea);
        case BSIM4_PS:
            return VariantValue(inst.BSIM4sourcePerimeter);
        case BSIM4_PD:
            return VariantValue(inst.BSIM4drainPerimeter);
        case BSIM4_NRS:
            return VariantValue(inst.BSIM4sourceSquares);
        case BSIM4_NRD:
            return VariantValue(inst.BSIM4drainSquares);
        case BSIM4_OFF:
            return VariantValue(inst.BSIM4off);
        case BSIM4_SA:
            return VariantValue(inst.BSIM4sa);
        case BSIM4_SB:
            return VariantValue(inst.BSIM4sb);
        case BSIM4_SD:
            return VariantValue(inst.BSIM4sd);
        case BSIM4_SCA:
            return VariantValue(inst.BSIM4sca);
        case BSIM4_SCB:
            return VariantValue(inst.BSIM4scb);
        case BSIM4_SCC:
            return VariantValue(inst.BSIM4scc);
        case BSIM4_SC:
            return VariantValue(inst.BSIM4sc);

        case BSIM4_RBSB:
            return VariantValue(inst.BSIM4rbsb);
        case BSIM4_RBDB:
            return VariantValue(inst.BSIM4rbdb);
        case BSIM4_RBPB:
            return VariantValue(inst.BSIM4rbpb);
        case BSIM4_RBPS:
            return VariantValue(inst.BSIM4rbps);
        case BSIM4_RBPD:
            return VariantValue(inst.BSIM4rbpd);
        case BSIM4_DELVTO:
            return VariantValue(inst.BSIM4delvto);
        // case BSIM4_MULU0:
        //     return VariantValue(inst.BSIM4mulu0);
        // case BSIM4_WNFLAG:
        //     return VariantValue(inst.BSIM4wnflag);
        case BSIM4_XGW:
            return VariantValue(inst.BSIM4xgw);
        case BSIM4_NGCON:
            return VariantValue(inst.BSIM4ngcon);
        case BSIM4_TRNQSMOD:
            return VariantValue(inst.BSIM4trnqsMod);
        case BSIM4_ACNQSMOD:
            return VariantValue(inst.BSIM4acnqsMod);
        case BSIM4_RBODYMOD:
            return VariantValue(inst.BSIM4rbodyMod);
        case BSIM4_RGATEMOD:
            return VariantValue(inst.BSIM4rgateMod);
        case BSIM4_GEOMOD:
            return VariantValue(inst.BSIM4geoMod);
        case BSIM4_RGEOMOD:
            return VariantValue(inst.BSIM4rgeoMod);
        case BSIM4_IC_VDS:
            return VariantValue(inst.BSIM4icVDS);
        case BSIM4_IC_VGS:
            return VariantValue(inst.BSIM4icVGS);
        case BSIM4_IC_VBS:
            return VariantValue(inst.BSIM4icVBS);
        case BSIM4_DNODE:
            return VariantValue(inst.BSIM4dNode);
        case BSIM4_GNODEEXT:
            return VariantValue(inst.BSIM4gNodeExt);
        case BSIM4_SNODE:
            return VariantValue(inst.BSIM4sNode);
        case BSIM4_BNODE:
            return VariantValue(inst.BSIM4bNode);
        case BSIM4_DNODEPRIME:
            return VariantValue(inst.BSIM4dNodePrime);
        case BSIM4_GNODEPRIME:
            return VariantValue(inst.BSIM4gNodePrime);
        case BSIM4_GNODEMID:
            return VariantValue(inst.BSIM4gNodeMid);
        case BSIM4_SNODEPRIME:
            return VariantValue(inst.BSIM4sNodePrime);
        case BSIM4_DBNODE:
            return VariantValue(inst.BSIM4dbNode);
        case BSIM4_BNODEPRIME:
            return VariantValue(inst.BSIM4bNodePrime);
        case BSIM4_SBNODE:
            return VariantValue(inst.BSIM4sbNode);
        case BSIM4_SOURCECONDUCT:
            // return VariantValue(inst.BSIM4sourceConductance * inst.BSIM4m);
            return VariantValue(inst.BSIM4sourceConductance);
        case BSIM4_DRAINCONDUCT:
            // return VariantValue(inst.BSIM4drainConductance * inst.BSIM4m);
            return VariantValue(inst.BSIM4drainConductance);
        case BSIM4_VBD:
            return VariantValue(inst.BSIM4states0[BSIM4vbd]);
        case BSIM4_VBS:
            return VariantValue(inst.BSIM4states0[BSIM4vbs]);
        case BSIM4_VGS:
            return VariantValue(inst.BSIM4states0[BSIM4vgs]);
        case BSIM4_VDS:
            return VariantValue(inst.BSIM4states0[BSIM4vds]);
        case BSIM4_CD:
            // return VariantValue(inst.BSIM4cd * inst.BSIM4m);
            return VariantValue(inst.BSIM4cd);
        case BSIM4_CBS:
            // return VariantValue(inst.BSIM4cbs * inst.BSIM4m);
            return VariantValue(inst.BSIM4cbs);
        case BSIM4_CBD:
            // return VariantValue(inst.BSIM4cbd * inst.BSIM4m);
            return VariantValue(inst.BSIM4cbd);
        case BSIM4_CSUB:
            // return VariantValue(inst.BSIM4csub * inst.BSIM4m);
            return VariantValue(inst.BSIM4csub);
        case BSIM4_QINV:
            // return VariantValue(inst.BSIM4qinv * inst.BSIM4m);
            return VariantValue(inst.BSIM4qinv);
        case BSIM4_IGIDL:
            // return VariantValue(inst.BSIM4Igidl * inst.BSIM4m);
            return VariantValue(inst.BSIM4Igidl);
        case BSIM4_IGISL:
            // return VariantValue(inst.BSIM4Igisl * inst.BSIM4m);
            return VariantValue(inst.BSIM4Igisl);
        case BSIM4_IGS:
            // return VariantValue(inst.BSIM4Igs * inst.BSIM4m);
            return VariantValue(inst.BSIM4Igs);
        case BSIM4_IGD:
            // return VariantValue(inst.BSIM4Igd * inst.BSIM4m);
            return VariantValue(inst.BSIM4Igd);
        case BSIM4_IGB:
            // return VariantValue(inst.BSIM4Igb * inst.BSIM4m);
            return VariantValue(inst.BSIM4Igb);
        case BSIM4_IGCS:
            // return VariantValue(inst.BSIM4Igcs * inst.BSIM4m);
            return VariantValue(inst.BSIM4Igcs);
        case BSIM4_IGCD:
            // return VariantValue(inst.BSIM4Igcd * inst.BSIM4m);
            return VariantValue(inst.BSIM4Igcd);
        case BSIM4_GM:
            // return VariantValue(inst.BSIM4gm * inst.BSIM4m);
            return VariantValue(inst.BSIM4gm);
        case BSIM4_GDS:
            // return VariantValue(inst.BSIM4gds * inst.BSIM4m);
            return VariantValue(inst.BSIM4gds);
        case BSIM4_GMBS:
            // return VariantValue(inst.BSIM4gmbs * inst.BSIM4m);
            return VariantValue(inst.BSIM4gmbs);
        case BSIM4_GBD:
            // return VariantValue(inst.BSIM4gbd * inst.BSIM4m);
            return VariantValue(inst.BSIM4gbd);
        case BSIM4_GBS:
            // return VariantValue(inst.BSIM4gbs * inst.BSIM4m);
            return VariantValue(inst.BSIM4gbs);
/*        case BSIM4_QB:
            value->rValue = *(ckt->CKTstate0 + inst.BSIM4qb); 
            return(OK); */
        case BSIM4_CQB:
            return VariantValue(inst.BSIM4states0[BSIM4cqb]);
/*        case BSIM4_QG:
            value->rValue = *(ckt->CKTstate0 + inst.BSIM4qg); 
            return(OK); */
        case BSIM4_CQG:
            return VariantValue(inst.BSIM4states0[BSIM4cqg]);
/*        case BSIM4_QD:
            value->rValue = *(ckt->CKTstate0 + inst.BSIM4qd); 
            return(OK); */
        case BSIM4_CQD:
            return VariantValue(inst.BSIM4states0[BSIM4cqd]);
/*        case BSIM4_QS:
            value->rValue = *(ckt->CKTstate0 + inst.BSIM4qs); 
            return(OK); */
        case BSIM4_QB:
            // return VariantValue(inst.BSIM4qbulk * inst.BSIM4m);
            return VariantValue(inst.BSIM4qbulk);
        case BSIM4_QG:
            // return VariantValue(inst.BSIM4qgate * inst.BSIM4m);
            return VariantValue(inst.BSIM4qgate);
        case BSIM4_QS:
            // return VariantValue(inst.BSIM4qsource * inst.BSIM4m);
            return VariantValue(inst.BSIM4qsrc);
        case BSIM4_QD:
            // return VariantValue(inst.BSIM4qdrn * inst.BSIM4m);
            return VariantValue(inst.BSIM4qdrn);
        case BSIM4_QDEF:
            return VariantValue(inst.BSIM4states0[BSIM4qdef]);
        case BSIM4_GCRG:
            // return VariantValue(inst.BSIM4gcrg * inst.BSIM4m);
            return VariantValue(inst.BSIM4gcrg);
        case BSIM4_GTAU:
            // return VariantValue(inst.BSIM4gtau * inst.BSIM4m);
            return VariantValue(inst.BSIM4gtau);
        case BSIM4_CGGB:
            // return VariantValue(inst.BSIM4cggb * inst.BSIM4m);
            return VariantValue(inst.BSIM4cggb);
        case BSIM4_CGDB:
            // return VariantValue(inst.BSIM4cgdb * inst.BSIM4m);
            return VariantValue(inst.BSIM4cgdb);
        case BSIM4_CGSB:
            // return VariantValue(inst.BSIM4cgsb * inst.BSIM4m);
            return VariantValue(inst.BSIM4cgsb);
        case BSIM4_CDGB:
            // return VariantValue(inst.BSIM4cdgb * inst.BSIM4m);
            return VariantValue(inst.BSIM4cdgb);
        case BSIM4_CDDB:
            // return VariantValue(inst.BSIM4cddb * inst.BSIM4m);
            return VariantValue(inst.BSIM4cddb);
        case BSIM4_CDSB:
            // return VariantValue(inst.BSIM4cdsb * inst.BSIM4m);
            return VariantValue(inst.BSIM4cdsb);
        case BSIM4_CBGB:
            // return VariantValue(inst.BSIM4cbgb * inst.BSIM4m);
            return VariantValue(inst.BSIM4cbgb);
        case BSIM4_CBDB:
            // return VariantValue(inst.BSIM4cbdb * inst.BSIM4m);
            return VariantValue(inst.BSIM4cbdb);
        case BSIM4_CBSB:
            // return VariantValue(inst.BSIM4cbsb * inst.BSIM4m);
            return VariantValue(inst.BSIM4cbsb);
        case BSIM4_CSGB:
            // return VariantValue(inst.BSIM4csgb * inst.BSIM4m);
            return VariantValue(inst.BSIM4csgb);
        case BSIM4_CSDB:
            // return VariantValue(inst.BSIM4csdb * inst.BSIM4m);
            return VariantValue(inst.BSIM4csdb);
        case BSIM4_CSSB:
            // return VariantValue(inst.BSIM4cssb * inst.BSIM4m);
            return VariantValue(inst.BSIM4cssb);
        case BSIM4_CGBB:
            // return VariantValue(inst.BSIM4cgbb * inst.BSIM4m);
            return VariantValue(inst.BSIM4cgbb);
        case BSIM4_CDBB:
            // return VariantValue(inst.BSIM4cdbb * inst.BSIM4m);
            return VariantValue(inst.BSIM4cdbb);
        case BSIM4_CSBB:
            // return VariantValue(inst.BSIM4csbb * inst.BSIM4m);
            return VariantValue(inst.BSIM4csbb);
        case BSIM4_CBBB:
            // return VariantValue(inst.BSIM4cbbb * inst.BSIM4m);
            return VariantValue(inst.BSIM4cbbb);
        case BSIM4_CAPBD:
            // return VariantValue(inst.BSIM4capbd * inst.BSIM4m);
            return VariantValue(inst.BSIM4capbd);
        case BSIM4_CAPBS:
            // return VariantValue(inst.BSIM4capbs * inst.BSIM4m);
            return VariantValue(inst.BSIM4capbs);
        case BSIM4_VON:
            return VariantValue(inst.BSIM4von);
        case BSIM4_VDSAT:
            return VariantValue(inst.BSIM4vdsat);
        case BSIM4_QBS:
            return VariantValue(inst.BSIM4states0[BSIM4qbs]);
        case BSIM4_QBD:
            return VariantValue(inst.BSIM4states0[BSIM4qbd]);
        // case BSIM4_VGSTEFF:
        //     return VariantValue(inst.BSIM4Vgsteff);
        // case BSIM4_VDSEFF:
        //     return VariantValue(inst.BSIM4Vdseff);
        // case BSIM4_CGSO:
        //     // return VariantValue(inst.BSIM4cgso * inst.BSIM4m);
        //     return VariantValue(inst.BSIM4cgso);
        // case BSIM4_CGDO:
        //     // return VariantValue(inst.BSIM4cgdo * inst.BSIM4m);
        //     return VariantValue(inst.BSIM4cgdo);
        // case BSIM4_CGBO:
        //     // return VariantValue(inst.pParam->BSIM4cgbo * inst.BSIM4m);
        //     return VariantValue(inst.BSIM4cgbo);
        // case BSIM4_WEFF:
        //     value->rValue = inst.pParam->BSIM4weff;
        //     return(OK);
        // case BSIM4_LEFF:
        //     value->rValue = inst.pParam->BSIM4leff;
        //     return(OK);
        default:
            std::cerr << "BSIM4ask: Invalid parameter index: " << which << std::endl;
            return VariantValue(); // Return an empty VariantValue or handle error appropriately
    }
    /* NOTREACHED */
}
} // namespace bsim4

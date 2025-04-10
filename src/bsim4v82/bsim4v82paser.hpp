/*
   b4mpar.c - BSIM4v4.8.2
   b4par.c  - BSIM4v4.8.2
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
#include <stdio.h>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <algorithm>
#include <cctype> // std::tolower
#include "bsim4v82.hpp"
#include "bsim4v82map.hpp"
#include "VariantValue.hpp"


namespace bsim4{

inline constexpr double CONSTCtoK = 273.15;

// 1. set model parameters in BSIM4model structure
int BSIM4mParam(int param, const VariantValue &value, BSIM4model &mod){

    switch(param)
    {   case  BSIM4_MOD_MOBMOD :
            mod.BSIM4mobMod = value.get_as<int>().value();
            mod.BSIM4mobModGiven = true;
            break;
        case  BSIM4_MOD_BINUNIT :
            mod.BSIM4binUnit = value.get_as<int>().value();
            mod.BSIM4binUnitGiven = true;
            break;
        case  BSIM4_MOD_PARAMCHK :
            mod.BSIM4paramChk = value.get_as<int>().value();
            mod.BSIM4paramChkGiven = true;
            break;
        case  BSIM4_MOD_CVCHARGEMOD :
            mod.BSIM4cvchargeMod = value.get_as<int>().value();
            mod.BSIM4cvchargeModGiven = true;
            break;
        case  BSIM4_MOD_CAPMOD :
            mod.BSIM4capMod = value.get_as<int>().value();
            mod.BSIM4capModGiven = true;
            break;
        case  BSIM4_MOD_DIOMOD :
            mod.BSIM4dioMod = value.get_as<int>().value();
            mod.BSIM4dioModGiven = true;
            break;
        case  BSIM4_MOD_RDSMOD :
            mod.BSIM4rdsMod = value.get_as<int>().value();
            mod.BSIM4rdsModGiven = true;
            break;
        case  BSIM4_MOD_TRNQSMOD :
            mod.BSIM4trnqsMod = value.get_as<int>().value();
            mod.BSIM4trnqsModGiven = true;
            break;
        case  BSIM4_MOD_ACNQSMOD :
            mod.BSIM4acnqsMod = value.get_as<int>().value();
            mod.BSIM4acnqsModGiven = true;
            break;
        case  BSIM4_MOD_RBODYMOD :
            mod.BSIM4rbodyMod = value.get_as<int>().value();
            mod.BSIM4rbodyModGiven = true;
            break;
        case  BSIM4_MOD_RGATEMOD :
            mod.BSIM4rgateMod = value.get_as<int>().value();
            mod.BSIM4rgateModGiven = true;
            break;
        case  BSIM4_MOD_PERMOD :
            mod.BSIM4perMod = value.get_as<int>().value();
            mod.BSIM4perModGiven = true;
            break;
        case  BSIM4_MOD_GEOMOD :
            mod.BSIM4geoMod = value.get_as<int>().value();
            mod.BSIM4geoModGiven = true;
            break;
        case  BSIM4_MOD_FNOIMOD :
            mod.BSIM4fnoiMod = value.get_as<int>().value();
            mod.BSIM4fnoiModGiven = true;
            break;
        case  BSIM4_MOD_TNOIMOD :
            mod.BSIM4tnoiMod = value.get_as<int>().value();
            mod.BSIM4tnoiModGiven = true;
            break;
        case  BSIM4_MOD_MTRLMOD :
            mod.BSIM4mtrlMod = value.get_as<int>().value();
            mod.BSIM4mtrlModGiven = true;
            break;
        case  BSIM4_MOD_MTRLCOMPATMOD :
            mod.BSIM4mtrlCompatMod = value.get_as<int>().value();
            mod.BSIM4mtrlCompatModGiven = true;
            break;
    case  BSIM4_MOD_GIDLMOD :   /* v4.7 New GIDL/GISL */
            mod.BSIM4gidlMod = value.get_as<int>().value();
            mod.BSIM4gidlModGiven = true;
            break;
        case  BSIM4_MOD_IGCMOD :
            mod.BSIM4igcMod = value.get_as<int>().value();
            mod.BSIM4igcModGiven = true;
            break;
        case  BSIM4_MOD_IGBMOD :
            mod.BSIM4igbMod = value.get_as<int>().value();
            mod.BSIM4igbModGiven = true;
            break;
        case  BSIM4_MOD_TEMPMOD :
            mod.BSIM4tempMod = value.get_as<int>().value();
            mod.BSIM4tempModGiven = true;
            break;

        case  BSIM4_MOD_VERSION :
            mod.BSIM4version = value.get<std::string>();
            mod.BSIM4versionGiven = true;
            break;
        case  BSIM4_MOD_TOXREF :
            mod.BSIM4toxref = value.get_as<double>().value();
            mod.BSIM4toxrefGiven = true;
            break;
        case  BSIM4_MOD_EOT :
            mod.BSIM4eot = value.get_as<double>().value();
            mod.BSIM4eotGiven = true;
            break;
        case  BSIM4_MOD_VDDEOT :
            mod.BSIM4vddeot = value.get_as<double>().value();
            mod.BSIM4vddeotGiven = true;
            break;
        case  BSIM4_MOD_TEMPEOT :
            mod.BSIM4tempeot = value.get_as<double>().value();
            mod.BSIM4tempeotGiven = true;
            break;
    case  BSIM4_MOD_LEFFEOT :
            mod.BSIM4leffeot = value.get_as<double>().value();
            mod.BSIM4leffeotGiven = true;
            break;
    case  BSIM4_MOD_WEFFEOT :
            mod.BSIM4weffeot = value.get_as<double>().value();
            mod.BSIM4weffeotGiven = true;
            break;
        case  BSIM4_MOD_ADOS :
            mod.BSIM4ados = value.get_as<double>().value();
            mod.BSIM4adosGiven = true;
            break;
        case  BSIM4_MOD_BDOS :
            mod.BSIM4bdos = value.get_as<double>().value();
            mod.BSIM4bdosGiven = true;
            break;
        case  BSIM4_MOD_TOXE :
            mod.BSIM4toxe = value.get_as<double>().value();
            mod.BSIM4toxeGiven = true;
            break;
        case  BSIM4_MOD_TOXP :
            mod.BSIM4toxp = value.get_as<double>().value();
            mod.BSIM4toxpGiven = true;
            break;
        case  BSIM4_MOD_TOXM :
            mod.BSIM4toxm = value.get_as<double>().value();
            mod.BSIM4toxmGiven = true;
            break;
        case  BSIM4_MOD_DTOX :
            mod.BSIM4dtox = value.get_as<double>().value();
            mod.BSIM4dtoxGiven = true;
            break;
        case  BSIM4_MOD_EPSROX :
            mod.BSIM4epsrox = value.get_as<double>().value();
            mod.BSIM4epsroxGiven = true;
            break;
        case  BSIM4_MOD_CDSC :
            mod.BSIM4cdsc = value.get_as<double>().value();
            mod.BSIM4cdscGiven = true;
            break;
        case  BSIM4_MOD_CDSCB :
            mod.BSIM4cdscb = value.get_as<double>().value();
            mod.BSIM4cdscbGiven = true;
            break;
        case  BSIM4_MOD_CDSCD :
            mod.BSIM4cdscd = value.get_as<double>().value();
            mod.BSIM4cdscdGiven = true;
            break;
        case  BSIM4_MOD_CIT :
            mod.BSIM4cit = value.get_as<double>().value();
            mod.BSIM4citGiven = true;
            break;
        case  BSIM4_MOD_NFACTOR :
            mod.BSIM4nfactor = value.get_as<double>().value();
            mod.BSIM4nfactorGiven = true;
            break;
        case BSIM4_MOD_XJ:
            mod.BSIM4xj = value.get_as<double>().value();
            mod.BSIM4xjGiven = true;
            break;
        case BSIM4_MOD_VSAT:
            mod.BSIM4vsat = value.get_as<double>().value();
            mod.BSIM4vsatGiven = true;
            break;
        case BSIM4_MOD_A0:
            mod.BSIM4a0 = value.get_as<double>().value();
            mod.BSIM4a0Given = true;
            break;

        case BSIM4_MOD_AGS:
            mod.BSIM4ags= value.get_as<double>().value();
            mod.BSIM4agsGiven = true;
            break;

        case BSIM4_MOD_A1:
            mod.BSIM4a1 = value.get_as<double>().value();
            mod.BSIM4a1Given = true;
            break;
        case BSIM4_MOD_A2:
            mod.BSIM4a2 = value.get_as<double>().value();
            mod.BSIM4a2Given = true;
            break;
        case BSIM4_MOD_AT:
            mod.BSIM4at = value.get_as<double>().value();
            mod.BSIM4atGiven = true;
            break;
        case BSIM4_MOD_KETA:
            mod.BSIM4keta = value.get_as<double>().value();
            mod.BSIM4ketaGiven = true;
            break;
        case BSIM4_MOD_NSUB:
            mod.BSIM4nsub = value.get_as<double>().value();
            mod.BSIM4nsubGiven = true;
            break;
        case BSIM4_MOD_PHIG:
        mod.BSIM4phig = value.get_as<double>().value();
        mod.BSIM4phigGiven = true;
        break;
        case BSIM4_MOD_EPSRGATE:
        mod.BSIM4epsrgate = value.get_as<double>().value();
        mod.BSIM4epsrgateGiven = true;
        break;
        case BSIM4_MOD_EASUB:
            mod.BSIM4easub = value.get_as<double>().value();
            mod.BSIM4easubGiven = true;
            break;
        case BSIM4_MOD_EPSRSUB:
            mod.BSIM4epsrsub = value.get_as<double>().value();
            mod.BSIM4epsrsubGiven = true;
            break;
        case BSIM4_MOD_NI0SUB:
            mod.BSIM4ni0sub = value.get_as<double>().value();
            mod.BSIM4ni0subGiven = true;
            break;
        case BSIM4_MOD_BG0SUB:
            mod.BSIM4bg0sub = value.get_as<double>().value();
            mod.BSIM4bg0subGiven = true;
            break;
        case BSIM4_MOD_TBGASUB:
            mod.BSIM4tbgasub = value.get_as<double>().value();
            mod.BSIM4tbgasubGiven = true;
            break;
        case BSIM4_MOD_TBGBSUB:
            mod.BSIM4tbgbsub = value.get_as<double>().value();
            mod.BSIM4tbgbsubGiven = true;
            break;
        case BSIM4_MOD_NDEP:
            mod.BSIM4ndep = value.get_as<double>().value();
            mod.BSIM4ndepGiven = true;
        if (mod.BSIM4ndep > 1.0e20)
        mod.BSIM4ndep *= 1.0e-6;
            break;
        case BSIM4_MOD_NSD:
            mod.BSIM4nsd = value.get_as<double>().value();
            mod.BSIM4nsdGiven = true;
            if (mod.BSIM4nsd > 1.0e23)
                mod.BSIM4nsd *= 1.0e-6;
            break;
        case BSIM4_MOD_NGATE:
            mod.BSIM4ngate = value.get_as<double>().value();
            mod.BSIM4ngateGiven = true;
        if (mod.BSIM4ngate > 1.0e23)
        mod.BSIM4ngate *= 1.0e-6;
            break;
        case BSIM4_MOD_GAMMA1:
            mod.BSIM4gamma1 = value.get_as<double>().value();
            mod.BSIM4gamma1Given = true;
            break;
        case BSIM4_MOD_GAMMA2:
            mod.BSIM4gamma2 = value.get_as<double>().value();
            mod.BSIM4gamma2Given = true;
            break;
        case BSIM4_MOD_VBX:
            mod.BSIM4vbx = value.get_as<double>().value();
            mod.BSIM4vbxGiven = true;
            break;
        case BSIM4_MOD_VBM:
            mod.BSIM4vbm = value.get_as<double>().value();
            mod.BSIM4vbmGiven = true;
            break;
        case BSIM4_MOD_XT:
            mod.BSIM4xt = value.get_as<double>().value();
            mod.BSIM4xtGiven = true;
            break;
        case  BSIM4_MOD_K1:
            mod.BSIM4k1 = value.get_as<double>().value();
            mod.BSIM4k1Given = true;
            break;
        case  BSIM4_MOD_KT1:
            mod.BSIM4kt1 = value.get_as<double>().value();
            mod.BSIM4kt1Given = true;
            break;
        case  BSIM4_MOD_KT1L:
            mod.BSIM4kt1l = value.get_as<double>().value();
            mod.BSIM4kt1lGiven = true;
            break;
        case  BSIM4_MOD_KT2:
            mod.BSIM4kt2 = value.get_as<double>().value();
            mod.BSIM4kt2Given = true;
            break;
        case  BSIM4_MOD_K2:
            mod.BSIM4k2 = value.get_as<double>().value();
            mod.BSIM4k2Given = true;
            break;
        case  BSIM4_MOD_K3:
            mod.BSIM4k3 = value.get_as<double>().value();
            mod.BSIM4k3Given = true;
            break;
        case  BSIM4_MOD_K3B:
            mod.BSIM4k3b = value.get_as<double>().value();
            mod.BSIM4k3bGiven = true;
            break;
        case  BSIM4_MOD_LPE0:
            mod.BSIM4lpe0 = value.get_as<double>().value();
            mod.BSIM4lpe0Given = true;
            break;
        case  BSIM4_MOD_LPEB:
            mod.BSIM4lpeb = value.get_as<double>().value();
            mod.BSIM4lpebGiven = true;
            break;
        case  BSIM4_MOD_DVTP0:
            mod.BSIM4dvtp0 = value.get_as<double>().value();
            mod.BSIM4dvtp0Given = true;
            break;
        case  BSIM4_MOD_DVTP1:
            mod.BSIM4dvtp1 = value.get_as<double>().value();
            mod.BSIM4dvtp1Given = true;
            break;
        case  BSIM4_MOD_DVTP2:     /* New DIBL/Rout */
            mod.BSIM4dvtp2 = value.get_as<double>().value();
            mod.BSIM4dvtp2Given = true;
            break;
    case  BSIM4_MOD_DVTP3:
            mod.BSIM4dvtp3 = value.get_as<double>().value();
            mod.BSIM4dvtp3Given = true;
            break;
    case  BSIM4_MOD_DVTP4:
            mod.BSIM4dvtp4 = value.get_as<double>().value();
            mod.BSIM4dvtp4Given = true;
            break;
    case  BSIM4_MOD_DVTP5:
            mod.BSIM4dvtp5 = value.get_as<double>().value();
            mod.BSIM4dvtp5Given = true;
            break;
        case  BSIM4_MOD_W0:
            mod.BSIM4w0 = value.get_as<double>().value();
            mod.BSIM4w0Given = true;
            break;
        case  BSIM4_MOD_DVT0:
            mod.BSIM4dvt0 = value.get_as<double>().value();
            mod.BSIM4dvt0Given = true;
            break;
        case  BSIM4_MOD_DVT1:
            mod.BSIM4dvt1 = value.get_as<double>().value();
            mod.BSIM4dvt1Given = true;
            break;
        case  BSIM4_MOD_DVT2:
            mod.BSIM4dvt2 = value.get_as<double>().value();
            mod.BSIM4dvt2Given = true;
            break;
        case  BSIM4_MOD_DVT0W:
            mod.BSIM4dvt0w = value.get_as<double>().value();
            mod.BSIM4dvt0wGiven = true;
            break;
        case  BSIM4_MOD_DVT1W:
            mod.BSIM4dvt1w = value.get_as<double>().value();
            mod.BSIM4dvt1wGiven = true;
            break;
        case  BSIM4_MOD_DVT2W:
            mod.BSIM4dvt2w = value.get_as<double>().value();
            mod.BSIM4dvt2wGiven = true;
            break;
        case  BSIM4_MOD_DROUT:
            mod.BSIM4drout = value.get_as<double>().value();
            mod.BSIM4droutGiven = true;
            break;
        case  BSIM4_MOD_DSUB:
            mod.BSIM4dsub = value.get_as<double>().value();
            mod.BSIM4dsubGiven = true;
            break;
        case BSIM4_MOD_VTH0:
            mod.BSIM4vth0 = value.get_as<double>().value();
            mod.BSIM4vth0Given = true;
            break;
        case BSIM4_MOD_EU:
            mod.BSIM4eu = value.get_as<double>().value();
            mod.BSIM4euGiven = true;
            break;
        case BSIM4_MOD_UCS:
            mod.BSIM4ucs = value.get_as<double>().value();
            mod.BSIM4ucsGiven = true;
            break;
        case BSIM4_MOD_UA:
            mod.BSIM4ua = value.get_as<double>().value();
            mod.BSIM4uaGiven = true;
            break;
        case BSIM4_MOD_UA1:
            mod.BSIM4ua1 = value.get_as<double>().value();
            mod.BSIM4ua1Given = true;
            break;
        case BSIM4_MOD_UB:
            mod.BSIM4ub = value.get_as<double>().value();
            mod.BSIM4ubGiven = true;
            break;
        case BSIM4_MOD_UB1:
            mod.BSIM4ub1 = value.get_as<double>().value();
            mod.BSIM4ub1Given = true;
            break;
        case BSIM4_MOD_UC:
            mod.BSIM4uc = value.get_as<double>().value();
            mod.BSIM4ucGiven = true;
            break;
        case BSIM4_MOD_UC1:
            mod.BSIM4uc1 = value.get_as<double>().value();
            mod.BSIM4uc1Given = true;
            break;
        case  BSIM4_MOD_U0 :
            mod.BSIM4u0 = value.get_as<double>().value();
            mod.BSIM4u0Given = true;
            break;
        case  BSIM4_MOD_UTE :
            mod.BSIM4ute = value.get_as<double>().value();
            mod.BSIM4uteGiven = true;
            break;
        case  BSIM4_MOD_UCSTE :
            mod.BSIM4ucste = value.get_as<double>().value();
            mod.BSIM4ucsteGiven = true;
            break;
        case BSIM4_MOD_UD:
            mod.BSIM4ud = value.get_as<double>().value();
            mod.BSIM4udGiven = true;
            break;
        case BSIM4_MOD_UD1:
            mod.BSIM4ud1 = value.get_as<double>().value();
            mod.BSIM4ud1Given = true;
            break;
        case BSIM4_MOD_UP:
            mod.BSIM4up = value.get_as<double>().value();
            mod.BSIM4upGiven = true;
            break;
        case BSIM4_MOD_LP:
            mod.BSIM4lp = value.get_as<double>().value();
            mod.BSIM4lpGiven = true;
            break;
        case BSIM4_MOD_LUD:
            mod.BSIM4lud = value.get_as<double>().value();
            mod.BSIM4ludGiven = true;
            break;
        case BSIM4_MOD_LUD1:
            mod.BSIM4lud1 = value.get_as<double>().value();
            mod.BSIM4lud1Given = true;
            break;
        case BSIM4_MOD_LUP:
            mod.BSIM4lup = value.get_as<double>().value();
            mod.BSIM4lupGiven = true;
            break;
        case BSIM4_MOD_LLP:
            mod.BSIM4llp = value.get_as<double>().value();
            mod.BSIM4llpGiven = true;
            break;
        case BSIM4_MOD_WUD:
            mod.BSIM4wud = value.get_as<double>().value();
            mod.BSIM4wudGiven = true;
            break;
        case BSIM4_MOD_WUD1:
            mod.BSIM4wud1 = value.get_as<double>().value();
            mod.BSIM4wud1Given = true;
            break;
        case BSIM4_MOD_WUP:
            mod.BSIM4wup = value.get_as<double>().value();
            mod.BSIM4wupGiven = true;
            break;
        case BSIM4_MOD_WLP:
            mod.BSIM4wlp = value.get_as<double>().value();
            mod.BSIM4wlpGiven = true;
            break;
        case BSIM4_MOD_PUD:
            mod.BSIM4pud = value.get_as<double>().value();
            mod.BSIM4pudGiven = true;
            break;
        case BSIM4_MOD_PUD1:
            mod.BSIM4pud1 = value.get_as<double>().value();
            mod.BSIM4pud1Given = true;
            break;
        case BSIM4_MOD_PUP:
            mod.BSIM4pup = value.get_as<double>().value();
            mod.BSIM4pupGiven = true;
            break;
        case BSIM4_MOD_PLP:
            mod.BSIM4plp = value.get_as<double>().value();
            mod.BSIM4plpGiven = true;
            break;
        case BSIM4_MOD_VOFF:
            mod.BSIM4voff = value.get_as<double>().value();
            mod.BSIM4voffGiven = true;
            break;
        case BSIM4_MOD_TVOFF:
            mod.BSIM4tvoff = value.get_as<double>().value();
            mod.BSIM4tvoffGiven = true;
            break;
        case BSIM4_MOD_TNFACTOR:    /* v4.7 temp dep of leakage current  */
            mod.BSIM4tnfactor = value.get_as<double>().value();
            mod.BSIM4tnfactorGiven = true;
            break;
        case BSIM4_MOD_TETA0:       /* v4.7 temp dep of leakage current  */
            mod.BSIM4teta0 = value.get_as<double>().value();
            mod.BSIM4teta0Given = true;
            break;
        case BSIM4_MOD_TVOFFCV:     /* v4.7 temp dep of leakage current  */
            mod.BSIM4tvoffcv = value.get_as<double>().value();
            mod.BSIM4tvoffcvGiven = true;
            break;
        case BSIM4_MOD_VOFFL:
            mod.BSIM4voffl = value.get_as<double>().value();
            mod.BSIM4vofflGiven = true;
            break;
        case BSIM4_MOD_VOFFCVL:
            mod.BSIM4voffcvl = value.get_as<double>().value();
            mod.BSIM4voffcvlGiven = true;
            break;
        case BSIM4_MOD_MINV:
            mod.BSIM4minv = value.get_as<double>().value();
            mod.BSIM4minvGiven = true;
            break;
        case BSIM4_MOD_MINVCV:
            mod.BSIM4minvcv = value.get_as<double>().value();
            mod.BSIM4minvcvGiven = true;
            break;
        case BSIM4_MOD_FPROUT:
            mod.BSIM4fprout = value.get_as<double>().value();
            mod.BSIM4fproutGiven = true;
            break;
        case BSIM4_MOD_PDITS:
            mod.BSIM4pdits = value.get_as<double>().value();
            mod.BSIM4pditsGiven = true;
            break;
        case BSIM4_MOD_PDITSD:
            mod.BSIM4pditsd = value.get_as<double>().value();
            mod.BSIM4pditsdGiven = true;
            break;
        case BSIM4_MOD_PDITSL:
            mod.BSIM4pditsl = value.get_as<double>().value();
            mod.BSIM4pditslGiven = true;
            break;
        case  BSIM4_MOD_DELTA :
            mod.BSIM4delta = value.get_as<double>().value();
            mod.BSIM4deltaGiven = true;
            break;
        case BSIM4_MOD_RDSW:
            mod.BSIM4rdsw = value.get_as<double>().value();
            mod.BSIM4rdswGiven = true;
            break;
        case BSIM4_MOD_RDSWMIN:
            mod.BSIM4rdswmin = value.get_as<double>().value();
            mod.BSIM4rdswminGiven = true;
            break;
        case BSIM4_MOD_RDWMIN:
            mod.BSIM4rdwmin = value.get_as<double>().value();
            mod.BSIM4rdwminGiven = true;
            break;
        case BSIM4_MOD_RSWMIN:
            mod.BSIM4rswmin = value.get_as<double>().value();
            mod.BSIM4rswminGiven = true;
            break;
        case BSIM4_MOD_RDW:
            mod.BSIM4rdw = value.get_as<double>().value();
            mod.BSIM4rdwGiven = true;
            break;
        case BSIM4_MOD_RSW:
            mod.BSIM4rsw = value.get_as<double>().value();
            mod.BSIM4rswGiven = true;
            break;
        case BSIM4_MOD_PRWG:
            mod.BSIM4prwg = value.get_as<double>().value();
            mod.BSIM4prwgGiven = true;
            break;
        case BSIM4_MOD_PRWB:
            mod.BSIM4prwb = value.get_as<double>().value();
            mod.BSIM4prwbGiven = true;
            break;
        case BSIM4_MOD_PRT:
            mod.BSIM4prt = value.get_as<double>().value();
            mod.BSIM4prtGiven = true;
            break;
        case BSIM4_MOD_ETA0:
            mod.BSIM4eta0 = value.get_as<double>().value();
            mod.BSIM4eta0Given = true;
            break;
        case BSIM4_MOD_ETAB:
            mod.BSIM4etab = value.get_as<double>().value();
            mod.BSIM4etabGiven = true;
            break;
        case BSIM4_MOD_PCLM:
            mod.BSIM4pclm = value.get_as<double>().value();
            mod.BSIM4pclmGiven = true;
            break;
        case BSIM4_MOD_PDIBL1:
            mod.BSIM4pdibl1 = value.get_as<double>().value();
            mod.BSIM4pdibl1Given = true;
            break;
        case BSIM4_MOD_PDIBL2:
            mod.BSIM4pdibl2 = value.get_as<double>().value();
            mod.BSIM4pdibl2Given = true;
            break;
        case BSIM4_MOD_PDIBLB:
            mod.BSIM4pdiblb = value.get_as<double>().value();
            mod.BSIM4pdiblbGiven = true;
            break;
        case BSIM4_MOD_PSCBE1:
            mod.BSIM4pscbe1 = value.get_as<double>().value();
            mod.BSIM4pscbe1Given = true;
            break;
        case BSIM4_MOD_PSCBE2:
            mod.BSIM4pscbe2 = value.get_as<double>().value();
            mod.BSIM4pscbe2Given = true;
            break;
        case BSIM4_MOD_PVAG:
            mod.BSIM4pvag = value.get_as<double>().value();
            mod.BSIM4pvagGiven = true;
            break;
        case  BSIM4_MOD_WR :
            mod.BSIM4wr = value.get_as<double>().value();
            mod.BSIM4wrGiven = true;
            break;
        case  BSIM4_MOD_DWG :
            mod.BSIM4dwg = value.get_as<double>().value();
            mod.BSIM4dwgGiven = true;
            break;
        case  BSIM4_MOD_DWB :
            mod.BSIM4dwb = value.get_as<double>().value();
            mod.BSIM4dwbGiven = true;
            break;
        case  BSIM4_MOD_B0 :
            mod.BSIM4b0 = value.get_as<double>().value();
            mod.BSIM4b0Given = true;
            break;
        case  BSIM4_MOD_B1 :
            mod.BSIM4b1 = value.get_as<double>().value();
            mod.BSIM4b1Given = true;
            break;
        case  BSIM4_MOD_ALPHA0 :
            mod.BSIM4alpha0 = value.get_as<double>().value();
            mod.BSIM4alpha0Given = true;
            break;
        case  BSIM4_MOD_ALPHA1 :
            mod.BSIM4alpha1 = value.get_as<double>().value();
            mod.BSIM4alpha1Given = true;
            break;
        case  BSIM4_MOD_PHIN :
            mod.BSIM4phin = value.get_as<double>().value();
            mod.BSIM4phinGiven = true;
            break;
        case  BSIM4_MOD_AGIDL :
            mod.BSIM4agidl = value.get_as<double>().value();
            mod.BSIM4agidlGiven = true;
            break;
        case  BSIM4_MOD_BGIDL :
            mod.BSIM4bgidl = value.get_as<double>().value();
            mod.BSIM4bgidlGiven = true;
            break;
        case  BSIM4_MOD_CGIDL :
            mod.BSIM4cgidl = value.get_as<double>().value();
            mod.BSIM4cgidlGiven = true;
            break;
        case  BSIM4_MOD_EGIDL :
            mod.BSIM4egidl = value.get_as<double>().value();
            mod.BSIM4egidlGiven = true;
            break;
    case  BSIM4_MOD_FGIDL :         /* v4.7 New GIDL/GISL */
            mod.BSIM4fgidl = value.get_as<double>().value();
            mod.BSIM4fgidlGiven = true;
            break;
    case  BSIM4_MOD_KGIDL :         /* v4.7 New GIDL/GISL */
            mod.BSIM4kgidl = value.get_as<double>().value();
            mod.BSIM4kgidlGiven = true;
            break;
    case  BSIM4_MOD_RGIDL :         /* v4.7 New GIDL/GISL */
            mod.BSIM4rgidl = value.get_as<double>().value();
            mod.BSIM4rgidlGiven = true;
            break;
        case  BSIM4_MOD_AGISL :
            mod.BSIM4agisl = value.get_as<double>().value();
            mod.BSIM4agislGiven = true;
            break;
        case  BSIM4_MOD_BGISL :
            mod.BSIM4bgisl = value.get_as<double>().value();
            mod.BSIM4bgislGiven = true;
            break;
        case  BSIM4_MOD_CGISL :
            mod.BSIM4cgisl = value.get_as<double>().value();
            mod.BSIM4cgislGiven = true;
            break;
        case  BSIM4_MOD_EGISL :
            mod.BSIM4egisl = value.get_as<double>().value();
            mod.BSIM4egislGiven = true;
            break;
    case  BSIM4_MOD_FGISL :         /* v4.7 New GIDL/GISL */
            mod.BSIM4fgisl = value.get_as<double>().value();
            mod.BSIM4fgislGiven = true;
            break;
    case  BSIM4_MOD_KGISL :         /* v4.7 New GIDL/GISL */
            mod.BSIM4kgisl = value.get_as<double>().value();
            mod.BSIM4kgislGiven = true;
            break;
    case  BSIM4_MOD_RGISL :         /* v4.7 New GIDL/GISL */
            mod.BSIM4rgisl = value.get_as<double>().value();
            mod.BSIM4rgislGiven = true;
            break;
        case  BSIM4_MOD_AIGC :
            mod.BSIM4aigc = value.get_as<double>().value();
            mod.BSIM4aigcGiven = true;
            break;
        case  BSIM4_MOD_BIGC :
            mod.BSIM4bigc = value.get_as<double>().value();
            mod.BSIM4bigcGiven = true;
            break;
        case  BSIM4_MOD_CIGC :
            mod.BSIM4cigc = value.get_as<double>().value();
            mod.BSIM4cigcGiven = true;
            break;
        case  BSIM4_MOD_AIGSD :
            mod.BSIM4aigsd = value.get_as<double>().value();
            mod.BSIM4aigsdGiven = true;
            break;
        case  BSIM4_MOD_BIGSD :
            mod.BSIM4bigsd = value.get_as<double>().value();
            mod.BSIM4bigsdGiven = true;
            break;
        case  BSIM4_MOD_CIGSD :
            mod.BSIM4cigsd = value.get_as<double>().value();
            mod.BSIM4cigsdGiven = true;
            break;
        case  BSIM4_MOD_AIGS :
            mod.BSIM4aigs = value.get_as<double>().value();
            mod.BSIM4aigsGiven = true;
            break;
        case  BSIM4_MOD_BIGS :
            mod.BSIM4bigs = value.get_as<double>().value();
            mod.BSIM4bigsGiven = true;
            break;
        case  BSIM4_MOD_CIGS :
            mod.BSIM4cigs = value.get_as<double>().value();
            mod.BSIM4cigsGiven = true;
            break;
        case  BSIM4_MOD_AIGD :
            mod.BSIM4aigd = value.get_as<double>().value();
            mod.BSIM4aigdGiven = true;
            break;
        case  BSIM4_MOD_BIGD :
            mod.BSIM4bigd = value.get_as<double>().value();
            mod.BSIM4bigdGiven = true;
            break;
        case  BSIM4_MOD_CIGD :
            mod.BSIM4cigd = value.get_as<double>().value();
            mod.BSIM4cigdGiven = true;
            break;
        case  BSIM4_MOD_AIGBACC :
            mod.BSIM4aigbacc = value.get_as<double>().value();
            mod.BSIM4aigbaccGiven = true;
            break;
        case  BSIM4_MOD_BIGBACC :
            mod.BSIM4bigbacc = value.get_as<double>().value();
            mod.BSIM4bigbaccGiven = true;
            break;
        case  BSIM4_MOD_CIGBACC :
            mod.BSIM4cigbacc = value.get_as<double>().value();
            mod.BSIM4cigbaccGiven = true;
            break;
        case  BSIM4_MOD_AIGBINV :
            mod.BSIM4aigbinv = value.get_as<double>().value();
            mod.BSIM4aigbinvGiven = true;
            break;
        case  BSIM4_MOD_BIGBINV :
            mod.BSIM4bigbinv = value.get_as<double>().value();
            mod.BSIM4bigbinvGiven = true;
            break;
        case  BSIM4_MOD_CIGBINV :
            mod.BSIM4cigbinv = value.get_as<double>().value();
            mod.BSIM4cigbinvGiven = true;
            break;
        case  BSIM4_MOD_NIGC :
            mod.BSIM4nigc = value.get_as<double>().value();
            mod.BSIM4nigcGiven = true;
            break;
        case  BSIM4_MOD_NIGBINV :
            mod.BSIM4nigbinv = value.get_as<double>().value();
            mod.BSIM4nigbinvGiven = true;
            break;
        case  BSIM4_MOD_NIGBACC :
            mod.BSIM4nigbacc = value.get_as<double>().value();
            mod.BSIM4nigbaccGiven = true;
            break;
        case  BSIM4_MOD_NTOX :
            mod.BSIM4ntox = value.get_as<double>().value();
            mod.BSIM4ntoxGiven = true;
            break;
        case  BSIM4_MOD_EIGBINV :
            mod.BSIM4eigbinv = value.get_as<double>().value();
            mod.BSIM4eigbinvGiven = true;
            break;
        case  BSIM4_MOD_PIGCD :
            mod.BSIM4pigcd = value.get_as<double>().value();
            mod.BSIM4pigcdGiven = true;
            break;
        case  BSIM4_MOD_POXEDGE :
            mod.BSIM4poxedge = value.get_as<double>().value();
            mod.BSIM4poxedgeGiven = true;
            break;
        case  BSIM4_MOD_XRCRG1 :
            mod.BSIM4xrcrg1 = value.get_as<double>().value();
            mod.BSIM4xrcrg1Given = true;
            break;
        case  BSIM4_MOD_XRCRG2 :
            mod.BSIM4xrcrg2 = value.get_as<double>().value();
            mod.BSIM4xrcrg2Given = true;
            break;
        case  BSIM4_MOD_LAMBDA :
            mod.BSIM4lambda = value.get_as<double>().value();
            mod.BSIM4lambdaGiven = true;
            break;
        case  BSIM4_MOD_VTL :
            mod.BSIM4vtl = value.get_as<double>().value();
            mod.BSIM4vtlGiven = true;
            break;
        case  BSIM4_MOD_XN:
            mod.BSIM4xn = value.get_as<double>().value();
            mod.BSIM4xnGiven = true;
            break;
        case  BSIM4_MOD_LC:
            mod.BSIM4lc = value.get_as<double>().value();
            mod.BSIM4lcGiven = true;
            break;
        case  BSIM4_MOD_TNOIA :
            mod.BSIM4tnoia = value.get_as<double>().value();
            mod.BSIM4tnoiaGiven = true;
            break;
        case  BSIM4_MOD_TNOIB :
            mod.BSIM4tnoib = value.get_as<double>().value();
            mod.BSIM4tnoibGiven = true;
            break;
        case  BSIM4_MOD_TNOIC :
            mod.BSIM4tnoic = value.get_as<double>().value();
            mod.BSIM4tnoicGiven = true;
            break;
        case  BSIM4_MOD_RNOIA :
            mod.BSIM4rnoia = value.get_as<double>().value();
            mod.BSIM4rnoiaGiven = true;
            break;
        case  BSIM4_MOD_RNOIB :
            mod.BSIM4rnoib = value.get_as<double>().value();
            mod.BSIM4rnoibGiven = true;
            break;
        case  BSIM4_MOD_RNOIC :
            mod.BSIM4rnoic = value.get_as<double>().value();
            mod.BSIM4rnoicGiven = true;
            break;
        case  BSIM4_MOD_NTNOI :
            mod.BSIM4ntnoi = value.get_as<double>().value();
            mod.BSIM4ntnoiGiven = true;
            break;
        case  BSIM4_MOD_VFBSDOFF:
            mod.BSIM4vfbsdoff = value.get_as<double>().value();
            mod.BSIM4vfbsdoffGiven = true;
            break;
        case  BSIM4_MOD_TVFBSDOFF:
            mod.BSIM4tvfbsdoff = value.get_as<double>().value();
            mod.BSIM4tvfbsdoffGiven = true;
            break;
        case  BSIM4_MOD_LINTNOI:
            mod.BSIM4lintnoi = value.get_as<double>().value();
            mod.BSIM4lintnoiGiven = true;
            break;

        /* stress effect */
        case  BSIM4_MOD_SAREF :
            mod.BSIM4saref = value.get_as<double>().value();
            mod.BSIM4sarefGiven = true;
            break;
        case  BSIM4_MOD_SBREF :
            mod.BSIM4sbref = value.get_as<double>().value();
            mod.BSIM4sbrefGiven = true;
            break;
        case  BSIM4_MOD_WLOD :
            mod.BSIM4wlod = value.get_as<double>().value();
            mod.BSIM4wlodGiven = true;
            break;
        case  BSIM4_MOD_KU0 :
            mod.BSIM4ku0 = value.get_as<double>().value();
            mod.BSIM4ku0Given = true;
            break;
        case  BSIM4_MOD_KVSAT :
            mod.BSIM4kvsat = value.get_as<double>().value();
            mod.BSIM4kvsatGiven = true;
            break;
        case  BSIM4_MOD_KVTH0 :
            mod.BSIM4kvth0 = value.get_as<double>().value();
            mod.BSIM4kvth0Given = true;
            break;
        case  BSIM4_MOD_TKU0 :
            mod.BSIM4tku0 = value.get_as<double>().value();
            mod.BSIM4tku0Given = true;
            break;
        case  BSIM4_MOD_LLODKU0 :
            mod.BSIM4llodku0 = value.get_as<double>().value();
            mod.BSIM4llodku0Given = true;
            break;
        case  BSIM4_MOD_WLODKU0 :
            mod.BSIM4wlodku0 = value.get_as<double>().value();
            mod.BSIM4wlodku0Given = true;
            break;
        case  BSIM4_MOD_LLODVTH :
            mod.BSIM4llodvth = value.get_as<double>().value();
            mod.BSIM4llodvthGiven = true;
            break;
        case  BSIM4_MOD_WLODVTH :
            mod.BSIM4wlodvth = value.get_as<double>().value();
            mod.BSIM4wlodvthGiven = true;
            break;
        case  BSIM4_MOD_LKU0 :
            mod.BSIM4lku0 = value.get_as<double>().value();
            mod.BSIM4lku0Given = true;
            break;
        case  BSIM4_MOD_WKU0 :
            mod.BSIM4wku0 = value.get_as<double>().value();
            mod.BSIM4wku0Given = true;
            break;
        case  BSIM4_MOD_PKU0 :
            mod.BSIM4pku0 = value.get_as<double>().value();
            mod.BSIM4pku0Given = true;
            break;
        case  BSIM4_MOD_LKVTH0 :
            mod.BSIM4lkvth0 = value.get_as<double>().value();
            mod.BSIM4lkvth0Given = true;
            break;
        case  BSIM4_MOD_WKVTH0 :
            mod.BSIM4wkvth0 = value.get_as<double>().value();
            mod.BSIM4wkvth0Given = true;
            break;
        case  BSIM4_MOD_PKVTH0 :
            mod.BSIM4pkvth0 = value.get_as<double>().value();
            mod.BSIM4pkvth0Given = true;
            break;
        case  BSIM4_MOD_STK2 :
            mod.BSIM4stk2 = value.get_as<double>().value();
            mod.BSIM4stk2Given = true;
            break;
        case  BSIM4_MOD_LODK2 :
            mod.BSIM4lodk2 = value.get_as<double>().value();
            mod.BSIM4lodk2Given = true;
            break;
        case  BSIM4_MOD_STETA0 :
            mod.BSIM4steta0 = value.get_as<double>().value();
            mod.BSIM4steta0Given = true;
            break;
        case  BSIM4_MOD_LODETA0 :
            mod.BSIM4lodeta0 = value.get_as<double>().value();
            mod.BSIM4lodeta0Given = true;
            break;

        case  BSIM4_MOD_WEB :
            mod.BSIM4web = value.get_as<double>().value();
            mod.BSIM4webGiven = true;
            break;
    case BSIM4_MOD_WEC :
            mod.BSIM4wec = value.get_as<double>().value();
            mod.BSIM4wecGiven = true;
            break;
        case  BSIM4_MOD_KVTH0WE :
            mod.BSIM4kvth0we = value.get_as<double>().value();
            mod.BSIM4kvth0weGiven = true;
            break;
        case  BSIM4_MOD_K2WE :
            mod.BSIM4k2we = value.get_as<double>().value();
            mod.BSIM4k2weGiven = true;
            break;
        case  BSIM4_MOD_KU0WE :
            mod.BSIM4ku0we = value.get_as<double>().value();
            mod.BSIM4ku0weGiven = true;
            break;
        case  BSIM4_MOD_SCREF :
            mod.BSIM4scref = value.get_as<double>().value();
            mod.BSIM4screfGiven = true;
            break;
        case  BSIM4_MOD_WPEMOD :
            mod.BSIM4wpemod = value.get_as<double>().value();
            mod.BSIM4wpemodGiven = true;
            break;
        case  BSIM4_MOD_LKVTH0WE :
            mod.BSIM4lkvth0we = value.get_as<double>().value();
            mod.BSIM4lkvth0weGiven = true;
            break;
        case  BSIM4_MOD_LK2WE :
            mod.BSIM4lk2we = value.get_as<double>().value();
            mod.BSIM4lk2weGiven = true;
            break;
        case  BSIM4_MOD_LKU0WE :
            mod.BSIM4lku0we = value.get_as<double>().value();
            mod.BSIM4lku0weGiven = true;
            break;
        case  BSIM4_MOD_WKVTH0WE :
            mod.BSIM4wkvth0we = value.get_as<double>().value();
            mod.BSIM4wkvth0weGiven = true;
            break;
        case  BSIM4_MOD_WK2WE :
            mod.BSIM4wk2we = value.get_as<double>().value();
            mod.BSIM4wk2weGiven = true;
            break;
        case  BSIM4_MOD_WKU0WE :
            mod.BSIM4wku0we = value.get_as<double>().value();
            mod.BSIM4wku0weGiven = true;
            break;
        case  BSIM4_MOD_PKVTH0WE :
            mod.BSIM4pkvth0we = value.get_as<double>().value();
            mod.BSIM4pkvth0weGiven = true;
            break;
        case  BSIM4_MOD_PK2WE :
            mod.BSIM4pk2we = value.get_as<double>().value();
            mod.BSIM4pk2weGiven = true;
            break;
        case  BSIM4_MOD_PKU0WE :
            mod.BSIM4pku0we = value.get_as<double>().value();
            mod.BSIM4pku0weGiven = true;
            break;

        case  BSIM4_MOD_BETA0 :
            mod.BSIM4beta0 = value.get_as<double>().value();
            mod.BSIM4beta0Given = true;
            break;
        case  BSIM4_MOD_IJTHDFWD :
            mod.BSIM4ijthdfwd = value.get_as<double>().value();
            mod.BSIM4ijthdfwdGiven = true;
            break;
        case  BSIM4_MOD_IJTHSFWD :
            mod.BSIM4ijthsfwd = value.get_as<double>().value();
            mod.BSIM4ijthsfwdGiven = true;
            break;
        case  BSIM4_MOD_IJTHDREV :
            mod.BSIM4ijthdrev = value.get_as<double>().value();
            mod.BSIM4ijthdrevGiven = true;
            break;
        case  BSIM4_MOD_IJTHSREV :
            mod.BSIM4ijthsrev = value.get_as<double>().value();
            mod.BSIM4ijthsrevGiven = true;
            break;
        case  BSIM4_MOD_XJBVD :
            mod.BSIM4xjbvd = value.get_as<double>().value();
            mod.BSIM4xjbvdGiven = true;
            break;
        case  BSIM4_MOD_XJBVS :
            mod.BSIM4xjbvs = value.get_as<double>().value();
            mod.BSIM4xjbvsGiven = true;
            break;
        case  BSIM4_MOD_BVD :
            mod.BSIM4bvd = value.get_as<double>().value();
            mod.BSIM4bvdGiven = true;
            break;
        case  BSIM4_MOD_BVS :
            mod.BSIM4bvs = value.get_as<double>().value();
            mod.BSIM4bvsGiven = true;
            break;

        /* reverse diode */
        case  BSIM4_MOD_JTSS :
            mod.BSIM4jtss = value.get_as<double>().value();
            mod.BSIM4jtssGiven = true;
            break;
        case  BSIM4_MOD_JTSD :
            mod.BSIM4jtsd = value.get_as<double>().value();
            mod.BSIM4jtsdGiven = true;
            break;
        case  BSIM4_MOD_JTSSWS :
            mod.BSIM4jtssws = value.get_as<double>().value();
            mod.BSIM4jtsswsGiven = true;
            break;
        case  BSIM4_MOD_JTSSWD :
            mod.BSIM4jtsswd = value.get_as<double>().value();
            mod.BSIM4jtsswdGiven = true;
            break;
        case  BSIM4_MOD_JTSSWGS :
            mod.BSIM4jtsswgs = value.get_as<double>().value();
            mod.BSIM4jtsswgsGiven = true;
            break;
        case  BSIM4_MOD_JTSSWGD :
            mod.BSIM4jtsswgd = value.get_as<double>().value();
            mod.BSIM4jtsswgdGiven = true;
            break;

        case  BSIM4_MOD_JTWEFF :
         mod.BSIM4jtweff = value.get_as<double>().value();
         mod.BSIM4jtweffGiven = true;
             break;
        case  BSIM4_MOD_NJTS :
            mod.BSIM4njts = value.get_as<double>().value();
            mod.BSIM4njtsGiven = true;
            break;
        case  BSIM4_MOD_NJTSSW :
            mod.BSIM4njtssw = value.get_as<double>().value();
            mod.BSIM4njtsswGiven = true;
            break;
        case  BSIM4_MOD_NJTSSWG :
            mod.BSIM4njtsswg = value.get_as<double>().value();
            mod.BSIM4njtsswgGiven = true;
            break;
        case  BSIM4_MOD_NJTSD :
            mod.BSIM4njtsd = value.get_as<double>().value();
            mod.BSIM4njtsdGiven = true;
            break;
        case  BSIM4_MOD_NJTSSWD :
            mod.BSIM4njtsswd = value.get_as<double>().value();
            mod.BSIM4njtsswdGiven = true;
            break;
        case  BSIM4_MOD_NJTSSWGD :
            mod.BSIM4njtsswgd = value.get_as<double>().value();
            mod.BSIM4njtsswgdGiven = true;
            break;
        case  BSIM4_MOD_XTSS :
            mod.BSIM4xtss = value.get_as<double>().value();
            mod.BSIM4xtssGiven = true;
            break;
        case  BSIM4_MOD_XTSD :
            mod.BSIM4xtsd = value.get_as<double>().value();
            mod.BSIM4xtsdGiven = true;
            break;
        case  BSIM4_MOD_XTSSWS :
            mod.BSIM4xtssws = value.get_as<double>().value();
            mod.BSIM4xtsswsGiven = true;
            break;
        case  BSIM4_MOD_XTSSWD :
            mod.BSIM4xtsswd = value.get_as<double>().value();
            mod.BSIM4xtsswdGiven = true;
            break;
        case  BSIM4_MOD_XTSSWGS :
            mod.BSIM4xtsswgs = value.get_as<double>().value();
            mod.BSIM4xtsswgsGiven = true;
            break;
        case  BSIM4_MOD_XTSSWGD :
            mod.BSIM4xtsswgd = value.get_as<double>().value();
            mod.BSIM4xtsswgdGiven = true;
            break;
        case  BSIM4_MOD_TNJTS :
            mod.BSIM4tnjts = value.get_as<double>().value();
            mod.BSIM4tnjtsGiven = true;
            break;
        case  BSIM4_MOD_TNJTSSW :
            mod.BSIM4tnjtssw = value.get_as<double>().value();
            mod.BSIM4tnjtsswGiven = true;
            break;
        case  BSIM4_MOD_TNJTSSWG :
            mod.BSIM4tnjtsswg = value.get_as<double>().value();
            mod.BSIM4tnjtsswgGiven = true;
            break;
        case  BSIM4_MOD_TNJTSD :
            mod.BSIM4tnjtsd = value.get_as<double>().value();
            mod.BSIM4tnjtsdGiven = true;
            break;
        case  BSIM4_MOD_TNJTSSWD :
            mod.BSIM4tnjtsswd = value.get_as<double>().value();
            mod.BSIM4tnjtsswdGiven = true;
            break;
        case  BSIM4_MOD_TNJTSSWGD :
            mod.BSIM4tnjtsswgd = value.get_as<double>().value();
            mod.BSIM4tnjtsswgdGiven = true;
            break;
        case  BSIM4_MOD_VTSS :
            mod.BSIM4vtss = value.get_as<double>().value();
            mod.BSIM4vtssGiven = true;
            break;
        case  BSIM4_MOD_VTSD :
            mod.BSIM4vtsd = value.get_as<double>().value();
            mod.BSIM4vtsdGiven = true;
            break;
        case  BSIM4_MOD_VTSSWS :
            mod.BSIM4vtssws = value.get_as<double>().value();
            mod.BSIM4vtsswsGiven = true;
            break;
        case  BSIM4_MOD_VTSSWD :
            mod.BSIM4vtsswd = value.get_as<double>().value();
            mod.BSIM4vtsswdGiven = true;
            break;
        case  BSIM4_MOD_VTSSWGS :
            mod.BSIM4vtsswgs = value.get_as<double>().value();
            mod.BSIM4vtsswgsGiven = true;
            break;
        case  BSIM4_MOD_VTSSWGD :
            mod.BSIM4vtsswgd = value.get_as<double>().value();
            mod.BSIM4vtsswgdGiven = true;
            break;

        case  BSIM4_MOD_VFB :
            mod.BSIM4vfb = value.get_as<double>().value();
            mod.BSIM4vfbGiven = true;
            break;

        case  BSIM4_MOD_GBMIN :
            mod.BSIM4gbmin = value.get_as<double>().value();
            mod.BSIM4gbminGiven = true;
            break;
        case  BSIM4_MOD_RBDB :
            mod.BSIM4rbdb = value.get_as<double>().value();
            mod.BSIM4rbdbGiven = true;
            break;
        case  BSIM4_MOD_RBPB :
            mod.BSIM4rbpb = value.get_as<double>().value();
            mod.BSIM4rbpbGiven = true;
            break;
        case  BSIM4_MOD_RBSB :
            mod.BSIM4rbsb = value.get_as<double>().value();
            mod.BSIM4rbsbGiven = true;
            break;
        case  BSIM4_MOD_RBPS :
            mod.BSIM4rbps = value.get_as<double>().value();
            mod.BSIM4rbpsGiven = true;
            break;
        case  BSIM4_MOD_RBPD :
            mod.BSIM4rbpd = value.get_as<double>().value();
            mod.BSIM4rbpdGiven = true;
            break;

        case  BSIM4_MOD_RBPS0 :
            mod.BSIM4rbps0 = value.get_as<double>().value();
            mod.BSIM4rbps0Given = true;
            break;
        case  BSIM4_MOD_RBPSL :
            mod.BSIM4rbpsl = value.get_as<double>().value();
            mod.BSIM4rbpslGiven = true;
            break;
        case  BSIM4_MOD_RBPSW :
            mod.BSIM4rbpsw = value.get_as<double>().value();
            mod.BSIM4rbpswGiven = true;
            break;
        case  BSIM4_MOD_RBPSNF :
            mod.BSIM4rbpsnf = value.get_as<double>().value();
            mod.BSIM4rbpsnfGiven = true;
            break;

        case  BSIM4_MOD_RBPD0 :
            mod.BSIM4rbpd0 = value.get_as<double>().value();
            mod.BSIM4rbpd0Given = true;
            break;
        case  BSIM4_MOD_RBPDL :
            mod.BSIM4rbpdl = value.get_as<double>().value();
            mod.BSIM4rbpdlGiven = true;
            break;
        case  BSIM4_MOD_RBPDW :
            mod.BSIM4rbpdw = value.get_as<double>().value();
            mod.BSIM4rbpdwGiven = true;
            break;
        case  BSIM4_MOD_RBPDNF :
            mod.BSIM4rbpdnf = value.get_as<double>().value();
            mod.BSIM4rbpdnfGiven = true;
            break;

        case  BSIM4_MOD_RBPBX0 :
            mod.BSIM4rbpbx0 = value.get_as<double>().value();
            mod.BSIM4rbpbx0Given = true;
            break;
        case  BSIM4_MOD_RBPBXL :
            mod.BSIM4rbpbxl = value.get_as<double>().value();
            mod.BSIM4rbpbxlGiven = true;
            break;
        case  BSIM4_MOD_RBPBXW :
            mod.BSIM4rbpbxw = value.get_as<double>().value();
            mod.BSIM4rbpbxwGiven = true;
            break;
        case  BSIM4_MOD_RBPBXNF :
            mod.BSIM4rbpbxnf = value.get_as<double>().value();
            mod.BSIM4rbpbxnfGiven = true;
            break;
        case  BSIM4_MOD_RBPBY0 :
            mod.BSIM4rbpby0 = value.get_as<double>().value();
            mod.BSIM4rbpby0Given = true;
            break;
        case  BSIM4_MOD_RBPBYL :
            mod.BSIM4rbpbyl = value.get_as<double>().value();
            mod.BSIM4rbpbylGiven = true;
            break;
        case  BSIM4_MOD_RBPBYW :
            mod.BSIM4rbpbyw = value.get_as<double>().value();
            mod.BSIM4rbpbywGiven = true;
            break;
        case  BSIM4_MOD_RBPBYNF :
            mod.BSIM4rbpbynf = value.get_as<double>().value();
            mod.BSIM4rbpbynfGiven = true;
            break;
       case  BSIM4_MOD_RBSBX0 :
            mod.BSIM4rbsbx0 = value.get_as<double>().value();
            mod.BSIM4rbsbx0Given = true;
            break;
       case  BSIM4_MOD_RBSBY0 :
            mod.BSIM4rbsby0 = value.get_as<double>().value();
            mod.BSIM4rbsby0Given = true;
            break;
       case  BSIM4_MOD_RBDBX0 :
            mod.BSIM4rbdbx0 = value.get_as<double>().value();
            mod.BSIM4rbdbx0Given = true;
            break;
       case  BSIM4_MOD_RBDBY0 :
            mod.BSIM4rbdby0 = value.get_as<double>().value();
            mod.BSIM4rbdby0Given = true;
            break;


       case  BSIM4_MOD_RBSDBXL :
            mod.BSIM4rbsdbxl = value.get_as<double>().value();
            mod.BSIM4rbsdbxlGiven = true;
            break;
       case  BSIM4_MOD_RBSDBXW :
            mod.BSIM4rbsdbxw = value.get_as<double>().value();
            mod.BSIM4rbsdbxwGiven = true;
            break;
       case  BSIM4_MOD_RBSDBXNF :
            mod.BSIM4rbsdbxnf = value.get_as<double>().value();
            mod.BSIM4rbsdbxnfGiven = true;
            break;
       case  BSIM4_MOD_RBSDBYL :
            mod.BSIM4rbsdbyl = value.get_as<double>().value();
            mod.BSIM4rbsdbylGiven = true;
            break;
       case  BSIM4_MOD_RBSDBYW :
            mod.BSIM4rbsdbyw = value.get_as<double>().value();
            mod.BSIM4rbsdbywGiven = true;
            break;
       case  BSIM4_MOD_RBSDBYNF :
            mod.BSIM4rbsdbynf = value.get_as<double>().value();
            mod.BSIM4rbsdbynfGiven = true;
            break;

        case  BSIM4_MOD_CGSL :
            mod.BSIM4cgsl = value.get_as<double>().value();
            mod.BSIM4cgslGiven = true;
            break;
        case  BSIM4_MOD_CGDL :
            mod.BSIM4cgdl = value.get_as<double>().value();
            mod.BSIM4cgdlGiven = true;
            break;
        case  BSIM4_MOD_CKAPPAS :
            mod.BSIM4ckappas = value.get_as<double>().value();
            mod.BSIM4ckappasGiven = true;
            break;
        case  BSIM4_MOD_CKAPPAD :
            mod.BSIM4ckappad = value.get_as<double>().value();
            mod.BSIM4ckappadGiven = true;
            break;
        case  BSIM4_MOD_CF :
            mod.BSIM4cf = value.get_as<double>().value();
            mod.BSIM4cfGiven = true;
            break;
        case  BSIM4_MOD_CLC :
            mod.BSIM4clc = value.get_as<double>().value();
            mod.BSIM4clcGiven = true;
            break;
        case  BSIM4_MOD_CLE :
            mod.BSIM4cle = value.get_as<double>().value();
            mod.BSIM4cleGiven = true;
            break;
        case  BSIM4_MOD_DWC :
            mod.BSIM4dwc = value.get_as<double>().value();
            mod.BSIM4dwcGiven = true;
            break;
        case  BSIM4_MOD_DLC :
            mod.BSIM4dlc = value.get_as<double>().value();
            mod.BSIM4dlcGiven = true;
            break;
        case  BSIM4_MOD_XW :
            mod.BSIM4xw = value.get_as<double>().value();
            mod.BSIM4xwGiven = true;
            break;
        case  BSIM4_MOD_XL :
            mod.BSIM4xl = value.get_as<double>().value();
            mod.BSIM4xlGiven = true;
            break;
        case  BSIM4_MOD_DLCIG :
            mod.BSIM4dlcig = value.get_as<double>().value();
            mod.BSIM4dlcigGiven = true;
            break;
        case  BSIM4_MOD_DLCIGD :
            mod.BSIM4dlcigd = value.get_as<double>().value();
            mod.BSIM4dlcigdGiven = true;
            break;
        case  BSIM4_MOD_DWJ :
            mod.BSIM4dwj = value.get_as<double>().value();
            mod.BSIM4dwjGiven = true;
            break;
        case  BSIM4_MOD_VFBCV :
            mod.BSIM4vfbcv = value.get_as<double>().value();
            mod.BSIM4vfbcvGiven = true;
            break;
        case  BSIM4_MOD_ACDE :
            mod.BSIM4acde = value.get_as<double>().value();
            mod.BSIM4acdeGiven = true;
            break;
        case  BSIM4_MOD_MOIN :
            mod.BSIM4moin = value.get_as<double>().value();
            mod.BSIM4moinGiven = true;
            break;
        case  BSIM4_MOD_NOFF :
            mod.BSIM4noff = value.get_as<double>().value();
            mod.BSIM4noffGiven = true;
            break;
        case  BSIM4_MOD_VOFFCV :
            mod.BSIM4voffcv = value.get_as<double>().value();
            mod.BSIM4voffcvGiven = true;
            break;
        case  BSIM4_MOD_DMCG :
            mod.BSIM4dmcg = value.get_as<double>().value();
            mod.BSIM4dmcgGiven = true;
            break;
        case  BSIM4_MOD_DMCI :
            mod.BSIM4dmci = value.get_as<double>().value();
            mod.BSIM4dmciGiven = true;
            break;
        case  BSIM4_MOD_DMDG :
            mod.BSIM4dmdg = value.get_as<double>().value();
            mod.BSIM4dmdgGiven = true;
            break;
        case  BSIM4_MOD_DMCGT :
            mod.BSIM4dmcgt = value.get_as<double>().value();
            mod.BSIM4dmcgtGiven = true;
            break;
        case  BSIM4_MOD_XGW :
            mod.BSIM4xgw = value.get_as<double>().value();
            mod.BSIM4xgwGiven = true;
            break;
        case  BSIM4_MOD_XGL :
            mod.BSIM4xgl = value.get_as<double>().value();
            mod.BSIM4xglGiven = true;
            break;
        case  BSIM4_MOD_RSHG :
            mod.BSIM4rshg = value.get_as<double>().value();
            mod.BSIM4rshgGiven = true;
            break;
        case  BSIM4_MOD_NGCON :
            mod.BSIM4ngcon = value.get_as<double>().value();
            mod.BSIM4ngconGiven = true;
            break;
        case  BSIM4_MOD_TCJ :
            mod.BSIM4tcj = value.get_as<double>().value();
            mod.BSIM4tcjGiven = true;
            break;
        case  BSIM4_MOD_TPB :
            mod.BSIM4tpb = value.get_as<double>().value();
            mod.BSIM4tpbGiven = true;
            break;
        case  BSIM4_MOD_TCJSW :
            mod.BSIM4tcjsw = value.get_as<double>().value();
            mod.BSIM4tcjswGiven = true;
            break;
        case  BSIM4_MOD_TPBSW :
            mod.BSIM4tpbsw = value.get_as<double>().value();
            mod.BSIM4tpbswGiven = true;
            break;
        case  BSIM4_MOD_TCJSWG :
            mod.BSIM4tcjswg = value.get_as<double>().value();
            mod.BSIM4tcjswgGiven = true;
            break;
        case  BSIM4_MOD_TPBSWG :
            mod.BSIM4tpbswg = value.get_as<double>().value();
            mod.BSIM4tpbswgGiven = true;
            break;

    /* Length dependence */
        case  BSIM4_MOD_LCDSC :
            mod.BSIM4lcdsc = value.get_as<double>().value();
            mod.BSIM4lcdscGiven = true;
            break;
        case  BSIM4_MOD_LCDSCB :
            mod.BSIM4lcdscb = value.get_as<double>().value();
            mod.BSIM4lcdscbGiven = true;
            break;
        case  BSIM4_MOD_LCDSCD :
            mod.BSIM4lcdscd = value.get_as<double>().value();
            mod.BSIM4lcdscdGiven = true;
            break;
        case  BSIM4_MOD_LCIT :
            mod.BSIM4lcit = value.get_as<double>().value();
            mod.BSIM4lcitGiven = true;
            break;
        case  BSIM4_MOD_LNFACTOR :
            mod.BSIM4lnfactor = value.get_as<double>().value();
            mod.BSIM4lnfactorGiven = true;
            break;
        case BSIM4_MOD_LXJ:
            mod.BSIM4lxj = value.get_as<double>().value();
            mod.BSIM4lxjGiven = true;
            break;
        case BSIM4_MOD_LVSAT:
            mod.BSIM4lvsat = value.get_as<double>().value();
            mod.BSIM4lvsatGiven = true;
            break;


        case BSIM4_MOD_LA0:
            mod.BSIM4la0 = value.get_as<double>().value();
            mod.BSIM4la0Given = true;
            break;
        case BSIM4_MOD_LAGS:
            mod.BSIM4lags = value.get_as<double>().value();
            mod.BSIM4lagsGiven = true;
            break;
        case BSIM4_MOD_LA1:
            mod.BSIM4la1 = value.get_as<double>().value();
            mod.BSIM4la1Given = true;
            break;
        case BSIM4_MOD_LA2:
            mod.BSIM4la2 = value.get_as<double>().value();
            mod.BSIM4la2Given = true;
            break;
        case BSIM4_MOD_LAT:
            mod.BSIM4lat = value.get_as<double>().value();
            mod.BSIM4latGiven = true;
            break;
        case BSIM4_MOD_LKETA:
            mod.BSIM4lketa = value.get_as<double>().value();
            mod.BSIM4lketaGiven = true;
            break;
        case BSIM4_MOD_LNSUB:
            mod.BSIM4lnsub = value.get_as<double>().value();
            mod.BSIM4lnsubGiven = true;
            break;
        case BSIM4_MOD_LNDEP:
            mod.BSIM4lndep = value.get_as<double>().value();
            mod.BSIM4lndepGiven = true;
        if (mod.BSIM4lndep > 1.0e20)
        mod.BSIM4lndep *= 1.0e-6;
            break;
        case BSIM4_MOD_LNSD:
            mod.BSIM4lnsd = value.get_as<double>().value();
            mod.BSIM4lnsdGiven = true;
            if (mod.BSIM4lnsd > 1.0e23)
                mod.BSIM4lnsd *= 1.0e-6;
            break;
        case BSIM4_MOD_LNGATE:
            mod.BSIM4lngate = value.get_as<double>().value();
            mod.BSIM4lngateGiven = true;
        if (mod.BSIM4lngate > 1.0e23)
        mod.BSIM4lngate *= 1.0e-6;
            break;
        case BSIM4_MOD_LGAMMA1:
            mod.BSIM4lgamma1 = value.get_as<double>().value();
            mod.BSIM4lgamma1Given = true;
            break;
        case BSIM4_MOD_LGAMMA2:
            mod.BSIM4lgamma2 = value.get_as<double>().value();
            mod.BSIM4lgamma2Given = true;
            break;
        case BSIM4_MOD_LVBX:
            mod.BSIM4lvbx = value.get_as<double>().value();
            mod.BSIM4lvbxGiven = true;
            break;
        case BSIM4_MOD_LVBM:
            mod.BSIM4lvbm = value.get_as<double>().value();
            mod.BSIM4lvbmGiven = true;
            break;
        case BSIM4_MOD_LXT:
            mod.BSIM4lxt = value.get_as<double>().value();
            mod.BSIM4lxtGiven = true;
            break;
        case  BSIM4_MOD_LK1:
            mod.BSIM4lk1 = value.get_as<double>().value();
            mod.BSIM4lk1Given = true;
            break;
        case  BSIM4_MOD_LKT1:
            mod.BSIM4lkt1 = value.get_as<double>().value();
            mod.BSIM4lkt1Given = true;
            break;
        case  BSIM4_MOD_LKT1L:
            mod.BSIM4lkt1l = value.get_as<double>().value();
            mod.BSIM4lkt1lGiven = true;
            break;
        case  BSIM4_MOD_LKT2:
            mod.BSIM4lkt2 = value.get_as<double>().value();
            mod.BSIM4lkt2Given = true;
            break;
        case  BSIM4_MOD_LK2:
            mod.BSIM4lk2 = value.get_as<double>().value();
            mod.BSIM4lk2Given = true;
            break;
        case  BSIM4_MOD_LK3:
            mod.BSIM4lk3 = value.get_as<double>().value();
            mod.BSIM4lk3Given = true;
            break;
        case  BSIM4_MOD_LK3B:
            mod.BSIM4lk3b = value.get_as<double>().value();
            mod.BSIM4lk3bGiven = true;
            break;
        case  BSIM4_MOD_LLPE0:
            mod.BSIM4llpe0 = value.get_as<double>().value();
            mod.BSIM4llpe0Given = true;
            break;
        case  BSIM4_MOD_LLPEB:
            mod.BSIM4llpeb = value.get_as<double>().value();
            mod.BSIM4llpebGiven = true;
            break;
        case  BSIM4_MOD_LDVTP0:
            mod.BSIM4ldvtp0 = value.get_as<double>().value();
            mod.BSIM4ldvtp0Given = true;
            break;
        case  BSIM4_MOD_LDVTP1:
            mod.BSIM4ldvtp1 = value.get_as<double>().value();
            mod.BSIM4ldvtp1Given = true;
            break;
    case  BSIM4_MOD_LDVTP2:     /* New DIBL/Rout */
            mod.BSIM4ldvtp2 = value.get_as<double>().value();
            mod.BSIM4ldvtp2Given = true;
            break;
    case  BSIM4_MOD_LDVTP3:
            mod.BSIM4ldvtp3 = value.get_as<double>().value();
            mod.BSIM4ldvtp3Given = true;
            break;
    case  BSIM4_MOD_LDVTP4:
            mod.BSIM4ldvtp4 = value.get_as<double>().value();
            mod.BSIM4ldvtp4Given = true;
            break;
     case  BSIM4_MOD_LDVTP5:
            mod.BSIM4ldvtp5 = value.get_as<double>().value();
            mod.BSIM4ldvtp5Given = true;
            break;
       case  BSIM4_MOD_LW0:
            mod.BSIM4lw0 = value.get_as<double>().value();
            mod.BSIM4lw0Given = true;
            break;
        case  BSIM4_MOD_LDVT0:
            mod.BSIM4ldvt0 = value.get_as<double>().value();
            mod.BSIM4ldvt0Given = true;
            break;
        case  BSIM4_MOD_LDVT1:
            mod.BSIM4ldvt1 = value.get_as<double>().value();
            mod.BSIM4ldvt1Given = true;
            break;
        case  BSIM4_MOD_LDVT2:
            mod.BSIM4ldvt2 = value.get_as<double>().value();
            mod.BSIM4ldvt2Given = true;
            break;
        case  BSIM4_MOD_LDVT0W:
            mod.BSIM4ldvt0w = value.get_as<double>().value();
            mod.BSIM4ldvt0wGiven = true;
            break;
        case  BSIM4_MOD_LDVT1W:
            mod.BSIM4ldvt1w = value.get_as<double>().value();
            mod.BSIM4ldvt1wGiven = true;
            break;
        case  BSIM4_MOD_LDVT2W:
            mod.BSIM4ldvt2w = value.get_as<double>().value();
            mod.BSIM4ldvt2wGiven = true;
            break;
        case  BSIM4_MOD_LDROUT:
            mod.BSIM4ldrout = value.get_as<double>().value();
            mod.BSIM4ldroutGiven = true;
            break;
        case  BSIM4_MOD_LDSUB:
            mod.BSIM4ldsub = value.get_as<double>().value();
            mod.BSIM4ldsubGiven = true;
            break;
        case BSIM4_MOD_LVTH0:
            mod.BSIM4lvth0 = value.get_as<double>().value();
            mod.BSIM4lvth0Given = true;
            break;
        case BSIM4_MOD_LUA:
            mod.BSIM4lua = value.get_as<double>().value();
            mod.BSIM4luaGiven = true;
            break;
        case BSIM4_MOD_LUA1:
            mod.BSIM4lua1 = value.get_as<double>().value();
            mod.BSIM4lua1Given = true;
            break;
        case BSIM4_MOD_LUB:
            mod.BSIM4lub = value.get_as<double>().value();
            mod.BSIM4lubGiven = true;
            break;
        case BSIM4_MOD_LUB1:
            mod.BSIM4lub1 = value.get_as<double>().value();
            mod.BSIM4lub1Given = true;
            break;
        case BSIM4_MOD_LUC:
            mod.BSIM4luc = value.get_as<double>().value();
            mod.BSIM4lucGiven = true;
            break;
        case BSIM4_MOD_LUC1:
            mod.BSIM4luc1 = value.get_as<double>().value();
            mod.BSIM4luc1Given = true;
            break;
        case  BSIM4_MOD_LU0 :
            mod.BSIM4lu0 = value.get_as<double>().value();
            mod.BSIM4lu0Given = true;
            break;
        case  BSIM4_MOD_LUTE :
            mod.BSIM4lute = value.get_as<double>().value();
            mod.BSIM4luteGiven = true;
            break;
        case  BSIM4_MOD_LUCSTE :
            mod.BSIM4lucste = value.get_as<double>().value();
            mod.BSIM4lucsteGiven = true;
            break;
        case BSIM4_MOD_LVOFF:
            mod.BSIM4lvoff = value.get_as<double>().value();
            mod.BSIM4lvoffGiven = true;
            break;
        case BSIM4_MOD_LTVOFF:
            mod.BSIM4ltvoff = value.get_as<double>().value();
            mod.BSIM4ltvoffGiven = true;
            break;
        case BSIM4_MOD_LTNFACTOR:       /* v4.7 temp dep of leakage current  */
            mod.BSIM4ltnfactor = value.get_as<double>().value();
            mod.BSIM4ltnfactorGiven = true;
            break;
        case BSIM4_MOD_LTETA0:      /* v4.7 temp dep of leakage current  */
            mod.BSIM4lteta0 = value.get_as<double>().value();
            mod.BSIM4lteta0Given = true;
            break;
        case BSIM4_MOD_LTVOFFCV:    /* v4.7 temp dep of leakage current  */
            mod.BSIM4ltvoffcv = value.get_as<double>().value();
            mod.BSIM4ltvoffcvGiven = true;
            break;
        case BSIM4_MOD_LMINV:
            mod.BSIM4lminv = value.get_as<double>().value();
            mod.BSIM4lminvGiven = true;
            break;
        case BSIM4_MOD_LMINVCV:
            mod.BSIM4lminvcv = value.get_as<double>().value();
            mod.BSIM4lminvcvGiven = true;
            break;
        case BSIM4_MOD_LFPROUT:
            mod.BSIM4lfprout = value.get_as<double>().value();
            mod.BSIM4lfproutGiven = true;
            break;
        case BSIM4_MOD_LPDITS:
            mod.BSIM4lpdits = value.get_as<double>().value();
            mod.BSIM4lpditsGiven = true;
            break;
        case BSIM4_MOD_LPDITSD:
            mod.BSIM4lpditsd = value.get_as<double>().value();
            mod.BSIM4lpditsdGiven = true;
            break;
        case  BSIM4_MOD_LDELTA :
            mod.BSIM4ldelta = value.get_as<double>().value();
            mod.BSIM4ldeltaGiven = true;
            break;
        case BSIM4_MOD_LRDSW:
            mod.BSIM4lrdsw = value.get_as<double>().value();
            mod.BSIM4lrdswGiven = true;
            break;
        case BSIM4_MOD_LRDW:
            mod.BSIM4lrdw = value.get_as<double>().value();
            mod.BSIM4lrdwGiven = true;
            break;
        case BSIM4_MOD_LRSW:
            mod.BSIM4lrsw = value.get_as<double>().value();
            mod.BSIM4lrswGiven = true;
            break;
        case BSIM4_MOD_LPRWB:
            mod.BSIM4lprwb = value.get_as<double>().value();
            mod.BSIM4lprwbGiven = true;
            break;
        case BSIM4_MOD_LPRWG:
            mod.BSIM4lprwg = value.get_as<double>().value();
            mod.BSIM4lprwgGiven = true;
            break;
        case BSIM4_MOD_LPRT:
            mod.BSIM4lprt = value.get_as<double>().value();
            mod.BSIM4lprtGiven = true;
            break;
        case BSIM4_MOD_LETA0:
            mod.BSIM4leta0 = value.get_as<double>().value();
            mod.BSIM4leta0Given = true;
            break;
        case BSIM4_MOD_LETAB:
            mod.BSIM4letab = value.get_as<double>().value();
            mod.BSIM4letabGiven = true;
            break;
        case BSIM4_MOD_LPCLM:
            mod.BSIM4lpclm = value.get_as<double>().value();
            mod.BSIM4lpclmGiven = true;
            break;
        case BSIM4_MOD_LPDIBL1:
            mod.BSIM4lpdibl1 = value.get_as<double>().value();
            mod.BSIM4lpdibl1Given = true;
            break;
        case BSIM4_MOD_LPDIBL2:
            mod.BSIM4lpdibl2 = value.get_as<double>().value();
            mod.BSIM4lpdibl2Given = true;
            break;
        case BSIM4_MOD_LPDIBLB:
            mod.BSIM4lpdiblb = value.get_as<double>().value();
            mod.BSIM4lpdiblbGiven = true;
            break;
        case BSIM4_MOD_LPSCBE1:
            mod.BSIM4lpscbe1 = value.get_as<double>().value();
            mod.BSIM4lpscbe1Given = true;
            break;
        case BSIM4_MOD_LPSCBE2:
            mod.BSIM4lpscbe2 = value.get_as<double>().value();
            mod.BSIM4lpscbe2Given = true;
            break;
        case BSIM4_MOD_LPVAG:
            mod.BSIM4lpvag = value.get_as<double>().value();
            mod.BSIM4lpvagGiven = true;
            break;
        case  BSIM4_MOD_LWR :
            mod.BSIM4lwr = value.get_as<double>().value();
            mod.BSIM4lwrGiven = true;
            break;
        case  BSIM4_MOD_LDWG :
            mod.BSIM4ldwg = value.get_as<double>().value();
            mod.BSIM4ldwgGiven = true;
            break;
        case  BSIM4_MOD_LDWB :
            mod.BSIM4ldwb = value.get_as<double>().value();
            mod.BSIM4ldwbGiven = true;
            break;
        case  BSIM4_MOD_LB0 :
            mod.BSIM4lb0 = value.get_as<double>().value();
            mod.BSIM4lb0Given = true;
            break;
        case  BSIM4_MOD_LB1 :
            mod.BSIM4lb1 = value.get_as<double>().value();
            mod.BSIM4lb1Given = true;
            break;
        case  BSIM4_MOD_LALPHA0 :
            mod.BSIM4lalpha0 = value.get_as<double>().value();
            mod.BSIM4lalpha0Given = true;
            break;
        case  BSIM4_MOD_LALPHA1 :
            mod.BSIM4lalpha1 = value.get_as<double>().value();
            mod.BSIM4lalpha1Given = true;
            break;
        case  BSIM4_MOD_LBETA0 :
            mod.BSIM4lbeta0 = value.get_as<double>().value();
            mod.BSIM4lbeta0Given = true;
            break;
        case  BSIM4_MOD_LPHIN :
            mod.BSIM4lphin = value.get_as<double>().value();
            mod.BSIM4lphinGiven = true;
            break;
        case  BSIM4_MOD_LAGIDL :
            mod.BSIM4lagidl = value.get_as<double>().value();
            mod.BSIM4lagidlGiven = true;
            break;
        case  BSIM4_MOD_LBGIDL :
            mod.BSIM4lbgidl = value.get_as<double>().value();
            mod.BSIM4lbgidlGiven = true;
            break;
        case  BSIM4_MOD_LCGIDL :
            mod.BSIM4lcgidl = value.get_as<double>().value();
            mod.BSIM4lcgidlGiven = true;
            break;
        case  BSIM4_MOD_LEGIDL :
            mod.BSIM4legidl = value.get_as<double>().value();
            mod.BSIM4legidlGiven = true;
            break;
    case  BSIM4_MOD_LFGIDL :            /* v4.7 New GIDL/GISL */
            mod.BSIM4lfgidl = value.get_as<double>().value();
            mod.BSIM4lfgidlGiven = true;
            break;
    case  BSIM4_MOD_LKGIDL :            /* v4.7 New GIDL/GISL */
            mod.BSIM4lkgidl = value.get_as<double>().value();
            mod.BSIM4lkgidlGiven = true;
            break;
    case  BSIM4_MOD_LRGIDL :            /* v4.7 New GIDL/GISL */
            mod.BSIM4lrgidl = value.get_as<double>().value();
            mod.BSIM4lrgidlGiven = true;
            break;
        case  BSIM4_MOD_LAGISL :
            mod.BSIM4lagisl = value.get_as<double>().value();
            mod.BSIM4lagislGiven = true;
            break;
        case  BSIM4_MOD_LBGISL :
            mod.BSIM4lbgisl = value.get_as<double>().value();
            mod.BSIM4lbgislGiven = true;
            break;
        case  BSIM4_MOD_LCGISL :
            mod.BSIM4lcgisl = value.get_as<double>().value();
            mod.BSIM4lcgislGiven = true;
            break;
        case  BSIM4_MOD_LEGISL :
            mod.BSIM4legisl = value.get_as<double>().value();
            mod.BSIM4legislGiven = true;
            break;
    case  BSIM4_MOD_LFGISL :            /* v4.7 New GIDL/GISL */
            mod.BSIM4lfgisl = value.get_as<double>().value();
            mod.BSIM4lfgislGiven = true;
            break;
    case  BSIM4_MOD_LKGISL :            /* v4.7 New GIDL/GISL */
            mod.BSIM4lkgisl = value.get_as<double>().value();
            mod.BSIM4lkgislGiven = true;
            break;
    case  BSIM4_MOD_LRGISL :            /* v4.7 New GIDL/GISL */
            mod.BSIM4lrgisl = value.get_as<double>().value();
            mod.BSIM4lrgislGiven = true;
            break;
        case  BSIM4_MOD_LAIGC :
            mod.BSIM4laigc = value.get_as<double>().value();
            mod.BSIM4laigcGiven = true;
            break;
        case  BSIM4_MOD_LBIGC :
            mod.BSIM4lbigc = value.get_as<double>().value();
            mod.BSIM4lbigcGiven = true;
            break;
        case  BSIM4_MOD_LCIGC :
            mod.BSIM4lcigc = value.get_as<double>().value();
            mod.BSIM4lcigcGiven = true;
            break;
        case  BSIM4_MOD_LAIGSD :
            mod.BSIM4laigsd = value.get_as<double>().value();
            mod.BSIM4laigsdGiven = true;
            break;
        case  BSIM4_MOD_LBIGSD :
            mod.BSIM4lbigsd = value.get_as<double>().value();
            mod.BSIM4lbigsdGiven = true;
            break;
        case  BSIM4_MOD_LCIGSD :
            mod.BSIM4lcigsd = value.get_as<double>().value();
            mod.BSIM4lcigsdGiven = true;
            break;
        case  BSIM4_MOD_LAIGS :
            mod.BSIM4laigs = value.get_as<double>().value();
            mod.BSIM4laigsGiven = true;
            break;
        case  BSIM4_MOD_LBIGS :
            mod.BSIM4lbigs = value.get_as<double>().value();
            mod.BSIM4lbigsGiven = true;
            break;
        case  BSIM4_MOD_LCIGS :
            mod.BSIM4lcigs = value.get_as<double>().value();
            mod.BSIM4lcigsGiven = true;
            break;
        case  BSIM4_MOD_LAIGD :
            mod.BSIM4laigd = value.get_as<double>().value();
            mod.BSIM4laigdGiven = true;
            break;
        case  BSIM4_MOD_LBIGD :
            mod.BSIM4lbigd = value.get_as<double>().value();
            mod.BSIM4lbigdGiven = true;
            break;
        case  BSIM4_MOD_LCIGD :
            mod.BSIM4lcigd = value.get_as<double>().value();
            mod.BSIM4lcigdGiven = true;
            break;
        case  BSIM4_MOD_LAIGBACC :
            mod.BSIM4laigbacc = value.get_as<double>().value();
            mod.BSIM4laigbaccGiven = true;
            break;
        case  BSIM4_MOD_LBIGBACC :
            mod.BSIM4lbigbacc = value.get_as<double>().value();
            mod.BSIM4lbigbaccGiven = true;
            break;
        case  BSIM4_MOD_LCIGBACC :
            mod.BSIM4lcigbacc = value.get_as<double>().value();
            mod.BSIM4lcigbaccGiven = true;
            break;
        case  BSIM4_MOD_LAIGBINV :
            mod.BSIM4laigbinv = value.get_as<double>().value();
            mod.BSIM4laigbinvGiven = true;
            break;
        case  BSIM4_MOD_LBIGBINV :
            mod.BSIM4lbigbinv = value.get_as<double>().value();
            mod.BSIM4lbigbinvGiven = true;
            break;
        case  BSIM4_MOD_LCIGBINV :
            mod.BSIM4lcigbinv = value.get_as<double>().value();
            mod.BSIM4lcigbinvGiven = true;
            break;
        case  BSIM4_MOD_LNIGC :
            mod.BSIM4lnigc = value.get_as<double>().value();
            mod.BSIM4lnigcGiven = true;
            break;
        case  BSIM4_MOD_LNIGBINV :
            mod.BSIM4lnigbinv = value.get_as<double>().value();
            mod.BSIM4lnigbinvGiven = true;
            break;
        case  BSIM4_MOD_LNIGBACC :
            mod.BSIM4lnigbacc = value.get_as<double>().value();
            mod.BSIM4lnigbaccGiven = true;
            break;
        case  BSIM4_MOD_LNTOX :
            mod.BSIM4lntox = value.get_as<double>().value();
            mod.BSIM4lntoxGiven = true;
            break;
        case  BSIM4_MOD_LEIGBINV :
            mod.BSIM4leigbinv = value.get_as<double>().value();
            mod.BSIM4leigbinvGiven = true;
            break;
        case  BSIM4_MOD_LPIGCD :
            mod.BSIM4lpigcd = value.get_as<double>().value();
            mod.BSIM4lpigcdGiven = true;
            break;
        case  BSIM4_MOD_LPOXEDGE :
            mod.BSIM4lpoxedge = value.get_as<double>().value();
            mod.BSIM4lpoxedgeGiven = true;
            break;
        case  BSIM4_MOD_LXRCRG1 :
            mod.BSIM4lxrcrg1 = value.get_as<double>().value();
            mod.BSIM4lxrcrg1Given = true;
            break;
        case  BSIM4_MOD_LXRCRG2 :
            mod.BSIM4lxrcrg2 = value.get_as<double>().value();
            mod.BSIM4lxrcrg2Given = true;
            break;
        case  BSIM4_MOD_LLAMBDA :
            mod.BSIM4llambda = value.get_as<double>().value();
            mod.BSIM4llambdaGiven = true;
            break;
        case  BSIM4_MOD_LVTL :
            mod.BSIM4lvtl = value.get_as<double>().value();
            mod.BSIM4lvtlGiven = true;
            break;
        case  BSIM4_MOD_LXN:
            mod.BSIM4lxn = value.get_as<double>().value();
            mod.BSIM4lxnGiven = true;
            break;
        case  BSIM4_MOD_LVFBSDOFF:
            mod.BSIM4lvfbsdoff = value.get_as<double>().value();
            mod.BSIM4lvfbsdoffGiven = true;
            break;
        case  BSIM4_MOD_LTVFBSDOFF:
            mod.BSIM4ltvfbsdoff = value.get_as<double>().value();
            mod.BSIM4ltvfbsdoffGiven = true;
            break;
        case  BSIM4_MOD_LEU :
            mod.BSIM4leu = value.get_as<double>().value();
            mod.BSIM4leuGiven = true;
            break;
        case  BSIM4_MOD_LUCS :
            mod.BSIM4lucs = value.get_as<double>().value();
            mod.BSIM4lucsGiven = true;
            break;
        case  BSIM4_MOD_LVFB :
            mod.BSIM4lvfb = value.get_as<double>().value();
            mod.BSIM4lvfbGiven = true;
            break;
        case  BSIM4_MOD_LCGSL :
            mod.BSIM4lcgsl = value.get_as<double>().value();
            mod.BSIM4lcgslGiven = true;
            break;
        case  BSIM4_MOD_LCGDL :
            mod.BSIM4lcgdl = value.get_as<double>().value();
            mod.BSIM4lcgdlGiven = true;
            break;
        case  BSIM4_MOD_LCKAPPAS :
            mod.BSIM4lckappas = value.get_as<double>().value();
            mod.BSIM4lckappasGiven = true;
            break;
        case  BSIM4_MOD_LCKAPPAD :
            mod.BSIM4lckappad = value.get_as<double>().value();
            mod.BSIM4lckappadGiven = true;
            break;
        case  BSIM4_MOD_LCF :
            mod.BSIM4lcf = value.get_as<double>().value();
            mod.BSIM4lcfGiven = true;
            break;
        case  BSIM4_MOD_LCLC :
            mod.BSIM4lclc = value.get_as<double>().value();
            mod.BSIM4lclcGiven = true;
            break;
        case  BSIM4_MOD_LCLE :
            mod.BSIM4lcle = value.get_as<double>().value();
            mod.BSIM4lcleGiven = true;
            break;
        case  BSIM4_MOD_LVFBCV :
            mod.BSIM4lvfbcv = value.get_as<double>().value();
            mod.BSIM4lvfbcvGiven = true;
            break;
        case  BSIM4_MOD_LACDE :
            mod.BSIM4lacde = value.get_as<double>().value();
            mod.BSIM4lacdeGiven = true;
            break;
        case  BSIM4_MOD_LMOIN :
            mod.BSIM4lmoin = value.get_as<double>().value();
            mod.BSIM4lmoinGiven = true;
            break;
        case  BSIM4_MOD_LNOFF :
            mod.BSIM4lnoff = value.get_as<double>().value();
            mod.BSIM4lnoffGiven = true;
            break;
        case  BSIM4_MOD_LVOFFCV :
            mod.BSIM4lvoffcv = value.get_as<double>().value();
            mod.BSIM4lvoffcvGiven = true;
            break;

    /* Width dependence */
        case  BSIM4_MOD_WCDSC :
            mod.BSIM4wcdsc = value.get_as<double>().value();
            mod.BSIM4wcdscGiven = true;
            break;


         case  BSIM4_MOD_WCDSCB :
            mod.BSIM4wcdscb = value.get_as<double>().value();
            mod.BSIM4wcdscbGiven = true;
            break;
         case  BSIM4_MOD_WCDSCD :
            mod.BSIM4wcdscd = value.get_as<double>().value();
            mod.BSIM4wcdscdGiven = true;
            break;
        case  BSIM4_MOD_WCIT :
            mod.BSIM4wcit = value.get_as<double>().value();
            mod.BSIM4wcitGiven = true;
            break;
        case  BSIM4_MOD_WNFACTOR :
            mod.BSIM4wnfactor = value.get_as<double>().value();
            mod.BSIM4wnfactorGiven = true;
            break;
        case BSIM4_MOD_WXJ:
            mod.BSIM4wxj = value.get_as<double>().value();
            mod.BSIM4wxjGiven = true;
            break;
        case BSIM4_MOD_WVSAT:
            mod.BSIM4wvsat = value.get_as<double>().value();
            mod.BSIM4wvsatGiven = true;
            break;


        case BSIM4_MOD_WA0:
            mod.BSIM4wa0 = value.get_as<double>().value();
            mod.BSIM4wa0Given = true;
            break;
        case BSIM4_MOD_WAGS:
            mod.BSIM4wags = value.get_as<double>().value();
            mod.BSIM4wagsGiven = true;
            break;
        case BSIM4_MOD_WA1:
            mod.BSIM4wa1 = value.get_as<double>().value();
            mod.BSIM4wa1Given = true;
            break;
        case BSIM4_MOD_WA2:
            mod.BSIM4wa2 = value.get_as<double>().value();
            mod.BSIM4wa2Given = true;
            break;
        case BSIM4_MOD_WAT:
            mod.BSIM4wat = value.get_as<double>().value();
            mod.BSIM4watGiven = true;
            break;
        case BSIM4_MOD_WKETA:
            mod.BSIM4wketa = value.get_as<double>().value();
            mod.BSIM4wketaGiven = true;
            break;
        case BSIM4_MOD_WNSUB:
            mod.BSIM4wnsub = value.get_as<double>().value();
            mod.BSIM4wnsubGiven = true;
            break;
        case BSIM4_MOD_WNDEP:
            mod.BSIM4wndep = value.get_as<double>().value();
            mod.BSIM4wndepGiven = true;
        if (mod.BSIM4wndep > 1.0e20)
        mod.BSIM4wndep *= 1.0e-6;
            break;
        case BSIM4_MOD_WNSD:
            mod.BSIM4wnsd = value.get_as<double>().value();
            mod.BSIM4wnsdGiven = true;
            if (mod.BSIM4wnsd > 1.0e23)
                mod.BSIM4wnsd *= 1.0e-6;
            break;
        case BSIM4_MOD_WNGATE:
            mod.BSIM4wngate = value.get_as<double>().value();
            mod.BSIM4wngateGiven = true;
        if (mod.BSIM4wngate > 1.0e23)
        mod.BSIM4wngate *= 1.0e-6;
            break;
        case BSIM4_MOD_WGAMMA1:
            mod.BSIM4wgamma1 = value.get_as<double>().value();
            mod.BSIM4wgamma1Given = true;
            break;
        case BSIM4_MOD_WGAMMA2:
            mod.BSIM4wgamma2 = value.get_as<double>().value();
            mod.BSIM4wgamma2Given = true;
            break;
        case BSIM4_MOD_WVBX:
            mod.BSIM4wvbx = value.get_as<double>().value();
            mod.BSIM4wvbxGiven = true;
            break;
        case BSIM4_MOD_WVBM:
            mod.BSIM4wvbm = value.get_as<double>().value();
            mod.BSIM4wvbmGiven = true;
            break;
        case BSIM4_MOD_WXT:
            mod.BSIM4wxt = value.get_as<double>().value();
            mod.BSIM4wxtGiven = true;
            break;
        case  BSIM4_MOD_WK1:
            mod.BSIM4wk1 = value.get_as<double>().value();
            mod.BSIM4wk1Given = true;
            break;
        case  BSIM4_MOD_WKT1:
            mod.BSIM4wkt1 = value.get_as<double>().value();
            mod.BSIM4wkt1Given = true;
            break;
        case  BSIM4_MOD_WKT1L:
            mod.BSIM4wkt1l = value.get_as<double>().value();
            mod.BSIM4wkt1lGiven = true;
            break;
        case  BSIM4_MOD_WKT2:
            mod.BSIM4wkt2 = value.get_as<double>().value();
            mod.BSIM4wkt2Given = true;
            break;
        case  BSIM4_MOD_WK2:
            mod.BSIM4wk2 = value.get_as<double>().value();
            mod.BSIM4wk2Given = true;
            break;
        case  BSIM4_MOD_WK3:
            mod.BSIM4wk3 = value.get_as<double>().value();
            mod.BSIM4wk3Given = true;
            break;
        case  BSIM4_MOD_WK3B:
            mod.BSIM4wk3b = value.get_as<double>().value();
            mod.BSIM4wk3bGiven = true;
            break;
        case  BSIM4_MOD_WLPE0:
            mod.BSIM4wlpe0 = value.get_as<double>().value();
            mod.BSIM4wlpe0Given = true;
            break;
        case  BSIM4_MOD_WLPEB:
            mod.BSIM4wlpeb = value.get_as<double>().value();
            mod.BSIM4wlpebGiven = true;
            break;
        case  BSIM4_MOD_WDVTP0:
            mod.BSIM4wdvtp0 = value.get_as<double>().value();
            mod.BSIM4wdvtp0Given = true;
            break;
        case  BSIM4_MOD_WDVTP1:
            mod.BSIM4wdvtp1 = value.get_as<double>().value();
            mod.BSIM4wdvtp1Given = true;
            break;
    case  BSIM4_MOD_WDVTP2:     /* New DIBL/Rout */
            mod.BSIM4wdvtp2 = value.get_as<double>().value();
            mod.BSIM4wdvtp2Given = true;
            break;
    case  BSIM4_MOD_WDVTP3:
            mod.BSIM4wdvtp3 = value.get_as<double>().value();
            mod.BSIM4wdvtp3Given = true;
            break;
    case  BSIM4_MOD_WDVTP4:
            mod.BSIM4wdvtp4 = value.get_as<double>().value();
            mod.BSIM4wdvtp4Given = true;
            break;
    case  BSIM4_MOD_WDVTP5:
            mod.BSIM4wdvtp5 = value.get_as<double>().value();
            mod.BSIM4wdvtp5Given = true;
            break;
        case  BSIM4_MOD_WW0:
            mod.BSIM4ww0 = value.get_as<double>().value();
            mod.BSIM4ww0Given = true;
            break;
        case  BSIM4_MOD_WDVT0:
            mod.BSIM4wdvt0 = value.get_as<double>().value();
            mod.BSIM4wdvt0Given = true;
            break;
        case  BSIM4_MOD_WDVT1:
            mod.BSIM4wdvt1 = value.get_as<double>().value();
            mod.BSIM4wdvt1Given = true;
            break;
        case  BSIM4_MOD_WDVT2:
            mod.BSIM4wdvt2 = value.get_as<double>().value();
            mod.BSIM4wdvt2Given = true;
            break;
        case  BSIM4_MOD_WDVT0W:
            mod.BSIM4wdvt0w = value.get_as<double>().value();
            mod.BSIM4wdvt0wGiven = true;
            break;
        case  BSIM4_MOD_WDVT1W:
            mod.BSIM4wdvt1w = value.get_as<double>().value();
            mod.BSIM4wdvt1wGiven = true;
            break;
        case  BSIM4_MOD_WDVT2W:
            mod.BSIM4wdvt2w = value.get_as<double>().value();
            mod.BSIM4wdvt2wGiven = true;
            break;
        case  BSIM4_MOD_WDROUT:
            mod.BSIM4wdrout = value.get_as<double>().value();
            mod.BSIM4wdroutGiven = true;
            break;
        case  BSIM4_MOD_WDSUB:
            mod.BSIM4wdsub = value.get_as<double>().value();
            mod.BSIM4wdsubGiven = true;
            break;
        case BSIM4_MOD_WVTH0:
            mod.BSIM4wvth0 = value.get_as<double>().value();
            mod.BSIM4wvth0Given = true;
            break;
        case BSIM4_MOD_WUA:
            mod.BSIM4wua = value.get_as<double>().value();
            mod.BSIM4wuaGiven = true;
            break;
        case BSIM4_MOD_WUA1:
            mod.BSIM4wua1 = value.get_as<double>().value();
            mod.BSIM4wua1Given = true;
            break;
        case BSIM4_MOD_WUB:
            mod.BSIM4wub = value.get_as<double>().value();
            mod.BSIM4wubGiven = true;
            break;
        case BSIM4_MOD_WUB1:
            mod.BSIM4wub1 = value.get_as<double>().value();
            mod.BSIM4wub1Given = true;
            break;
        case BSIM4_MOD_WUC:
            mod.BSIM4wuc = value.get_as<double>().value();
            mod.BSIM4wucGiven = true;
            break;
        case BSIM4_MOD_WUC1:
            mod.BSIM4wuc1 = value.get_as<double>().value();
            mod.BSIM4wuc1Given = true;
            break;
        case  BSIM4_MOD_WU0 :
            mod.BSIM4wu0 = value.get_as<double>().value();
            mod.BSIM4wu0Given = true;
            break;
        case  BSIM4_MOD_WUTE :
            mod.BSIM4wute = value.get_as<double>().value();
            mod.BSIM4wuteGiven = true;
            break;
        case  BSIM4_MOD_WUCSTE :
            mod.BSIM4wucste = value.get_as<double>().value();
            mod.BSIM4wucsteGiven = true;
            break;
        case BSIM4_MOD_WVOFF:
            mod.BSIM4wvoff = value.get_as<double>().value();
            mod.BSIM4wvoffGiven = true;
            break;
        case BSIM4_MOD_WTVOFF:
            mod.BSIM4wtvoff = value.get_as<double>().value();
            mod.BSIM4wtvoffGiven = true;
            break;
        case BSIM4_MOD_WTNFACTOR:       /* v4.7 temp dep of leakage current  */
            mod.BSIM4wtnfactor = value.get_as<double>().value();
            mod.BSIM4wtnfactorGiven = true;
            break;
        case BSIM4_MOD_WTETA0:      /* v4.7 temp dep of leakage current  */
            mod.BSIM4wteta0 = value.get_as<double>().value();
            mod.BSIM4wteta0Given = true;
            break;
        case BSIM4_MOD_WTVOFFCV:    /* v4.7 temp dep of leakage current  */
            mod.BSIM4wtvoffcv = value.get_as<double>().value();
            mod.BSIM4wtvoffcvGiven = true;
            break;
        case BSIM4_MOD_WMINV:
            mod.BSIM4wminv = value.get_as<double>().value();
            mod.BSIM4wminvGiven = true;
            break;
        case BSIM4_MOD_WMINVCV:
            mod.BSIM4wminvcv = value.get_as<double>().value();
            mod.BSIM4wminvcvGiven = true;
            break;
        case BSIM4_MOD_WFPROUT:
            mod.BSIM4wfprout = value.get_as<double>().value();
            mod.BSIM4wfproutGiven = true;
            break;
        case BSIM4_MOD_WPDITS:
            mod.BSIM4wpdits = value.get_as<double>().value();
            mod.BSIM4wpditsGiven = true;
            break;
        case BSIM4_MOD_WPDITSD:
            mod.BSIM4wpditsd = value.get_as<double>().value();
            mod.BSIM4wpditsdGiven = true;
            break;
        case  BSIM4_MOD_WDELTA :
            mod.BSIM4wdelta = value.get_as<double>().value();
            mod.BSIM4wdeltaGiven = true;
            break;
        case BSIM4_MOD_WRDSW:
            mod.BSIM4wrdsw = value.get_as<double>().value();
            mod.BSIM4wrdswGiven = true;
            break;
        case BSIM4_MOD_WRDW:
            mod.BSIM4wrdw = value.get_as<double>().value();
            mod.BSIM4wrdwGiven = true;
            break;
        case BSIM4_MOD_WRSW:
            mod.BSIM4wrsw = value.get_as<double>().value();
            mod.BSIM4wrswGiven = true;
            break;
        case BSIM4_MOD_WPRWB:
            mod.BSIM4wprwb = value.get_as<double>().value();
            mod.BSIM4wprwbGiven = true;
            break;
        case BSIM4_MOD_WPRWG:
            mod.BSIM4wprwg = value.get_as<double>().value();
            mod.BSIM4wprwgGiven = true;
            break;
        case BSIM4_MOD_WPRT:
            mod.BSIM4wprt = value.get_as<double>().value();
            mod.BSIM4wprtGiven = true;
            break;
        case BSIM4_MOD_WETA0:
            mod.BSIM4weta0 = value.get_as<double>().value();
            mod.BSIM4weta0Given = true;
            break;
        case BSIM4_MOD_WETAB:
            mod.BSIM4wetab = value.get_as<double>().value();
            mod.BSIM4wetabGiven = true;
            break;
        case BSIM4_MOD_WPCLM:
            mod.BSIM4wpclm = value.get_as<double>().value();
            mod.BSIM4wpclmGiven = true;
            break;
        case BSIM4_MOD_WPDIBL1:
            mod.BSIM4wpdibl1 = value.get_as<double>().value();
            mod.BSIM4wpdibl1Given = true;
            break;
        case BSIM4_MOD_WPDIBL2:
            mod.BSIM4wpdibl2 = value.get_as<double>().value();
            mod.BSIM4wpdibl2Given = true;
            break;
        case BSIM4_MOD_WPDIBLB:
            mod.BSIM4wpdiblb = value.get_as<double>().value();
            mod.BSIM4wpdiblbGiven = true;
            break;
        case BSIM4_MOD_WPSCBE1:
            mod.BSIM4wpscbe1 = value.get_as<double>().value();
            mod.BSIM4wpscbe1Given = true;
            break;
        case BSIM4_MOD_WPSCBE2:
            mod.BSIM4wpscbe2 = value.get_as<double>().value();
            mod.BSIM4wpscbe2Given = true;
            break;
        case BSIM4_MOD_WPVAG:
            mod.BSIM4wpvag = value.get_as<double>().value();
            mod.BSIM4wpvagGiven = true;
            break;
        case  BSIM4_MOD_WWR :
            mod.BSIM4wwr = value.get_as<double>().value();
            mod.BSIM4wwrGiven = true;
            break;
        case  BSIM4_MOD_WDWG :
            mod.BSIM4wdwg = value.get_as<double>().value();
            mod.BSIM4wdwgGiven = true;
            break;
        case  BSIM4_MOD_WDWB :
            mod.BSIM4wdwb = value.get_as<double>().value();
            mod.BSIM4wdwbGiven = true;
            break;
        case  BSIM4_MOD_WB0 :
            mod.BSIM4wb0 = value.get_as<double>().value();
            mod.BSIM4wb0Given = true;
            break;
        case  BSIM4_MOD_WB1 :
            mod.BSIM4wb1 = value.get_as<double>().value();
            mod.BSIM4wb1Given = true;
            break;
        case  BSIM4_MOD_WALPHA0 :
            mod.BSIM4walpha0 = value.get_as<double>().value();
            mod.BSIM4walpha0Given = true;
            break;
        case  BSIM4_MOD_WALPHA1 :
            mod.BSIM4walpha1 = value.get_as<double>().value();
            mod.BSIM4walpha1Given = true;
            break;
        case  BSIM4_MOD_WBETA0 :
            mod.BSIM4wbeta0 = value.get_as<double>().value();
            mod.BSIM4wbeta0Given = true;
            break;
        case  BSIM4_MOD_WPHIN :
            mod.BSIM4wphin = value.get_as<double>().value();
            mod.BSIM4wphinGiven = true;
            break;
        case  BSIM4_MOD_WAGIDL :
            mod.BSIM4wagidl = value.get_as<double>().value();
            mod.BSIM4wagidlGiven = true;
            break;
        case  BSIM4_MOD_WBGIDL :
            mod.BSIM4wbgidl = value.get_as<double>().value();
            mod.BSIM4wbgidlGiven = true;
            break;
        case  BSIM4_MOD_WCGIDL :
            mod.BSIM4wcgidl = value.get_as<double>().value();
            mod.BSIM4wcgidlGiven = true;
            break;
        case  BSIM4_MOD_WEGIDL :
            mod.BSIM4wegidl = value.get_as<double>().value();
            mod.BSIM4wegidlGiven = true;
            break;
    case  BSIM4_MOD_WFGIDL :            /* v4.7 New GIDL/GISL */
            mod.BSIM4wfgidl = value.get_as<double>().value();
            mod.BSIM4wfgidlGiven = true;
            break;
    case  BSIM4_MOD_WKGIDL :            /* v4.7 New GIDL/GISL */
            mod.BSIM4wkgidl = value.get_as<double>().value();
            mod.BSIM4wkgidlGiven = true;
            break;
    case  BSIM4_MOD_WRGIDL :            /* v4.7 New GIDL/GISL */
            mod.BSIM4wrgidl = value.get_as<double>().value();
            mod.BSIM4wrgidlGiven = true;
            break;
        case  BSIM4_MOD_WAGISL :
            mod.BSIM4wagisl = value.get_as<double>().value();
            mod.BSIM4wagislGiven = true;
            break;
        case  BSIM4_MOD_WBGISL :
            mod.BSIM4wbgisl = value.get_as<double>().value();
            mod.BSIM4wbgislGiven = true;
            break;
        case  BSIM4_MOD_WCGISL :
            mod.BSIM4wcgisl = value.get_as<double>().value();
            mod.BSIM4wcgislGiven = true;
            break;
        case  BSIM4_MOD_WEGISL :
            mod.BSIM4wegisl = value.get_as<double>().value();
            mod.BSIM4wegislGiven = true;
            break;
    case  BSIM4_MOD_WFGISL :            /* v4.7 New GIDL/GISL */
            mod.BSIM4wfgisl = value.get_as<double>().value();
            mod.BSIM4wfgislGiven = true;
            break;
    case  BSIM4_MOD_WKGISL :            /* v4.7 New GIDL/GISL */
            mod.BSIM4wkgisl = value.get_as<double>().value();
            mod.BSIM4wkgislGiven = true;
            break;
    case  BSIM4_MOD_WRGISL :            /* v4.7 New GIDL/GISL */
            mod.BSIM4wrgisl = value.get_as<double>().value();
            mod.BSIM4wrgislGiven = true;
            break;
        case  BSIM4_MOD_WAIGC :
            mod.BSIM4waigc = value.get_as<double>().value();
            mod.BSIM4waigcGiven = true;
            break;
        case  BSIM4_MOD_WBIGC :
            mod.BSIM4wbigc = value.get_as<double>().value();
            mod.BSIM4wbigcGiven = true;
            break;
        case  BSIM4_MOD_WCIGC :
            mod.BSIM4wcigc = value.get_as<double>().value();
            mod.BSIM4wcigcGiven = true;
            break;
        case  BSIM4_MOD_WAIGSD :
            mod.BSIM4waigsd = value.get_as<double>().value();
            mod.BSIM4waigsdGiven = true;
            break;
        case  BSIM4_MOD_WBIGSD :
            mod.BSIM4wbigsd = value.get_as<double>().value();
            mod.BSIM4wbigsdGiven = true;
            break;
        case  BSIM4_MOD_WCIGSD :
            mod.BSIM4wcigsd = value.get_as<double>().value();
            mod.BSIM4wcigsdGiven = true;
            break;
        case  BSIM4_MOD_WAIGS :
            mod.BSIM4waigs = value.get_as<double>().value();
            mod.BSIM4waigsGiven = true;
            break;
        case  BSIM4_MOD_WBIGS :
            mod.BSIM4wbigs = value.get_as<double>().value();
            mod.BSIM4wbigsGiven = true;
            break;
        case  BSIM4_MOD_WCIGS :
            mod.BSIM4wcigs = value.get_as<double>().value();
            mod.BSIM4wcigsGiven = true;
            break;
        case  BSIM4_MOD_WAIGD :
            mod.BSIM4waigd = value.get_as<double>().value();
            mod.BSIM4waigdGiven = true;
            break;
        case  BSIM4_MOD_WBIGD :
            mod.BSIM4wbigd = value.get_as<double>().value();
            mod.BSIM4wbigdGiven = true;
            break;
        case  BSIM4_MOD_WCIGD :
            mod.BSIM4wcigd = value.get_as<double>().value();
            mod.BSIM4wcigdGiven = true;
            break;
        case  BSIM4_MOD_WAIGBACC :
            mod.BSIM4waigbacc = value.get_as<double>().value();
            mod.BSIM4waigbaccGiven = true;
            break;
        case  BSIM4_MOD_WBIGBACC :
            mod.BSIM4wbigbacc = value.get_as<double>().value();
            mod.BSIM4wbigbaccGiven = true;
            break;
        case  BSIM4_MOD_WCIGBACC :
            mod.BSIM4wcigbacc = value.get_as<double>().value();
            mod.BSIM4wcigbaccGiven = true;
            break;
        case  BSIM4_MOD_WAIGBINV :
            mod.BSIM4waigbinv = value.get_as<double>().value();
            mod.BSIM4waigbinvGiven = true;
            break;
        case  BSIM4_MOD_WBIGBINV :
            mod.BSIM4wbigbinv = value.get_as<double>().value();
            mod.BSIM4wbigbinvGiven = true;
            break;
        case  BSIM4_MOD_WCIGBINV :
            mod.BSIM4wcigbinv = value.get_as<double>().value();
            mod.BSIM4wcigbinvGiven = true;
            break;
        case  BSIM4_MOD_WNIGC :
            mod.BSIM4wnigc = value.get_as<double>().value();
            mod.BSIM4wnigcGiven = true;
            break;
        case  BSIM4_MOD_WNIGBINV :
            mod.BSIM4wnigbinv = value.get_as<double>().value();
            mod.BSIM4wnigbinvGiven = true;
            break;
        case  BSIM4_MOD_WNIGBACC :
            mod.BSIM4wnigbacc = value.get_as<double>().value();
            mod.BSIM4wnigbaccGiven = true;
            break;
        case  BSIM4_MOD_WNTOX :
            mod.BSIM4wntox = value.get_as<double>().value();
            mod.BSIM4wntoxGiven = true;
            break;
        case  BSIM4_MOD_WEIGBINV :
            mod.BSIM4weigbinv = value.get_as<double>().value();
            mod.BSIM4weigbinvGiven = true;
            break;
        case  BSIM4_MOD_WPIGCD :
            mod.BSIM4wpigcd = value.get_as<double>().value();
            mod.BSIM4wpigcdGiven = true;
            break;
        case  BSIM4_MOD_WPOXEDGE :
            mod.BSIM4wpoxedge = value.get_as<double>().value();
            mod.BSIM4wpoxedgeGiven = true;
            break;
        case  BSIM4_MOD_WXRCRG1 :
            mod.BSIM4wxrcrg1 = value.get_as<double>().value();
            mod.BSIM4wxrcrg1Given = true;
            break;
        case  BSIM4_MOD_WXRCRG2 :
            mod.BSIM4wxrcrg2 = value.get_as<double>().value();
            mod.BSIM4wxrcrg2Given = true;
            break;
        case  BSIM4_MOD_WLAMBDA :
            mod.BSIM4wlambda = value.get_as<double>().value();
            mod.BSIM4wlambdaGiven = true;
            break;
        case  BSIM4_MOD_WVTL :
            mod.BSIM4wvtl = value.get_as<double>().value();
            mod.BSIM4wvtlGiven = true;
            break;
        case  BSIM4_MOD_WXN:
            mod.BSIM4wxn = value.get_as<double>().value();
            mod.BSIM4wxnGiven = true;
            break;
        case  BSIM4_MOD_WVFBSDOFF:
            mod.BSIM4wvfbsdoff = value.get_as<double>().value();
            mod.BSIM4wvfbsdoffGiven = true;
            break;
        case  BSIM4_MOD_WTVFBSDOFF:
            mod.BSIM4wtvfbsdoff = value.get_as<double>().value();
            mod.BSIM4wtvfbsdoffGiven = true;
            break;
        case  BSIM4_MOD_WEU :
            mod.BSIM4weu = value.get_as<double>().value();
            mod.BSIM4weuGiven = true;
            break;
         case  BSIM4_MOD_WUCS :
            mod.BSIM4wucs = value.get_as<double>().value();
            mod.BSIM4wucsGiven = true;
            break;
        case  BSIM4_MOD_WVFB :
            mod.BSIM4wvfb = value.get_as<double>().value();
            mod.BSIM4wvfbGiven = true;
            break;
        case  BSIM4_MOD_WCGSL :
            mod.BSIM4wcgsl = value.get_as<double>().value();
            mod.BSIM4wcgslGiven = true;
            break;
        case  BSIM4_MOD_WCGDL :
            mod.BSIM4wcgdl = value.get_as<double>().value();
            mod.BSIM4wcgdlGiven = true;
            break;
        case  BSIM4_MOD_WCKAPPAS :
            mod.BSIM4wckappas = value.get_as<double>().value();
            mod.BSIM4wckappasGiven = true;
            break;
        case  BSIM4_MOD_WCKAPPAD :
            mod.BSIM4wckappad = value.get_as<double>().value();
            mod.BSIM4wckappadGiven = true;
            break;
        case  BSIM4_MOD_WCF :
            mod.BSIM4wcf = value.get_as<double>().value();
            mod.BSIM4wcfGiven = true;
            break;
        case  BSIM4_MOD_WCLC :
            mod.BSIM4wclc = value.get_as<double>().value();
            mod.BSIM4wclcGiven = true;
            break;
        case  BSIM4_MOD_WCLE :
            mod.BSIM4wcle = value.get_as<double>().value();
            mod.BSIM4wcleGiven = true;
            break;
        case  BSIM4_MOD_WVFBCV :
            mod.BSIM4wvfbcv = value.get_as<double>().value();
            mod.BSIM4wvfbcvGiven = true;
            break;
        case  BSIM4_MOD_WACDE :
            mod.BSIM4wacde = value.get_as<double>().value();
            mod.BSIM4wacdeGiven = true;
            break;
        case  BSIM4_MOD_WMOIN :
            mod.BSIM4wmoin = value.get_as<double>().value();
            mod.BSIM4wmoinGiven = true;
            break;
        case  BSIM4_MOD_WNOFF :
            mod.BSIM4wnoff = value.get_as<double>().value();
            mod.BSIM4wnoffGiven = true;
            break;
        case  BSIM4_MOD_WVOFFCV :
            mod.BSIM4wvoffcv = value.get_as<double>().value();
            mod.BSIM4wvoffcvGiven = true;
            break;

    /* Cross-term dependence */
        case  BSIM4_MOD_PCDSC :
            mod.BSIM4pcdsc = value.get_as<double>().value();
            mod.BSIM4pcdscGiven = true;
            break;


        case  BSIM4_MOD_PCDSCB :
            mod.BSIM4pcdscb = value.get_as<double>().value();
            mod.BSIM4pcdscbGiven = true;
            break;
        case  BSIM4_MOD_PCDSCD :
            mod.BSIM4pcdscd = value.get_as<double>().value();
            mod.BSIM4pcdscdGiven = true;
            break;
        case  BSIM4_MOD_PCIT :
            mod.BSIM4pcit = value.get_as<double>().value();
            mod.BSIM4pcitGiven = true;
            break;
        case  BSIM4_MOD_PNFACTOR :
            mod.BSIM4pnfactor = value.get_as<double>().value();
            mod.BSIM4pnfactorGiven = true;
            break;
        case BSIM4_MOD_PXJ:
            mod.BSIM4pxj = value.get_as<double>().value();
            mod.BSIM4pxjGiven = true;
            break;
        case BSIM4_MOD_PVSAT:
            mod.BSIM4pvsat = value.get_as<double>().value();
            mod.BSIM4pvsatGiven = true;
            break;


        case BSIM4_MOD_PA0:
            mod.BSIM4pa0 = value.get_as<double>().value();
            mod.BSIM4pa0Given = true;
            break;
        case BSIM4_MOD_PAGS:
            mod.BSIM4pags = value.get_as<double>().value();
            mod.BSIM4pagsGiven = true;
            break;
        case BSIM4_MOD_PA1:
            mod.BSIM4pa1 = value.get_as<double>().value();
            mod.BSIM4pa1Given = true;
            break;
        case BSIM4_MOD_PA2:
            mod.BSIM4pa2 = value.get_as<double>().value();
            mod.BSIM4pa2Given = true;
            break;
        case BSIM4_MOD_PAT:
            mod.BSIM4pat = value.get_as<double>().value();
            mod.BSIM4patGiven = true;
            break;
        case BSIM4_MOD_PKETA:
            mod.BSIM4pketa = value.get_as<double>().value();
            mod.BSIM4pketaGiven = true;
            break;
        case BSIM4_MOD_PNSUB:
            mod.BSIM4pnsub = value.get_as<double>().value();
            mod.BSIM4pnsubGiven = true;
            break;
        case BSIM4_MOD_PNDEP:
            mod.BSIM4pndep = value.get_as<double>().value();
            mod.BSIM4pndepGiven = true;
        if (mod.BSIM4pndep > 1.0e20)
        mod.BSIM4pndep *= 1.0e-6;
            break;
        case BSIM4_MOD_PNSD:
            mod.BSIM4pnsd = value.get_as<double>().value();
            mod.BSIM4pnsdGiven = true;
            if (mod.BSIM4pnsd > 1.0e23)
                mod.BSIM4pnsd *= 1.0e-6;
            break;
        case BSIM4_MOD_PNGATE:
            mod.BSIM4pngate = value.get_as<double>().value();
            mod.BSIM4pngateGiven = true;
        if (mod.BSIM4pngate > 1.0e23)
        mod.BSIM4pngate *= 1.0e-6;
            break;
        case BSIM4_MOD_PGAMMA1:
            mod.BSIM4pgamma1 = value.get_as<double>().value();
            mod.BSIM4pgamma1Given = true;
            break;
        case BSIM4_MOD_PGAMMA2:
            mod.BSIM4pgamma2 = value.get_as<double>().value();
            mod.BSIM4pgamma2Given = true;
            break;
        case BSIM4_MOD_PVBX:
            mod.BSIM4pvbx = value.get_as<double>().value();
            mod.BSIM4pvbxGiven = true;
            break;
        case BSIM4_MOD_PVBM:
            mod.BSIM4pvbm = value.get_as<double>().value();
            mod.BSIM4pvbmGiven = true;
            break;
        case BSIM4_MOD_PXT:
            mod.BSIM4pxt = value.get_as<double>().value();
            mod.BSIM4pxtGiven = true;
            break;
        case  BSIM4_MOD_PK1:
            mod.BSIM4pk1 = value.get_as<double>().value();
            mod.BSIM4pk1Given = true;
            break;
        case  BSIM4_MOD_PKT1:
            mod.BSIM4pkt1 = value.get_as<double>().value();
            mod.BSIM4pkt1Given = true;
            break;
        case  BSIM4_MOD_PKT1L:
            mod.BSIM4pkt1l = value.get_as<double>().value();
            mod.BSIM4pkt1lGiven = true;
            break;
        case  BSIM4_MOD_PKT2:
            mod.BSIM4pkt2 = value.get_as<double>().value();
            mod.BSIM4pkt2Given = true;
            break;
        case  BSIM4_MOD_PK2:
            mod.BSIM4pk2 = value.get_as<double>().value();
            mod.BSIM4pk2Given = true;
            break;
        case  BSIM4_MOD_PK3:
            mod.BSIM4pk3 = value.get_as<double>().value();
            mod.BSIM4pk3Given = true;
            break;
        case  BSIM4_MOD_PK3B:
            mod.BSIM4pk3b = value.get_as<double>().value();
            mod.BSIM4pk3bGiven = true;
            break;
        case  BSIM4_MOD_PLPE0:
            mod.BSIM4plpe0 = value.get_as<double>().value();
            mod.BSIM4plpe0Given = true;
            break;
        case  BSIM4_MOD_PLPEB:
            mod.BSIM4plpeb = value.get_as<double>().value();
            mod.BSIM4plpebGiven = true;
            break;
        case  BSIM4_MOD_PDVTP0:
            mod.BSIM4pdvtp0 = value.get_as<double>().value();
            mod.BSIM4pdvtp0Given = true;
            break;
        case  BSIM4_MOD_PDVTP1:
            mod.BSIM4pdvtp1 = value.get_as<double>().value();
            mod.BSIM4pdvtp1Given = true;
            break;
    case  BSIM4_MOD_PDVTP2:     /* New DIBL/Rout */
            mod.BSIM4pdvtp2 = value.get_as<double>().value();
            mod.BSIM4pdvtp2Given = true;
            break;
    case  BSIM4_MOD_PDVTP3:
            mod.BSIM4pdvtp3 = value.get_as<double>().value();
            mod.BSIM4pdvtp3Given = true;
            break;
    case  BSIM4_MOD_PDVTP4:
            mod.BSIM4pdvtp4 = value.get_as<double>().value();
            mod.BSIM4pdvtp4Given = true;
            break;
    case  BSIM4_MOD_PDVTP5:
            mod.BSIM4pdvtp5 = value.get_as<double>().value();
            mod.BSIM4pdvtp5Given = true;
            break;
        case  BSIM4_MOD_PW0:
            mod.BSIM4pw0 = value.get_as<double>().value();
            mod.BSIM4pw0Given = true;
            break;
        case  BSIM4_MOD_PDVT0:
            mod.BSIM4pdvt0 = value.get_as<double>().value();
            mod.BSIM4pdvt0Given = true;
            break;
        case  BSIM4_MOD_PDVT1:
            mod.BSIM4pdvt1 = value.get_as<double>().value();
            mod.BSIM4pdvt1Given = true;
            break;
        case  BSIM4_MOD_PDVT2:
            mod.BSIM4pdvt2 = value.get_as<double>().value();
            mod.BSIM4pdvt2Given = true;
            break;
        case  BSIM4_MOD_PDVT0W:
            mod.BSIM4pdvt0w = value.get_as<double>().value();
            mod.BSIM4pdvt0wGiven = true;
            break;
        case  BSIM4_MOD_PDVT1W:
            mod.BSIM4pdvt1w = value.get_as<double>().value();
            mod.BSIM4pdvt1wGiven = true;
            break;
        case  BSIM4_MOD_PDVT2W:
            mod.BSIM4pdvt2w = value.get_as<double>().value();
            mod.BSIM4pdvt2wGiven = true;
            break;
        case  BSIM4_MOD_PDROUT:
            mod.BSIM4pdrout = value.get_as<double>().value();
            mod.BSIM4pdroutGiven = true;
            break;
        case  BSIM4_MOD_PDSUB:
            mod.BSIM4pdsub = value.get_as<double>().value();
            mod.BSIM4pdsubGiven = true;
            break;
        case BSIM4_MOD_PVTH0:
            mod.BSIM4pvth0 = value.get_as<double>().value();
            mod.BSIM4pvth0Given = true;
            break;
        case BSIM4_MOD_PUA:
            mod.BSIM4pua = value.get_as<double>().value();
            mod.BSIM4puaGiven = true;
            break;
        case BSIM4_MOD_PUA1:
            mod.BSIM4pua1 = value.get_as<double>().value();
            mod.BSIM4pua1Given = true;
            break;
        case BSIM4_MOD_PUB:
            mod.BSIM4pub = value.get_as<double>().value();
            mod.BSIM4pubGiven = true;
            break;
        case BSIM4_MOD_PUB1:
            mod.BSIM4pub1 = value.get_as<double>().value();
            mod.BSIM4pub1Given = true;
            break;
        case BSIM4_MOD_PUC:
            mod.BSIM4puc = value.get_as<double>().value();
            mod.BSIM4pucGiven = true;
            break;
        case BSIM4_MOD_PUC1:
            mod.BSIM4puc1 = value.get_as<double>().value();
            mod.BSIM4puc1Given = true;
            break;
        case  BSIM4_MOD_PU0 :
            mod.BSIM4pu0 = value.get_as<double>().value();
            mod.BSIM4pu0Given = true;
            break;
        case  BSIM4_MOD_PUTE :
            mod.BSIM4pute = value.get_as<double>().value();
            mod.BSIM4puteGiven = true;
            break;
         case  BSIM4_MOD_PUCSTE :
            mod.BSIM4pucste = value.get_as<double>().value();
            mod.BSIM4pucsteGiven = true;
            break;
        case BSIM4_MOD_PVOFF:
            mod.BSIM4pvoff = value.get_as<double>().value();
            mod.BSIM4pvoffGiven = true;
            break;
        case BSIM4_MOD_PTVOFF:
            mod.BSIM4ptvoff = value.get_as<double>().value();
            mod.BSIM4ptvoffGiven = true;
            break;
        case BSIM4_MOD_PTNFACTOR:       /* v4.7 temp dep of leakage current  */
            mod.BSIM4ptnfactor = value.get_as<double>().value();
            mod.BSIM4ptnfactorGiven = true;
            break;
        case BSIM4_MOD_PTETA0:      /* v4.7 temp dep of leakage current  */
            mod.BSIM4pteta0 = value.get_as<double>().value();
            mod.BSIM4pteta0Given = true;
            break;
        case BSIM4_MOD_PTVOFFCV:    /* v4.7 temp dep of leakage current  */
            mod.BSIM4ptvoffcv = value.get_as<double>().value();
            mod.BSIM4ptvoffcvGiven = true;
            break;
        case BSIM4_MOD_PMINV:
            mod.BSIM4pminv = value.get_as<double>().value();
            mod.BSIM4pminvGiven = true;
            break;
        case BSIM4_MOD_PMINVCV:
            mod.BSIM4pminvcv = value.get_as<double>().value();
            mod.BSIM4pminvcvGiven = true;
            break;
        case BSIM4_MOD_PFPROUT:
            mod.BSIM4pfprout = value.get_as<double>().value();
            mod.BSIM4pfproutGiven = true;
            break;
        case BSIM4_MOD_PPDITS:
            mod.BSIM4ppdits = value.get_as<double>().value();
            mod.BSIM4ppditsGiven = true;
            break;
        case BSIM4_MOD_PPDITSD:
            mod.BSIM4ppditsd = value.get_as<double>().value();
            mod.BSIM4ppditsdGiven = true;
            break;
        case  BSIM4_MOD_PDELTA :
            mod.BSIM4pdelta = value.get_as<double>().value();
            mod.BSIM4pdeltaGiven = true;
            break;
        case BSIM4_MOD_PRDSW:
            mod.BSIM4prdsw = value.get_as<double>().value();
            mod.BSIM4prdswGiven = true;
            break;
        case BSIM4_MOD_PRDW:
            mod.BSIM4prdw = value.get_as<double>().value();
            mod.BSIM4prdwGiven = true;
            break;
        case BSIM4_MOD_PRSW:
            mod.BSIM4prsw = value.get_as<double>().value();
            mod.BSIM4prswGiven = true;
            break;
        case BSIM4_MOD_PPRWB:
            mod.BSIM4pprwb = value.get_as<double>().value();
            mod.BSIM4pprwbGiven = true;
            break;
        case BSIM4_MOD_PPRWG:
            mod.BSIM4pprwg = value.get_as<double>().value();
            mod.BSIM4pprwgGiven = true;
            break;
        case BSIM4_MOD_PPRT:
            mod.BSIM4pprt = value.get_as<double>().value();
            mod.BSIM4pprtGiven = true;
            break;
        case BSIM4_MOD_PETA0:
            mod.BSIM4peta0 = value.get_as<double>().value();
            mod.BSIM4peta0Given = true;
            break;
        case BSIM4_MOD_PETAB:
            mod.BSIM4petab = value.get_as<double>().value();
            mod.BSIM4petabGiven = true;
            break;
        case BSIM4_MOD_PPCLM:
            mod.BSIM4ppclm = value.get_as<double>().value();
            mod.BSIM4ppclmGiven = true;
            break;
        case BSIM4_MOD_PPDIBL1:
            mod.BSIM4ppdibl1 = value.get_as<double>().value();
            mod.BSIM4ppdibl1Given = true;
            break;
        case BSIM4_MOD_PPDIBL2:
            mod.BSIM4ppdibl2 = value.get_as<double>().value();
            mod.BSIM4ppdibl2Given = true;
            break;
        case BSIM4_MOD_PPDIBLB:
            mod.BSIM4ppdiblb = value.get_as<double>().value();
            mod.BSIM4ppdiblbGiven = true;
            break;
        case BSIM4_MOD_PPSCBE1:
            mod.BSIM4ppscbe1 = value.get_as<double>().value();
            mod.BSIM4ppscbe1Given = true;
            break;
        case BSIM4_MOD_PPSCBE2:
            mod.BSIM4ppscbe2 = value.get_as<double>().value();
            mod.BSIM4ppscbe2Given = true;
            break;
        case BSIM4_MOD_PPVAG:
            mod.BSIM4ppvag = value.get_as<double>().value();
            mod.BSIM4ppvagGiven = true;
            break;
        case  BSIM4_MOD_PWR :
            mod.BSIM4pwr = value.get_as<double>().value();
            mod.BSIM4pwrGiven = true;
            break;
        case  BSIM4_MOD_PDWG :
            mod.BSIM4pdwg = value.get_as<double>().value();
            mod.BSIM4pdwgGiven = true;
            break;
        case  BSIM4_MOD_PDWB :
            mod.BSIM4pdwb = value.get_as<double>().value();
            mod.BSIM4pdwbGiven = true;
            break;
        case  BSIM4_MOD_PB0 :
            mod.BSIM4pb0 = value.get_as<double>().value();
            mod.BSIM4pb0Given = true;
            break;
        case  BSIM4_MOD_PB1 :
            mod.BSIM4pb1 = value.get_as<double>().value();
            mod.BSIM4pb1Given = true;
            break;
        case  BSIM4_MOD_PALPHA0 :
            mod.BSIM4palpha0 = value.get_as<double>().value();
            mod.BSIM4palpha0Given = true;
            break;
        case  BSIM4_MOD_PALPHA1 :
            mod.BSIM4palpha1 = value.get_as<double>().value();
            mod.BSIM4palpha1Given = true;
            break;
        case  BSIM4_MOD_PBETA0 :
            mod.BSIM4pbeta0 = value.get_as<double>().value();
            mod.BSIM4pbeta0Given = true;
            break;
        case  BSIM4_MOD_PPHIN :
            mod.BSIM4pphin = value.get_as<double>().value();
            mod.BSIM4pphinGiven = true;
            break;
        case  BSIM4_MOD_PAGIDL :
            mod.BSIM4pagidl = value.get_as<double>().value();
            mod.BSIM4pagidlGiven = true;
            break;
        case  BSIM4_MOD_PBGIDL :
            mod.BSIM4pbgidl = value.get_as<double>().value();
            mod.BSIM4pbgidlGiven = true;
            break;
        case  BSIM4_MOD_PCGIDL :
            mod.BSIM4pcgidl = value.get_as<double>().value();
            mod.BSIM4pcgidlGiven = true;
            break;
        case  BSIM4_MOD_PEGIDL :
            mod.BSIM4pegidl = value.get_as<double>().value();
            mod.BSIM4pegidlGiven = true;
            break;
    case  BSIM4_MOD_PFGIDL :            /* v4.7 New GIDL/GISL */
            mod.BSIM4pfgidl = value.get_as<double>().value();
            mod.BSIM4pfgidlGiven = true;
            break;
    case  BSIM4_MOD_PKGIDL :            /* v4.7 New GIDL/GISL */
            mod.BSIM4pkgidl = value.get_as<double>().value();
            mod.BSIM4pkgidlGiven = true;
            break;
    case  BSIM4_MOD_PRGIDL :            /* v4.7 New GIDL/GISL */
            mod.BSIM4prgidl = value.get_as<double>().value();
            mod.BSIM4prgidlGiven = true;
            break;
        case  BSIM4_MOD_PAGISL :
            mod.BSIM4pagisl = value.get_as<double>().value();
            mod.BSIM4pagislGiven = true;
            break;
        case  BSIM4_MOD_PBGISL :
            mod.BSIM4pbgisl = value.get_as<double>().value();
            mod.BSIM4pbgislGiven = true;
            break;
        case  BSIM4_MOD_PCGISL :
            mod.BSIM4pcgisl = value.get_as<double>().value();
            mod.BSIM4pcgislGiven = true;
            break;
        case  BSIM4_MOD_PEGISL :
            mod.BSIM4pegisl = value.get_as<double>().value();
            mod.BSIM4pegislGiven = true;
            break;
    case  BSIM4_MOD_PFGISL :            /* v4.7 New GIDL/GISL */
            mod.BSIM4pfgisl = value.get_as<double>().value();
            mod.BSIM4pfgislGiven = true;
            break;
    case  BSIM4_MOD_PKGISL :            /* v4.7 New GIDL/GISL */
            mod.BSIM4pkgisl = value.get_as<double>().value();
            mod.BSIM4pkgislGiven = true;
            break;
    case  BSIM4_MOD_PRGISL :            /* v4.7 New GIDL/GISL */
            mod.BSIM4prgisl = value.get_as<double>().value();
            mod.BSIM4prgislGiven = true;
            break;
        case  BSIM4_MOD_PAIGC :
            mod.BSIM4paigc = value.get_as<double>().value();
            mod.BSIM4paigcGiven = true;
            break;
        case  BSIM4_MOD_PBIGC :
            mod.BSIM4pbigc = value.get_as<double>().value();
            mod.BSIM4pbigcGiven = true;
            break;
        case  BSIM4_MOD_PCIGC :
            mod.BSIM4pcigc = value.get_as<double>().value();
            mod.BSIM4pcigcGiven = true;
            break;
        case  BSIM4_MOD_PAIGSD :
            mod.BSIM4paigsd = value.get_as<double>().value();
            mod.BSIM4paigsdGiven = true;
            break;
        case  BSIM4_MOD_PBIGSD :
            mod.BSIM4pbigsd = value.get_as<double>().value();
            mod.BSIM4pbigsdGiven = true;
            break;
        case  BSIM4_MOD_PCIGSD :
            mod.BSIM4pcigsd = value.get_as<double>().value();
            mod.BSIM4pcigsdGiven = true;
            break;
        case  BSIM4_MOD_PAIGS :
            mod.BSIM4paigs = value.get_as<double>().value();
            mod.BSIM4paigsGiven = true;
            break;
        case  BSIM4_MOD_PBIGS :
            mod.BSIM4pbigs = value.get_as<double>().value();
            mod.BSIM4pbigsGiven = true;
            break;
        case  BSIM4_MOD_PCIGS :
            mod.BSIM4pcigs = value.get_as<double>().value();
            mod.BSIM4pcigsGiven = true;
            break;
        case  BSIM4_MOD_PAIGD :
            mod.BSIM4paigd = value.get_as<double>().value();
            mod.BSIM4paigdGiven = true;
            break;
        case  BSIM4_MOD_PBIGD :
            mod.BSIM4pbigd = value.get_as<double>().value();
            mod.BSIM4pbigdGiven = true;
            break;
        case  BSIM4_MOD_PCIGD :
            mod.BSIM4pcigd = value.get_as<double>().value();
            mod.BSIM4pcigdGiven = true;
            break;
        case  BSIM4_MOD_PAIGBACC :
            mod.BSIM4paigbacc = value.get_as<double>().value();
            mod.BSIM4paigbaccGiven = true;
            break;
        case  BSIM4_MOD_PBIGBACC :
            mod.BSIM4pbigbacc = value.get_as<double>().value();
            mod.BSIM4pbigbaccGiven = true;
            break;
        case  BSIM4_MOD_PCIGBACC :
            mod.BSIM4pcigbacc = value.get_as<double>().value();
            mod.BSIM4pcigbaccGiven = true;
            break;
        case  BSIM4_MOD_PAIGBINV :
            mod.BSIM4paigbinv = value.get_as<double>().value();
            mod.BSIM4paigbinvGiven = true;
            break;
        case  BSIM4_MOD_PBIGBINV :
            mod.BSIM4pbigbinv = value.get_as<double>().value();
            mod.BSIM4pbigbinvGiven = true;
            break;
        case  BSIM4_MOD_PCIGBINV :
            mod.BSIM4pcigbinv = value.get_as<double>().value();
            mod.BSIM4pcigbinvGiven = true;
            break;
        case  BSIM4_MOD_PNIGC :
            mod.BSIM4pnigc = value.get_as<double>().value();
            mod.BSIM4pnigcGiven = true;
            break;
        case  BSIM4_MOD_PNIGBINV :
            mod.BSIM4pnigbinv = value.get_as<double>().value();
            mod.BSIM4pnigbinvGiven = true;
            break;
        case  BSIM4_MOD_PNIGBACC :
            mod.BSIM4pnigbacc = value.get_as<double>().value();
            mod.BSIM4pnigbaccGiven = true;
            break;
        case  BSIM4_MOD_PNTOX :
            mod.BSIM4pntox = value.get_as<double>().value();
            mod.BSIM4pntoxGiven = true;
            break;
        case  BSIM4_MOD_PEIGBINV :
            mod.BSIM4peigbinv = value.get_as<double>().value();
            mod.BSIM4peigbinvGiven = true;
            break;
        case  BSIM4_MOD_PPIGCD :
            mod.BSIM4ppigcd = value.get_as<double>().value();
            mod.BSIM4ppigcdGiven = true;
            break;
        case  BSIM4_MOD_PPOXEDGE :
            mod.BSIM4ppoxedge = value.get_as<double>().value();
            mod.BSIM4ppoxedgeGiven = true;
            break;
        case  BSIM4_MOD_PXRCRG1 :
            mod.BSIM4pxrcrg1 = value.get_as<double>().value();
            mod.BSIM4pxrcrg1Given = true;
            break;
        case  BSIM4_MOD_PXRCRG2 :
            mod.BSIM4pxrcrg2 = value.get_as<double>().value();
            mod.BSIM4pxrcrg2Given = true;
            break;
        case  BSIM4_MOD_PLAMBDA :
            mod.BSIM4plambda = value.get_as<double>().value();
            mod.BSIM4plambdaGiven = true;
            break;
        case  BSIM4_MOD_PVTL :
            mod.BSIM4pvtl = value.get_as<double>().value();
            mod.BSIM4pvtlGiven = true;
            break;
        case  BSIM4_MOD_PXN:
            mod.BSIM4pxn = value.get_as<double>().value();
            mod.BSIM4pxnGiven = true;
            break;
        case  BSIM4_MOD_PVFBSDOFF:
            mod.BSIM4pvfbsdoff = value.get_as<double>().value();
            mod.BSIM4pvfbsdoffGiven = true;
            break;
        case  BSIM4_MOD_PTVFBSDOFF:
            mod.BSIM4ptvfbsdoff = value.get_as<double>().value();
            mod.BSIM4ptvfbsdoffGiven = true;
            break;
        case  BSIM4_MOD_PEU :
            mod.BSIM4peu = value.get_as<double>().value();
            mod.BSIM4peuGiven = true;
            break;
        case  BSIM4_MOD_PUCS :
            mod.BSIM4pucs = value.get_as<double>().value();
            mod.BSIM4pucsGiven = true;
            break;
        case  BSIM4_MOD_PVFB :
            mod.BSIM4pvfb = value.get_as<double>().value();
            mod.BSIM4pvfbGiven = true;
            break;
        case  BSIM4_MOD_PCGSL :
            mod.BSIM4pcgsl = value.get_as<double>().value();
            mod.BSIM4pcgslGiven = true;
            break;
        case  BSIM4_MOD_PCGDL :
            mod.BSIM4pcgdl = value.get_as<double>().value();
            mod.BSIM4pcgdlGiven = true;
            break;
        case  BSIM4_MOD_PCKAPPAS :
            mod.BSIM4pckappas = value.get_as<double>().value();
            mod.BSIM4pckappasGiven = true;
            break;
        case  BSIM4_MOD_PCKAPPAD :
            mod.BSIM4pckappad = value.get_as<double>().value();
            mod.BSIM4pckappadGiven = true;
            break;
        case  BSIM4_MOD_PCF :
            mod.BSIM4pcf = value.get_as<double>().value();
            mod.BSIM4pcfGiven = true;
            break;
        case  BSIM4_MOD_PCLC :
            mod.BSIM4pclc = value.get_as<double>().value();
            mod.BSIM4pclcGiven = true;
            break;
        case  BSIM4_MOD_PCLE :
            mod.BSIM4pcle = value.get_as<double>().value();
            mod.BSIM4pcleGiven = true;
            break;
        case  BSIM4_MOD_PVFBCV :
            mod.BSIM4pvfbcv = value.get_as<double>().value();
            mod.BSIM4pvfbcvGiven = true;
            break;
        case  BSIM4_MOD_PACDE :
            mod.BSIM4pacde = value.get_as<double>().value();
            mod.BSIM4pacdeGiven = true;
            break;
        case  BSIM4_MOD_PMOIN :
            mod.BSIM4pmoin = value.get_as<double>().value();
            mod.BSIM4pmoinGiven = true;
            break;
        case  BSIM4_MOD_PNOFF :
            mod.BSIM4pnoff = value.get_as<double>().value();
            mod.BSIM4pnoffGiven = true;
            break;
        case  BSIM4_MOD_PVOFFCV :
            mod.BSIM4pvoffcv = value.get_as<double>().value();
            mod.BSIM4pvoffcvGiven = true;
            break;

        case  BSIM4_MOD_TNOM :
            mod.BSIM4tnom = value.get_as<double>().value() + CONSTCtoK;
            mod.BSIM4tnomGiven = true;
            break;
        case  BSIM4_MOD_CGSO :
            mod.BSIM4cgso = value.get_as<double>().value();
            mod.BSIM4cgsoGiven = true;
            break;
        case  BSIM4_MOD_CGDO :
            mod.BSIM4cgdo = value.get_as<double>().value();
            mod.BSIM4cgdoGiven = true;
            break;
        case  BSIM4_MOD_CGBO :
            mod.BSIM4cgbo = value.get_as<double>().value();
            mod.BSIM4cgboGiven = true;
            break;
        case  BSIM4_MOD_XPART :
            mod.BSIM4xpart = value.get_as<double>().value();
            mod.BSIM4xpartGiven = true;
            break;
        case  BSIM4_MOD_RSH :
            mod.BSIM4sheetResistance = value.get_as<double>().value();
            mod.BSIM4sheetResistanceGiven = true;
            break;
        case  BSIM4_MOD_JSS :
            mod.BSIM4SjctSatCurDensity = value.get_as<double>().value();
            mod.BSIM4SjctSatCurDensityGiven = true;
            break;
        case  BSIM4_MOD_JSWS :
            mod.BSIM4SjctSidewallSatCurDensity = value.get_as<double>().value();
            mod.BSIM4SjctSidewallSatCurDensityGiven = true;
            break;
        case  BSIM4_MOD_JSWGS :
            mod.BSIM4SjctGateSidewallSatCurDensity = value.get_as<double>().value();
            mod.BSIM4SjctGateSidewallSatCurDensityGiven = true;
            break;
        case  BSIM4_MOD_PBS :
            mod.BSIM4SbulkJctPotential = value.get_as<double>().value();
            mod.BSIM4SbulkJctPotentialGiven = true;
            break;
        case  BSIM4_MOD_MJS :
            mod.BSIM4SbulkJctBotGradingCoeff = value.get_as<double>().value();
            mod.BSIM4SbulkJctBotGradingCoeffGiven = true;
            break;
        case  BSIM4_MOD_PBSWS :
            mod.BSIM4SsidewallJctPotential = value.get_as<double>().value();
            mod.BSIM4SsidewallJctPotentialGiven = true;
            break;
        case  BSIM4_MOD_MJSWS :
            mod.BSIM4SbulkJctSideGradingCoeff = value.get_as<double>().value();
            mod.BSIM4SbulkJctSideGradingCoeffGiven = true;
            break;
        case  BSIM4_MOD_CJS :
            mod.BSIM4SunitAreaJctCap = value.get_as<double>().value();
            mod.BSIM4SunitAreaJctCapGiven = true;
            break;
        case  BSIM4_MOD_CJSWS :
            mod.BSIM4SunitLengthSidewallJctCap = value.get_as<double>().value();
            mod.BSIM4SunitLengthSidewallJctCapGiven = true;
            break;
        case  BSIM4_MOD_NJS :
            mod.BSIM4SjctEmissionCoeff = value.get_as<double>().value();
            mod.BSIM4SjctEmissionCoeffGiven = true;
            break;
        case  BSIM4_MOD_PBSWGS :
            mod.BSIM4SGatesidewallJctPotential = value.get_as<double>().value();
            mod.BSIM4SGatesidewallJctPotentialGiven = true;
            break;
        case  BSIM4_MOD_MJSWGS :
            mod.BSIM4SbulkJctGateSideGradingCoeff = value.get_as<double>().value();
            mod.BSIM4SbulkJctGateSideGradingCoeffGiven = true;
            break;
        case  BSIM4_MOD_CJSWGS :
            mod.BSIM4SunitLengthGateSidewallJctCap = value.get_as<double>().value();
            mod.BSIM4SunitLengthGateSidewallJctCapGiven = true;
            break;
        case  BSIM4_MOD_XTIS :
            mod.BSIM4SjctTempExponent = value.get_as<double>().value();
            mod.BSIM4SjctTempExponentGiven = true;
            break;
        case  BSIM4_MOD_JSD :
            mod.BSIM4DjctSatCurDensity = value.get_as<double>().value();
            mod.BSIM4DjctSatCurDensityGiven = true;
            break;
        case  BSIM4_MOD_JSWD :
            mod.BSIM4DjctSidewallSatCurDensity = value.get_as<double>().value();
            mod.BSIM4DjctSidewallSatCurDensityGiven = true;
            break;
        case  BSIM4_MOD_JSWGD :
            mod.BSIM4DjctGateSidewallSatCurDensity = value.get_as<double>().value();
            mod.BSIM4DjctGateSidewallSatCurDensityGiven = true;
            break;
        case  BSIM4_MOD_PBD :
            mod.BSIM4DbulkJctPotential = value.get_as<double>().value();
            mod.BSIM4DbulkJctPotentialGiven = true;
            break;
        case  BSIM4_MOD_MJD :
            mod.BSIM4DbulkJctBotGradingCoeff = value.get_as<double>().value();
            mod.BSIM4DbulkJctBotGradingCoeffGiven = true;
            break;
        case  BSIM4_MOD_PBSWD :
            mod.BSIM4DsidewallJctPotential = value.get_as<double>().value();
            mod.BSIM4DsidewallJctPotentialGiven = true;
            break;
        case  BSIM4_MOD_MJSWD :
            mod.BSIM4DbulkJctSideGradingCoeff = value.get_as<double>().value();
            mod.BSIM4DbulkJctSideGradingCoeffGiven = true;
            break;
        case  BSIM4_MOD_CJD :
            mod.BSIM4DunitAreaJctCap = value.get_as<double>().value();
            mod.BSIM4DunitAreaJctCapGiven = true;
            break;
        case  BSIM4_MOD_CJSWD :
            mod.BSIM4DunitLengthSidewallJctCap = value.get_as<double>().value();
            mod.BSIM4DunitLengthSidewallJctCapGiven = true;
            break;
        case  BSIM4_MOD_NJD :
            mod.BSIM4DjctEmissionCoeff = value.get_as<double>().value();
            mod.BSIM4DjctEmissionCoeffGiven = true;
            break;
        case  BSIM4_MOD_PBSWGD :
            mod.BSIM4DGatesidewallJctPotential = value.get_as<double>().value();
            mod.BSIM4DGatesidewallJctPotentialGiven = true;
            break;
        case  BSIM4_MOD_MJSWGD :
            mod.BSIM4DbulkJctGateSideGradingCoeff = value.get_as<double>().value();
            mod.BSIM4DbulkJctGateSideGradingCoeffGiven = true;
            break;
        case  BSIM4_MOD_CJSWGD :
            mod.BSIM4DunitLengthGateSidewallJctCap = value.get_as<double>().value();
            mod.BSIM4DunitLengthGateSidewallJctCapGiven = true;
            break;
        case  BSIM4_MOD_XTID :
            mod.BSIM4DjctTempExponent = value.get_as<double>().value();
            mod.BSIM4DjctTempExponentGiven = true;
            break;
        case  BSIM4_MOD_LINT :
            mod.BSIM4Lint = value.get_as<double>().value();
            mod.BSIM4LintGiven = true;
            break;
        case  BSIM4_MOD_LL :
            mod.BSIM4Ll = value.get_as<double>().value();
            mod.BSIM4LlGiven = true;
            break;
        case  BSIM4_MOD_LLC :
            mod.BSIM4Llc = value.get_as<double>().value();
            mod.BSIM4LlcGiven = true;
            break;
        case  BSIM4_MOD_LLN :
            mod.BSIM4Lln = value.get_as<double>().value();
            mod.BSIM4LlnGiven = true;
            break;
        case  BSIM4_MOD_LW :
            mod.BSIM4Lw = value.get_as<double>().value();
            mod.BSIM4LwGiven = true;
            break;
        case  BSIM4_MOD_LWC :
            mod.BSIM4Lwc = value.get_as<double>().value();
            mod.BSIM4LwcGiven = true;
            break;
        case  BSIM4_MOD_LWN :
            mod.BSIM4Lwn = value.get_as<double>().value();
            mod.BSIM4LwnGiven = true;
            break;
        case  BSIM4_MOD_LWL :
            mod.BSIM4Lwl = value.get_as<double>().value();
            mod.BSIM4LwlGiven = true;
            break;
        case  BSIM4_MOD_LWLC :
            mod.BSIM4Lwlc = value.get_as<double>().value();
            mod.BSIM4LwlcGiven = true;
            break;
        case  BSIM4_MOD_LMIN :
            mod.BSIM4Lmin = value.get_as<double>().value();
            mod.BSIM4LminGiven = true;
            break;
        case  BSIM4_MOD_LMAX :
            mod.BSIM4Lmax = value.get_as<double>().value();
            mod.BSIM4LmaxGiven = true;
            break;
        case  BSIM4_MOD_WINT :
            mod.BSIM4Wint = value.get_as<double>().value();
            mod.BSIM4WintGiven = true;
            break;
        case  BSIM4_MOD_WL :
            mod.BSIM4Wl = value.get_as<double>().value();
            mod.BSIM4WlGiven = true;
            break;
        case  BSIM4_MOD_WLC :
            mod.BSIM4Wlc = value.get_as<double>().value();
            mod.BSIM4WlcGiven = true;
            break;
        case  BSIM4_MOD_WLN :
            mod.BSIM4Wln = value.get_as<double>().value();
            mod.BSIM4WlnGiven = true;
            break;
        case  BSIM4_MOD_WW :
            mod.BSIM4Ww = value.get_as<double>().value();
            mod.BSIM4WwGiven = true;
            break;
        case  BSIM4_MOD_WWC :
            mod.BSIM4Wwc = value.get_as<double>().value();
            mod.BSIM4WwcGiven = true;
            break;
        case  BSIM4_MOD_WWN :
            mod.BSIM4Wwn = value.get_as<double>().value();
            mod.BSIM4WwnGiven = true;
            break;
        case  BSIM4_MOD_WWL :
            mod.BSIM4Wwl = value.get_as<double>().value();
            mod.BSIM4WwlGiven = true;
            break;
        case  BSIM4_MOD_WWLC :
            mod.BSIM4Wwlc = value.get_as<double>().value();
            mod.BSIM4WwlcGiven = true;
            break;
        case  BSIM4_MOD_WMIN :
            mod.BSIM4Wmin = value.get_as<double>().value();
            mod.BSIM4WminGiven = true;
            break;
        case  BSIM4_MOD_WMAX :
            mod.BSIM4Wmax = value.get_as<double>().value();
            mod.BSIM4WmaxGiven = true;
            break;

        case  BSIM4_MOD_NOIA :
            mod.BSIM4oxideTrapDensityA = value.get_as<double>().value();
            mod.BSIM4oxideTrapDensityAGiven = true;
            break;
        case  BSIM4_MOD_NOIB :
            mod.BSIM4oxideTrapDensityB = value.get_as<double>().value();
            mod.BSIM4oxideTrapDensityBGiven = true;
            break;
        case  BSIM4_MOD_NOIC :
            mod.BSIM4oxideTrapDensityC = value.get_as<double>().value();
            mod.BSIM4oxideTrapDensityCGiven = true;
            break;
        case  BSIM4_MOD_EM :
            mod.BSIM4em = value.get_as<double>().value();
            mod.BSIM4emGiven = true;
            break;
        case  BSIM4_MOD_EF :
            mod.BSIM4ef = value.get_as<double>().value();
            mod.BSIM4efGiven = true;
            break;
        case  BSIM4_MOD_AF :
            mod.BSIM4af = value.get_as<double>().value();
            mod.BSIM4afGiven = true;
            break;
        case  BSIM4_MOD_KF :
            mod.BSIM4kf = value.get_as<double>().value();
            mod.BSIM4kfGiven = true;
            break;
        case  BSIM4_MOD_NMOS  :
            mod.BSIM4type = 1;
            mod.BSIM4typeGiven = true;
            break;
        case  BSIM4_MOD_PMOS  :
            mod.BSIM4type = - 1;
            mod.BSIM4typeGiven = true;
            break;
        default:
            std::cerr << "BSIM4param: unknown parameter " << param << std::endl;
            return(false);
    }
    return(true);
}

// 2. set instance parameters
int BSIM4param(int param, const VariantValue &value, BSIM4V82 &inst){
    switch(param)
    {   case BSIM4_W:
            inst.BSIM4w = value.get_as<double>().value();
            inst.BSIM4wGiven = true;
            break;
        case BSIM4_L:
            inst.BSIM4l = value.get_as<double>().value();
            inst.BSIM4lGiven = true;
            break;
        case BSIM4_NF:
            inst.BSIM4nf = value.get_as<double>().value();
            inst.BSIM4nfGiven = true;
            break;
        case BSIM4_MIN:
            inst.BSIM4min = value.get_as<int>().value();
            inst.BSIM4minGiven = true;
            break;
        case BSIM4_AS:
            inst.BSIM4sourceArea = value.get_as<double>().value();
            inst.BSIM4sourceAreaGiven = true;
            break;
        case BSIM4_AD:
            inst.BSIM4drainArea = value.get_as<double>().value();
            inst.BSIM4drainAreaGiven = true;
            break;
        case BSIM4_PS:
            inst.BSIM4sourcePerimeter = value.get_as<double>().value();
            inst.BSIM4sourcePerimeterGiven = true;
            break;
        case BSIM4_PD:
            inst.BSIM4drainPerimeter = value.get_as<double>().value();
            inst.BSIM4drainPerimeterGiven = true;
            break;
        case BSIM4_NRS:
            inst.BSIM4sourceSquares = value.get_as<double>().value();
            inst.BSIM4sourceSquaresGiven = true;
            break;
        case BSIM4_NRD:
            inst.BSIM4drainSquares = value.get_as<double>().value();
            inst.BSIM4drainSquaresGiven = true;
            break;
        case BSIM4_OFF:
            inst.BSIM4off = value.get_as<int>().value();
            break;
        case BSIM4_SA:
            inst.BSIM4sa = value.get_as<double>().value();
            inst.BSIM4saGiven = true;
            break;
        case BSIM4_SB:
            inst.BSIM4sb = value.get_as<double>().value();
            inst.BSIM4sbGiven = true;
            break;
        case BSIM4_SD:
            inst.BSIM4sd = value.get_as<double>().value();
            inst.BSIM4sdGiven = true;
            break;
        case BSIM4_SCA:
            inst.BSIM4sca = value.get_as<double>().value();
            inst.BSIM4scaGiven = true;
            break;
        case BSIM4_SCB:
            inst.BSIM4scb = value.get_as<double>().value();
            inst.BSIM4scbGiven = true;
            break;
        case BSIM4_SCC:
            inst.BSIM4scc = value.get_as<double>().value();
            inst.BSIM4sccGiven = true;
            break;
        case BSIM4_SC:
            inst.BSIM4sc = value.get_as<double>().value();
            inst.BSIM4scGiven = true;
            break;
        case BSIM4_RBSB:
            inst.BSIM4rbsb = value.get_as<double>().value();
            inst.BSIM4rbsbGiven = true;
            break;
        case BSIM4_RBDB:
            inst.BSIM4rbdb = value.get_as<double>().value();
            inst.BSIM4rbdbGiven = true;
            break;
        case BSIM4_RBPB:
            inst.BSIM4rbpb = value.get_as<double>().value();
            inst.BSIM4rbpbGiven = true;
            break;
        case BSIM4_RBPS:
            inst.BSIM4rbps = value.get_as<double>().value();
            inst.BSIM4rbpsGiven = true;
            break;
        case BSIM4_RBPD:
            inst.BSIM4rbpd = value.get_as<double>().value();
            inst.BSIM4rbpdGiven = true;
            break;
        case BSIM4_DELVTO:
            inst.BSIM4delvto = value.get_as<double>().value();
            inst.BSIM4delvtoGiven = true;
            break;
        case BSIM4_XGW:
            inst.BSIM4xgw = value.get_as<double>().value();
            inst.BSIM4xgwGiven = true;
            break;
        case BSIM4_NGCON:
            inst.BSIM4ngcon = value.get_as<double>().value();
            inst.BSIM4ngconGiven = true;
            break;
        case BSIM4_TRNQSMOD:
            inst.BSIM4trnqsMod = value.get_as<int>().value();
            inst.BSIM4trnqsModGiven = true;
            break;
        case BSIM4_ACNQSMOD:
            inst.BSIM4acnqsMod = value.get_as<int>().value();
            inst.BSIM4acnqsModGiven = true;
            break;
        case BSIM4_RBODYMOD:
            inst.BSIM4rbodyMod = value.get_as<int>().value();
            inst.BSIM4rbodyModGiven = true;
            break;
        case BSIM4_RGATEMOD:
            inst.BSIM4rgateMod = value.get_as<int>().value();
            inst.BSIM4rgateModGiven = true;
            break;
        case BSIM4_GEOMOD:
            inst.BSIM4geoMod = value.get_as<int>().value();
            inst.BSIM4geoModGiven = true;
            break;
        case BSIM4_RGEOMOD:
            inst.BSIM4rgeoMod = value.get_as<int>().value();
            inst.BSIM4rgeoModGiven = true;
            break;
        case BSIM4_IC_VDS:
            inst.BSIM4icVDS = value.get_as<double>().value();
            inst.BSIM4icVDSGiven = true;
            break;
        case BSIM4_IC_VGS:
            inst.BSIM4icVGS = value.get_as<double>().value();
            inst.BSIM4icVGSGiven = true;
            break;
        case BSIM4_IC_VBS:
            inst.BSIM4icVBS = value.get_as<double>().value();
            inst.BSIM4icVBSGiven = true;
            break;
       //   case BSIM4_IC:
       //       switch(value->v.numValue)
       //       {   case 3:
       //               inst.BSIM4icVBS = *(value->v.vec.rVec+2);
       //               inst.BSIM4icVBSGiven = true;
       //           case 2:
       //               inst.BSIM4icVGS = *(value->v.vec.rVec+1);
       //               inst.BSIM4icVGSGiven = true;
       //           case 1:
       //               inst.BSIM4icVDS = *(value->v.vec.rVec);
       //               inst.BSIM4icVDSGiven = true;
       //               break;
       //           default:
       //                std::cerr << "BSIM4param: illegal number of values for IC" << std::endl;
       //                return(false);
       //       }
       //       break;
        default:
            std::cerr << "BSIM4param: unknown parameter " << param << std::endl;
            return(false);
    }
    return(true);
}

// Paser for BSIM4 model
std::shared_ptr<BSIM4model> paserBSIM4Model(const std::string& modelType, const std::string& modelName, const std::map<std::string, std::string>& kvMap){
    std::shared_ptr<BSIM4model> model = std::make_shared<BSIM4model>();
    // Convert modelType to lowercase for case-insensitive comparison
    std::string typeLower = modelType;
    std::transform(typeLower.begin(), typeLower.end(), typeLower.begin(),
        [](unsigned char c) { return std::tolower(c); });
    
    // Set the type based on "nmos" or "pmos"
    if (typeLower == "nmos") {
        model->BSIM4type = BSIM4_NMOS;
        model->BSIM4typeGiven = true;
    } else if (typeLower == "pmos") {
        model->BSIM4type = BSIM4_PMOS;
        model->BSIM4typeGiven = true;
    } else {
        std::cerr << "Error: Invalid model type for BSIM4: " << modelType << std::endl;
        exit(1);
    }
    // Set the model name
    model->BSIM4modName = modelName;

    // Parse parameters from kvMap
    for (const auto& [key, valStr] : kvMap) {
        auto iter = bsim4ParamMap.find(key);
        if (iter != bsim4ParamMap.end()){
            int paramIndex = iter->second;
            try {
                VariantValue value(valStr);
                BSIM4mParam(paramIndex, value, *model);
            } catch (const std::exception& e) {
                std::cerr << "Warning: Failed to parse BSIM4 parameter " << key << " = " << valStr
                          << ": " << e.what() << std::endl;
            }
        } 
        else {
            std::cerr << "Warning: Unknown BSIM4 parameter: " << key << std::endl;
        }
    }

    return model;
}


} //namespace bsim4

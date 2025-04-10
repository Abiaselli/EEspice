#pragma once
#include <string>
#include <iostream>
#include <cmath>
#include <memory>
#include "bsim4v82.hpp"
#include "bsim4v82geo.hpp"
#include "bsim4v82check.hpp"
#include "sim_variables.hpp"
#include "bsim4v82const.hpp"


// #define DEXP(A,B) {                                                        \
//         if (A > EXP_THRESHOLD) {                                           \
//             B = MAX_EXP*(1.0+(A)-EXP_THRESHOLD);                           \
//         } else if (A < -EXP_THRESHOLD)  {                                  \
//             B = MIN_EXP;                                                   \
//         } else   {                                                         \
//             B = std::exp(A);                                                    \
//         }                                                                  \
//     }
// #define ABS(A)   (((A) > 0) ? (A) : (-(A)))

namespace bsim4{

inline double DEXP(double A){
    if (A > EXP_THRESHOLD) {
        return MAX_EXP*(1.0+(A)-EXP_THRESHOLD);
    } else if (A < -EXP_THRESHOLD)  {
        return MIN_EXP;
    } else   {
        return std::exp(A);
    }
}

int BSIM4DioIjthVjmEval(double Nvtm, double Ijth, double Isb, double XExpBV, double *Vjm)
{
    double Tb, Tc, EVjmovNv;
       Tc = XExpBV;
       Tb = 1.0 + Ijth / Isb - Tc;
       EVjmovNv = 0.5 * (Tb + std::sqrt(Tb * Tb + 4.0 * Tc));
       *Vjm = Nvtm * std::log(EVjmovNv);
    return 0;
}


/*  loop through all the BSIM4 device models in eespice not here */
int modelTemp(BSIM4model &model, double CKTtemp){
    // All Temporary parameters for this model 
    bsim4v82temp &BSIM4temp = model.BSIM4temp;

    BSIM4temp.Temp = CKTtemp;
         if (model.BSIM4SbulkJctPotential < 0.1)
     {   model.BSIM4SbulkJctPotential = 0.1;
         fprintf(stderr, "Given pbs is less than 0.1. Pbs is set to 0.1.\n");
     }
         if (model.BSIM4SsidewallJctPotential < 0.1)
     {   model.BSIM4SsidewallJctPotential = 0.1;
         fprintf(stderr, "Given pbsws is less than 0.1. Pbsws is set to 0.1.\n");
     }
         if (model.BSIM4SGatesidewallJctPotential < 0.1)
     {   model.BSIM4SGatesidewallJctPotential = 0.1;
         fprintf(stderr, "Given pbswgs is less than 0.1. Pbswgs is set to 0.1.\n");
     }

         if (model.BSIM4DbulkJctPotential < 0.1)
         {   model.BSIM4DbulkJctPotential = 0.1;
             fprintf(stderr, "Given pbd is less than 0.1. Pbd is set to 0.1.\n");
         }
         if (model.BSIM4DsidewallJctPotential < 0.1)
         {   model.BSIM4DsidewallJctPotential = 0.1;
             fprintf(stderr, "Given pbswd is less than 0.1. Pbswd is set to 0.1.\n");
         }
         if (model.BSIM4DGatesidewallJctPotential < 0.1)
         {   model.BSIM4DGatesidewallJctPotential = 0.1;
             fprintf(stderr, "Given pbswgd is less than 0.1. Pbswgd is set to 0.1.\n");
         }

     if(model.BSIM4mtrlMod == 0)
     {
         if ((model.BSIM4toxeGiven) && (model.BSIM4toxpGiven) && (model.BSIM4dtoxGiven)
         && (model.BSIM4toxe != (model.BSIM4toxp + model.BSIM4dtox)))
         {   printf("Warning:  BSIM4temp.toxe, toxp and dtox all given and  BSIM4temp.toxe != toxp + dtox; dtox ignored.\n");
         }
         else if ((model.BSIM4toxeGiven) && (!model.BSIM4toxpGiven))
         { model.BSIM4toxp = model.BSIM4toxe - model.BSIM4dtox;
         }
         else if ((!model.BSIM4toxeGiven) && (model.BSIM4toxpGiven))
         {
             model.BSIM4toxe = model.BSIM4toxp + model.BSIM4dtox;
             if (!model.BSIM4toxmGiven)            /* v4.7 */
                 model.BSIM4toxm = model.BSIM4toxe; 				   
         }
             if (!model.BSIM4cfGiven)            /* v4.8.2 */
                  model.BSIM4cf = 2.0 * model.BSIM4epsrox * EPS0 / PI
                   * std::log(1.0 + 0.4e-6 / model.BSIM4toxe);
     }

     else if(model.BSIM4mtrlCompatMod != 0) /* v4.7 */
     {
         BSIM4temp.T0 = model.BSIM4epsrox / 3.9;
         if ((model.BSIM4eotGiven) && (model.BSIM4toxpGiven) && (model.BSIM4dtoxGiven)
         && (std::abs(model.BSIM4eot *  BSIM4temp.T0 - (model.BSIM4toxp + model.BSIM4dtox)) > 1.0e-20))
         {
             printf("Warning: eot, toxp and dtox all given and eot * EPSROX / 3.9 != toxp + dtox; dtox ignored.\n");
         }
         else if ((model.BSIM4eotGiven) && (!model.BSIM4toxpGiven))
           model.BSIM4toxp =  BSIM4temp.T0 * model.BSIM4eot - model.BSIM4dtox;
         else if ((!model.BSIM4eotGiven) && (model.BSIM4toxpGiven)){
           model.BSIM4eot = (model.BSIM4toxp + model.BSIM4dtox) /  BSIM4temp.T0;
             if (!model.BSIM4toxmGiven)
                 model.BSIM4toxm = model.BSIM4eot;
         }
     }

     if(model.BSIM4mtrlMod)
       {
         BSIM4temp.epsrox = 3.9;
         BSIM4temp.toxe = model.BSIM4eot;
         BSIM4temp.epssub = EPS0 * model.BSIM4epsrsub;
       }
     else
       {
         BSIM4temp.epsrox = model.BSIM4epsrox;
         BSIM4temp.toxe = model.BSIM4toxe;
         BSIM4temp.epssub = EPSSI;
       }


         model.BSIM4coxe = BSIM4temp.epsrox * EPS0 / BSIM4temp.toxe;
     if(model.BSIM4mtrlMod == 0 || model.BSIM4mtrlCompatMod != 0)
       model.BSIM4coxp = model.BSIM4epsrox * EPS0 / model.BSIM4toxp;

     if (!model.BSIM4cgdoGiven)
         {   if (model.BSIM4dlcGiven && (model.BSIM4dlc > 0.0))
                 model.BSIM4cgdo = model.BSIM4dlc * model.BSIM4coxe
                                  - model.BSIM4cgdl ;
             else
                 model.BSIM4cgdo = 0.6 * model.BSIM4xj * model.BSIM4coxe;
         }
         if (!model.BSIM4cgsoGiven)
         {   if (model.BSIM4dlcGiven && (model.BSIM4dlc > 0.0))
                 model.BSIM4cgso = model.BSIM4dlc * model.BSIM4coxe
                                  - model.BSIM4cgsl ;
             else
                 model.BSIM4cgso = 0.6 * model.BSIM4xj * model.BSIM4coxe;
         }
         if (!model.BSIM4cgboGiven)
             model.BSIM4cgbo = 2.0 * model.BSIM4dwc * model.BSIM4coxe;
        model.vSizeDependParamKnot.clear();
        //  model.pSizeDependParamKnot = NULL;
        //  pLastKnot = NULL; Not pLastKnot in model any more. use std::vector to store the pSizeDependParamKnot

        BSIM4temp.Tnom = model.BSIM4tnom;
        BSIM4temp.TRatio = BSIM4temp.Temp / BSIM4temp.Tnom;

     model.BSIM4vcrit = CONSTvt0 * std::log(CONSTvt0 / (CONSTroot2 * 1.0e-14));
         model.BSIM4factor1 = std::sqrt(BSIM4temp.epssub / (BSIM4temp.epsrox * EPS0)* BSIM4temp.toxe);

         BSIM4temp.Vtm0 = model.BSIM4vtm0 = KboQ * BSIM4temp.Tnom;

     if(model.BSIM4mtrlMod==0)
     {
        BSIM4temp.Eg0 = 1.16 - 7.02e-4 * BSIM4temp.Tnom * BSIM4temp.Tnom / (BSIM4temp.Tnom + 1108.0);
        BSIM4temp.ni = 1.45e10 * (BSIM4temp.Tnom / 300.15) * std::sqrt(BSIM4temp.Tnom / 300.15)
                 * std::exp(21.5565981 - BSIM4temp.Eg0 / (2.0 * BSIM4temp.Vtm0));
     }
     else
     {
        BSIM4temp.Eg0 = model.BSIM4bg0sub - model.BSIM4tbgasub * BSIM4temp.Tnom * BSIM4temp.Tnom
                                      / (BSIM4temp.Tnom + model.BSIM4tbgbsub);
        BSIM4temp.T0 =  model.BSIM4bg0sub - model.BSIM4tbgasub * 90090.0225
                                      / (300.15 + model.BSIM4tbgbsub);
                                      BSIM4temp.ni = model.BSIM4ni0sub * (BSIM4temp.Tnom / 300.15) * std::sqrt(BSIM4temp.Tnom / 300.15)
                 * std::exp(( BSIM4temp.T0 - BSIM4temp.Eg0) / (2.0 * BSIM4temp.Vtm0));
     }

     model.BSIM4Eg0 = BSIM4temp.Eg0;
         model.BSIM4vtm = KboQ * BSIM4temp.Temp;
     if(model.BSIM4mtrlMod == 0)
     BSIM4temp.Eg = 1.16 - 7.02e-4 * BSIM4temp.Temp * BSIM4temp.Temp / (BSIM4temp.Temp + 1108.0);
     else
     BSIM4temp.Eg = model.BSIM4bg0sub - model.BSIM4tbgasub * BSIM4temp.Temp * BSIM4temp.Temp
                                      / (BSIM4temp.Temp + model.BSIM4tbgbsub);
     if (BSIM4temp.Temp != BSIM4temp.Tnom)
     {    BSIM4temp.T0 = BSIM4temp.Eg0 / BSIM4temp.Vtm0 - BSIM4temp.Eg / model.BSIM4vtm;
         BSIM4temp.T1 = std::log(BSIM4temp.Temp / BSIM4temp.Tnom);
         BSIM4temp.T2 =  BSIM4temp.T0 + model.BSIM4SjctTempExponent * BSIM4temp.T1;
         BSIM4temp.T3 = std::exp(BSIM4temp.T2 / model.BSIM4SjctEmissionCoeff);
         model.BSIM4SjctTempSatCurDensity = model.BSIM4SjctSatCurDensity
                           * BSIM4temp.T3;
         model.BSIM4SjctSidewallTempSatCurDensity
             = model.BSIM4SjctSidewallSatCurDensity * BSIM4temp.T3;
         model.BSIM4SjctGateSidewallTempSatCurDensity
                         = model.BSIM4SjctGateSidewallSatCurDensity * BSIM4temp.T3;

            BSIM4temp.T2 =  BSIM4temp.T0 + model.BSIM4DjctTempExponent * BSIM4temp.T1;
            BSIM4temp.T3 = std::exp(BSIM4temp.T2 / model.BSIM4DjctEmissionCoeff);
             model.BSIM4DjctTempSatCurDensity = model.BSIM4DjctSatCurDensity
                                               * BSIM4temp.T3;
             model.BSIM4DjctSidewallTempSatCurDensity
                         = model.BSIM4DjctSidewallSatCurDensity * BSIM4temp.T3;
             model.BSIM4DjctGateSidewallTempSatCurDensity
                         = model.BSIM4DjctGateSidewallSatCurDensity * BSIM4temp.T3;
     }
     else
     {   model.BSIM4SjctTempSatCurDensity = model.BSIM4SjctSatCurDensity;
         model.BSIM4SjctSidewallTempSatCurDensity
            = model.BSIM4SjctSidewallSatCurDensity;
             model.BSIM4SjctGateSidewallTempSatCurDensity
                        = model.BSIM4SjctGateSidewallSatCurDensity;
             model.BSIM4DjctTempSatCurDensity = model.BSIM4DjctSatCurDensity;
             model.BSIM4DjctSidewallTempSatCurDensity
                        = model.BSIM4DjctSidewallSatCurDensity;
             model.BSIM4DjctGateSidewallTempSatCurDensity
                        = model.BSIM4DjctGateSidewallSatCurDensity;
     }

     if (model.BSIM4SjctTempSatCurDensity < 0.0)
         model.BSIM4SjctTempSatCurDensity = 0.0;
     if (model.BSIM4SjctSidewallTempSatCurDensity < 0.0)
         model.BSIM4SjctSidewallTempSatCurDensity = 0.0;
         if (model.BSIM4SjctGateSidewallTempSatCurDensity < 0.0)
             model.BSIM4SjctGateSidewallTempSatCurDensity = 0.0;
         if (model.BSIM4DjctTempSatCurDensity < 0.0)
             model.BSIM4DjctTempSatCurDensity = 0.0;
         if (model.BSIM4DjctSidewallTempSatCurDensity < 0.0)
             model.BSIM4DjctSidewallTempSatCurDensity = 0.0;
         if (model.BSIM4DjctGateSidewallTempSatCurDensity < 0.0)
             model.BSIM4DjctGateSidewallTempSatCurDensity = 0.0;

     /* Temperature dependence of D/B and S/B diode capacitance begins */
     BSIM4temp.delTemp = CKTtemp - model.BSIM4tnom;
      BSIM4temp.T0 = model.BSIM4tcj * BSIM4temp.delTemp;
     if ( BSIM4temp.T0 >= -1.0)
     {   model.BSIM4SunitAreaTempJctCap = model.BSIM4SunitAreaJctCap *(1.0 +  BSIM4temp.T0); /*bug_fix -JX */
             model.BSIM4DunitAreaTempJctCap = model.BSIM4DunitAreaJctCap *(1.0 +  BSIM4temp.T0);
     }
     else
     {   if (model.BSIM4SunitAreaJctCap > 0.0)
         {   model.BSIM4SunitAreaTempJctCap = 0.0;
             fprintf(stderr, "Temperature effect has caused cjs to be negative. Cjs is clamped to zero.\n");
             }
         if (model.BSIM4DunitAreaJctCap > 0.0)
             {   model.BSIM4DunitAreaTempJctCap = 0.0;
                 fprintf(stderr, "Temperature effect has caused cjd to be negative. Cjd is clamped to zero.\n");
             }
     }
          BSIM4temp.T0 = model.BSIM4tcjsw * BSIM4temp.delTemp;
           if (model.BSIM4SunitLengthSidewallJctCap < 0.0)/*4.6.2*/
              {model.BSIM4SunitLengthSidewallJctCap = 0.0;
               fprintf(stderr, "CJSWS is negative. Cjsws is clamped to zero.\n");}
          if (model.BSIM4DunitLengthSidewallJctCap < 0.0)
              {model.BSIM4DunitLengthSidewallJctCap = 0.0;
               fprintf(stderr, "CJSWD is negative. Cjswd is clamped to zero.\n");}
     if ( BSIM4temp.T0 >= -1.0)
     {   model.BSIM4SunitLengthSidewallTempJctCap = model.BSIM4SunitLengthSidewallJctCap *(1.0 +  BSIM4temp.T0);
             model.BSIM4DunitLengthSidewallTempJctCap = model.BSIM4DunitLengthSidewallJctCap *(1.0 +  BSIM4temp.T0);
     }
     else
     {   if (model.BSIM4SunitLengthSidewallJctCap > 0.0)
         {   model.BSIM4SunitLengthSidewallTempJctCap = 0.0;
             fprintf(stderr, "Temperature effect has caused cjsws to be negative. Cjsws is clamped to zero.\n");
         }
         if (model.BSIM4DunitLengthSidewallJctCap > 0.0)
             {   model.BSIM4DunitLengthSidewallTempJctCap = 0.0;
                 fprintf(stderr, "Temperature effect has caused cjswd to be negative. Cjswd is clamped to zero.\n");
             }
     }
          BSIM4temp.T0 = model.BSIM4tcjswg * BSIM4temp.delTemp;
     if ( BSIM4temp.T0 >= -1.0)
     {   model.BSIM4SunitLengthGateSidewallTempJctCap = model.BSIM4SunitLengthGateSidewallJctCap *(1.0 +  BSIM4temp.T0);
             model.BSIM4DunitLengthGateSidewallTempJctCap = model.BSIM4DunitLengthGateSidewallJctCap *(1.0 +  BSIM4temp.T0);
     }
     else
     {   if (model.BSIM4SunitLengthGateSidewallJctCap > 0.0)
         {   model.BSIM4SunitLengthGateSidewallTempJctCap = 0.0;
             fprintf(stderr, "Temperature effect has caused cjswgs to be negative. Cjswgs is clamped to zero.\n");
         }
         if (model.BSIM4DunitLengthGateSidewallJctCap > 0.0)
             {   model.BSIM4DunitLengthGateSidewallTempJctCap = 0.0;
                 fprintf(stderr, "Temperature effect has caused cjswgd to be negative. Cjswgd is clamped to zero.\n");
             }
     }

         model.BSIM4PhiBS = model.BSIM4SbulkJctPotential
               - model.BSIM4tpb * BSIM4temp.delTemp;
         if (model.BSIM4PhiBS < 0.01)
     {   model.BSIM4PhiBS = 0.01;
         fprintf(stderr, "Temperature effect has caused pbs to be less than 0.01. Pbs is clamped to 0.01.\n");
     }
         model.BSIM4PhiBD = model.BSIM4DbulkJctPotential
                           - model.BSIM4tpb * BSIM4temp.delTemp;
         if (model.BSIM4PhiBD < 0.01)
         {   model.BSIM4PhiBD = 0.01;
             fprintf(stderr, "Temperature effect has caused pbd to be less than 0.01. Pbd is clamped to 0.01.\n");
         }

         model.BSIM4PhiBSWS = model.BSIM4SsidewallJctPotential
                             - model.BSIM4tpbsw * BSIM4temp.delTemp;
         if (model.BSIM4PhiBSWS <= 0.01)
     {   model.BSIM4PhiBSWS = 0.01;
         fprintf(stderr, "Temperature effect has caused pbsws to be less than 0.01. Pbsws is clamped to 0.01.\n");
     }
         model.BSIM4PhiBSWD = model.BSIM4DsidewallJctPotential
                             - model.BSIM4tpbsw * BSIM4temp.delTemp;
         if (model.BSIM4PhiBSWD <= 0.01)
         {   model.BSIM4PhiBSWD = 0.01;
             fprintf(stderr, "Temperature effect has caused pbswd to be less than 0.01. Pbswd is clamped to 0.01.\n");
         }

     model.BSIM4PhiBSWGS = model.BSIM4SGatesidewallJctPotential
                              - model.BSIM4tpbswg * BSIM4temp.delTemp;
         if (model.BSIM4PhiBSWGS <= 0.01)
     {   model.BSIM4PhiBSWGS = 0.01;
         fprintf(stderr, "Temperature effect has caused pbswgs to be less than 0.01. Pbswgs is clamped to 0.01.\n");
     }
         model.BSIM4PhiBSWGD = model.BSIM4DGatesidewallJctPotential
                              - model.BSIM4tpbswg * BSIM4temp.delTemp;
         if (model.BSIM4PhiBSWGD <= 0.01)
         {   model.BSIM4PhiBSWGD = 0.01;
             fprintf(stderr, "Temperature effect has caused pbswgd to be less than 0.01. Pbswgd is clamped to 0.01.\n");
         } /* End of junction capacitance */


         if (model.BSIM4ijthdfwd <= 0.0)
         {   model.BSIM4ijthdfwd = 0.0;
             fprintf(stderr, "Ijthdfwd reset to %g.\n", model.BSIM4ijthdfwd);
         }
         if (model.BSIM4ijthsfwd <= 0.0)
         {   model.BSIM4ijthsfwd = 0.0;
             fprintf(stderr, "Ijthsfwd reset to %g.\n", model.BSIM4ijthsfwd);
         }
     if (model.BSIM4ijthdrev <= 0.0)
         {   model.BSIM4ijthdrev = 0.0;
             fprintf(stderr, "Ijthdrev reset to %g.\n", model.BSIM4ijthdrev);
         }
         if (model.BSIM4ijthsrev <= 0.0)
         {   model.BSIM4ijthsrev = 0.0;
             fprintf(stderr, "Ijthsrev reset to %g.\n", model.BSIM4ijthsrev);
         }

         if ((model.BSIM4xjbvd <= 0.0) && (model.BSIM4dioMod == 2))
         {   model.BSIM4xjbvd = 0.0;
             fprintf(stderr, "Xjbvd reset to %g.\n", model.BSIM4xjbvd);
         }
         else if ((model.BSIM4xjbvd < 0.0) && (model.BSIM4dioMod == 0))
         {   model.BSIM4xjbvd = 0.0;
             fprintf(stderr, "Xjbvd reset to %g.\n", model.BSIM4xjbvd);
         }

         if (model.BSIM4bvd <= 0.0)   /*4.6.2*/
         {   model.BSIM4bvd = 0.0;
             fprintf(stderr, "BVD reset to %g.\n", model.BSIM4bvd);
         }

         if ((model.BSIM4xjbvs <= 0.0) && (model.BSIM4dioMod == 2))
         {   model.BSIM4xjbvs = 0.0;
             fprintf(stderr, "Xjbvs reset to %g.\n", model.BSIM4xjbvs);
         }
         else if ((model.BSIM4xjbvs < 0.0) && (model.BSIM4dioMod == 0))
         {   model.BSIM4xjbvs = 0.0;
             fprintf(stderr, "Xjbvs reset to %g.\n", model.BSIM4xjbvs);
         }

         if (model.BSIM4bvs <= 0.0)
         {   model.BSIM4bvs = 0.0;
             fprintf(stderr, "BVS reset to %g.\n", model.BSIM4bvs);
         }
    
    if(BSIM4checkModel(model, CKTtemp))
    {
        std::cerr << "BSIM4: model " << model.BSIM4modName << ": Model check failed." << std::endl;
        exit(1);
    }

    return 0;
}


int instanceTemp(BSIM4V82 &BSIM4instance, BSIM4model &model){
    int Size_Not_Found, i;
    // All Temporary parameters for this model
    bsim4v82temp &BSIM4temp = model.BSIM4temp;
    
    Size_Not_Found = 1;
    for(auto &pSizeDependParamKnot : model.vSizeDependParamKnot){
        if ((BSIM4instance.BSIM4l == pSizeDependParamKnot->Length)
        && (BSIM4instance.BSIM4w == pSizeDependParamKnot->Width)
        && (BSIM4instance.BSIM4nf == pSizeDependParamKnot->NFinger))
        {   Size_Not_Found = 0;
            BSIM4instance.pParam = pSizeDependParamKnot;
            BSIM4temp.pParam = BSIM4instance.pParam; /*bug-fix  */
        }
    }

    /* stress effect */
     BSIM4temp.Ldrn = BSIM4instance.BSIM4l;
     BSIM4temp.Wdrn = BSIM4instance.BSIM4w / BSIM4instance.BSIM4nf;

    if (Size_Not_Found)
    {   
        BSIM4temp.pParam = std::make_shared<bsim4SizeDependParam>();
        model.vSizeDependParamKnot.push_back(BSIM4temp.pParam);
        BSIM4instance.pParam = BSIM4temp.pParam;

        BSIM4temp.pParam->Length = BSIM4instance.BSIM4l;
        BSIM4temp.pParam->Width = BSIM4instance.BSIM4w;
        BSIM4temp.pParam->NFinger = BSIM4instance.BSIM4nf;
        BSIM4temp.Lnew = BSIM4instance.BSIM4l  + model.BSIM4xl ;
        BSIM4temp.Wnew = BSIM4instance.BSIM4w / BSIM4instance.BSIM4nf + model.BSIM4xw;

         BSIM4temp.T0 = std::pow(BSIM4temp.Lnew, model.BSIM4Lln);
        BSIM4temp.T1 = std::pow(BSIM4temp.Wnew, model.BSIM4Lwn);
         BSIM4temp.tmp1 = model.BSIM4Ll /  BSIM4temp.T0 + model.BSIM4Lw / BSIM4temp.T1
                       + model.BSIM4Lwl / ( BSIM4temp.T0 * BSIM4temp.T1);
        BSIM4temp.pParam->BSIM4dl = model.BSIM4Lint +  BSIM4temp.tmp1;
         BSIM4temp.tmp2 = model.BSIM4Llc /  BSIM4temp.T0 + model.BSIM4Lwc / BSIM4temp.T1
                       + model.BSIM4Lwlc / ( BSIM4temp.T0 * BSIM4temp.T1);
        BSIM4temp.pParam->BSIM4dlc = model.BSIM4dlc +  BSIM4temp.tmp2;

        BSIM4temp.T2 = std::pow(BSIM4temp.Lnew, model.BSIM4Wln);
        BSIM4temp.T3 = std::pow(BSIM4temp.Wnew, model.BSIM4Wwn);
         BSIM4temp.tmp1 = model.BSIM4Wl / BSIM4temp.T2 + model.BSIM4Ww / BSIM4temp.T3
                       + model.BSIM4Wwl / (BSIM4temp.T2 * BSIM4temp.T3);
        BSIM4temp.pParam->BSIM4dw = model.BSIM4Wint +  BSIM4temp.tmp1;
         BSIM4temp.tmp2 = model.BSIM4Wlc / BSIM4temp.T2 + model.BSIM4Wwc / BSIM4temp.T3
                       + model.BSIM4Wwlc / (BSIM4temp.T2 * BSIM4temp.T3);
        BSIM4temp.pParam->BSIM4dwc = model.BSIM4dwc +  BSIM4temp.tmp2;
        BSIM4temp.pParam->BSIM4dwj = model.BSIM4dwj +  BSIM4temp.tmp2;

        BSIM4temp.pParam->BSIM4leff = BSIM4temp.Lnew - 2.0 * BSIM4temp.pParam->BSIM4dl;
        if (BSIM4temp.pParam->BSIM4leff <= 0.0)
        {   
            std::cerr << "BSIM4: mosfet " << BSIM4instance.BSIM4name << ", model "<< model.BSIM4modName 
                    << ": Effective channel length <= 0" << std::endl;
            exit(1);
        }

        BSIM4temp.pParam->BSIM4weff = BSIM4temp.Wnew - 2.0 * BSIM4temp.pParam->BSIM4dw;
        if (BSIM4temp.pParam->BSIM4weff <= 0.0)
        {   
            std::cerr << "BSIM4: mosfet " << BSIM4instance.BSIM4name << ", model "<< model.BSIM4modName
            << ": Effective channel width <= 0" << std::endl;
            exit(1);
        }

        BSIM4temp.pParam->BSIM4leffCV = BSIM4temp.Lnew - 2.0 * BSIM4temp.pParam->BSIM4dlc;
        if (BSIM4temp.pParam->BSIM4leffCV <= 0.0)
        { 
            std::cerr << "BSIM4: mosfet " << BSIM4instance.BSIM4name << ", model "<< model.BSIM4modName
            << ": Effective channel length for C-V <= 0" << std::endl;
            exit(1);
        }

        BSIM4temp.pParam->BSIM4weffCV = BSIM4temp.Wnew - 2.0 * BSIM4temp.pParam->BSIM4dwc;
        if (BSIM4temp.pParam->BSIM4weffCV <= 0.0)
        {   
            std::cerr << "BSIM4: mosfet " << BSIM4instance.BSIM4name << ", model "<< model.BSIM4modName
            << ": Effective channel width for C-V <= 0" << std::endl;
            exit(1);
        }

        BSIM4temp.pParam->BSIM4weffCJ = BSIM4temp.Wnew - 2.0 * BSIM4temp.pParam->BSIM4dwj;
        if (BSIM4temp.pParam->BSIM4weffCJ <= 0.0)
        {   
            std::cerr << "BSIM4: mosfet " << BSIM4instance.BSIM4name << ", model "<< model.BSIM4modName
            << ": Effective channel width for S/D junctions <= 0" << std::endl;
            exit(1);
        }


        if (model.BSIM4binUnit == 1)
        {   BSIM4temp.Inv_L = 1.0e-6   / BSIM4temp.pParam->BSIM4leff;
            BSIM4temp.Inv_W = 1.0e-6   / BSIM4temp.pParam->BSIM4weff;
            BSIM4temp.Inv_LW = 1.0e-12 / (BSIM4temp.pParam->BSIM4leff * BSIM4temp.pParam->BSIM4weff);
        }
        else
        {   BSIM4temp.Inv_L = 1.0 / BSIM4temp.pParam->BSIM4leff;
            BSIM4temp.Inv_W = 1.0 / BSIM4temp.pParam->BSIM4weff;
            BSIM4temp.Inv_LW = 1.0 / (BSIM4temp.pParam->BSIM4leff * BSIM4temp.pParam->BSIM4weff);
        }
        BSIM4temp.pParam->BSIM4cdsc = model.BSIM4cdsc
                    + model.BSIM4lcdsc * BSIM4temp.Inv_L
                    + model.BSIM4wcdsc * BSIM4temp.Inv_W
                    + model.BSIM4pcdsc * BSIM4temp.Inv_LW;
        BSIM4temp.pParam->BSIM4cdscb = model.BSIM4cdscb
                     + model.BSIM4lcdscb * BSIM4temp.Inv_L
                     + model.BSIM4wcdscb * BSIM4temp.Inv_W
                     + model.BSIM4pcdscb * BSIM4temp.Inv_LW;

        BSIM4temp.pParam->BSIM4cdscd = model.BSIM4cdscd
                     + model.BSIM4lcdscd * BSIM4temp.Inv_L
                     + model.BSIM4wcdscd * BSIM4temp.Inv_W
                     + model.BSIM4pcdscd * BSIM4temp.Inv_LW;

        BSIM4temp.pParam->BSIM4cit = model.BSIM4cit
                   + model.BSIM4lcit * BSIM4temp.Inv_L
                   + model.BSIM4wcit * BSIM4temp.Inv_W
                   + model.BSIM4pcit * BSIM4temp.Inv_LW;
        BSIM4temp.pParam->BSIM4nfactor = model.BSIM4nfactor
                       + model.BSIM4lnfactor * BSIM4temp.Inv_L
                       + model.BSIM4wnfactor * BSIM4temp.Inv_W
                       + model.BSIM4pnfactor * BSIM4temp.Inv_LW;
        BSIM4temp.pParam->BSIM4tnfactor = model.BSIM4tnfactor          /* v4.7 */
                       + model.BSIM4ltnfactor * BSIM4temp.Inv_L
                       + model.BSIM4wtnfactor * BSIM4temp.Inv_W
                       + model.BSIM4ptnfactor * BSIM4temp.Inv_LW;
        BSIM4temp.pParam->BSIM4xj = model.BSIM4xj
                  + model.BSIM4lxj * BSIM4temp.Inv_L
                  + model.BSIM4wxj * BSIM4temp.Inv_W
                  + model.BSIM4pxj * BSIM4temp.Inv_LW;
        BSIM4temp.pParam->BSIM4vsat = model.BSIM4vsat
                    + model.BSIM4lvsat * BSIM4temp.Inv_L
                    + model.BSIM4wvsat * BSIM4temp.Inv_W
                    + model.BSIM4pvsat * BSIM4temp.Inv_LW;
        BSIM4temp.pParam->BSIM4at = model.BSIM4at
                  + model.BSIM4lat * BSIM4temp.Inv_L
                  + model.BSIM4wat * BSIM4temp.Inv_W
                  + model.BSIM4pat * BSIM4temp.Inv_LW;
        BSIM4temp.pParam->BSIM4a0 = model.BSIM4a0
                  + model.BSIM4la0 * BSIM4temp.Inv_L
                  + model.BSIM4wa0 * BSIM4temp.Inv_W
                  + model.BSIM4pa0 * BSIM4temp.Inv_LW;

       BSIM4temp.pParam->BSIM4ags = model.BSIM4ags
                  + model.BSIM4lags * BSIM4temp.Inv_L
                  + model.BSIM4wags * BSIM4temp.Inv_W
                  + model.BSIM4pags * BSIM4temp.Inv_LW;

       BSIM4temp.pParam->BSIM4a1 = model.BSIM4a1
                  + model.BSIM4la1 * BSIM4temp.Inv_L
                  + model.BSIM4wa1 * BSIM4temp.Inv_W
                  + model.BSIM4pa1 * BSIM4temp.Inv_LW;
        BSIM4temp.pParam->BSIM4a2 = model.BSIM4a2
                  + model.BSIM4la2 * BSIM4temp.Inv_L
                  + model.BSIM4wa2 * BSIM4temp.Inv_W
                  + model.BSIM4pa2 * BSIM4temp.Inv_LW;
       BSIM4temp.pParam->BSIM4keta = model.BSIM4keta
                    + model.BSIM4lketa * BSIM4temp.Inv_L
                    + model.BSIM4wketa * BSIM4temp.Inv_W
                    + model.BSIM4pketa * BSIM4temp.Inv_LW;
       BSIM4temp.pParam->BSIM4nsub = model.BSIM4nsub
                    + model.BSIM4lnsub * BSIM4temp.Inv_L
                    + model.BSIM4wnsub * BSIM4temp.Inv_W
                    + model.BSIM4pnsub * BSIM4temp.Inv_LW;
       BSIM4temp.pParam->BSIM4ndep = model.BSIM4ndep
                    + model.BSIM4lndep * BSIM4temp.Inv_L
                    + model.BSIM4wndep * BSIM4temp.Inv_W
                    + model.BSIM4pndep * BSIM4temp.Inv_LW;
       BSIM4temp.pParam->BSIM4nsd = model.BSIM4nsd
                                   + model.BSIM4lnsd * BSIM4temp.Inv_L
                                   + model.BSIM4wnsd * BSIM4temp.Inv_W
                                   + model.BSIM4pnsd * BSIM4temp.Inv_LW;
       BSIM4temp.pParam->BSIM4phin = model.BSIM4phin
                                    + model.BSIM4lphin * BSIM4temp.Inv_L
                                    + model.BSIM4wphin * BSIM4temp.Inv_W
                                    + model.BSIM4pphin * BSIM4temp.Inv_LW;
       BSIM4temp.pParam->BSIM4ngate = model.BSIM4ngate
                     + model.BSIM4lngate * BSIM4temp.Inv_L
                     + model.BSIM4wngate * BSIM4temp.Inv_W
                     + model.BSIM4pngate * BSIM4temp.Inv_LW;
       BSIM4temp.pParam->BSIM4gamma1 = model.BSIM4gamma1
                      + model.BSIM4lgamma1 * BSIM4temp.Inv_L
                      + model.BSIM4wgamma1 * BSIM4temp.Inv_W
                      + model.BSIM4pgamma1 * BSIM4temp.Inv_LW;
       BSIM4temp.pParam->BSIM4gamma2 = model.BSIM4gamma2
                      + model.BSIM4lgamma2 * BSIM4temp.Inv_L
                      + model.BSIM4wgamma2 * BSIM4temp.Inv_W
                      + model.BSIM4pgamma2 * BSIM4temp.Inv_LW;
         BSIM4temp.pParam->BSIM4vbx = model.BSIM4vbx
                   + model.BSIM4lvbx * BSIM4temp.Inv_L
                   + model.BSIM4wvbx * BSIM4temp.Inv_W
                   + model.BSIM4pvbx * BSIM4temp.Inv_LW;
         BSIM4temp.pParam->BSIM4vbm = model.BSIM4vbm
                   + model.BSIM4lvbm * BSIM4temp.Inv_L
                   + model.BSIM4wvbm * BSIM4temp.Inv_W
                   + model.BSIM4pvbm * BSIM4temp.Inv_LW;
         BSIM4temp.pParam->BSIM4xt = model.BSIM4xt
                   + model.BSIM4lxt * BSIM4temp.Inv_L
                   + model.BSIM4wxt * BSIM4temp.Inv_W
                   + model.BSIM4pxt * BSIM4temp.Inv_LW;
                 BSIM4temp.pParam->BSIM4vfb = model.BSIM4vfb
                                   + model.BSIM4lvfb * BSIM4temp.Inv_L
                                   + model.BSIM4wvfb * BSIM4temp.Inv_W
                                   + model.BSIM4pvfb * BSIM4temp.Inv_LW;
         BSIM4temp.pParam->BSIM4k1 = model.BSIM4k1
                  + model.BSIM4lk1 * BSIM4temp.Inv_L
                  + model.BSIM4wk1 * BSIM4temp.Inv_W
                  + model.BSIM4pk1 * BSIM4temp.Inv_LW;
         BSIM4temp.pParam->BSIM4kt1 = model.BSIM4kt1
                   + model.BSIM4lkt1 * BSIM4temp.Inv_L
                   + model.BSIM4wkt1 * BSIM4temp.Inv_W
                   + model.BSIM4pkt1 * BSIM4temp.Inv_LW;
         BSIM4temp.pParam->BSIM4kt1l = model.BSIM4kt1l
                    + model.BSIM4lkt1l * BSIM4temp.Inv_L
                    + model.BSIM4wkt1l * BSIM4temp.Inv_W
                    + model.BSIM4pkt1l * BSIM4temp.Inv_LW;
         BSIM4temp.pParam->BSIM4k2 = model.BSIM4k2
                  + model.BSIM4lk2 * BSIM4temp.Inv_L
                  + model.BSIM4wk2 * BSIM4temp.Inv_W
                  + model.BSIM4pk2 * BSIM4temp.Inv_LW;
         BSIM4temp.pParam->BSIM4kt2 = model.BSIM4kt2
                   + model.BSIM4lkt2 * BSIM4temp.Inv_L
                   + model.BSIM4wkt2 * BSIM4temp.Inv_W
                   + model.BSIM4pkt2 * BSIM4temp.Inv_LW;
         BSIM4temp.pParam->BSIM4k3 = model.BSIM4k3
                  + model.BSIM4lk3 * BSIM4temp.Inv_L
                  + model.BSIM4wk3 * BSIM4temp.Inv_W
                  + model.BSIM4pk3 * BSIM4temp.Inv_LW;
         BSIM4temp.pParam->BSIM4k3b = model.BSIM4k3b
                   + model.BSIM4lk3b * BSIM4temp.Inv_L
                   + model.BSIM4wk3b * BSIM4temp.Inv_W
                   + model.BSIM4pk3b * BSIM4temp.Inv_LW;
         BSIM4temp.pParam->BSIM4w0 = model.BSIM4w0
                  + model.BSIM4lw0 * BSIM4temp.Inv_L
                  + model.BSIM4ww0 * BSIM4temp.Inv_W
                  + model.BSIM4pw0 * BSIM4temp.Inv_LW;
         BSIM4temp.pParam->BSIM4lpe0 = model.BSIM4lpe0
                    + model.BSIM4llpe0 * BSIM4temp.Inv_L
                    + model.BSIM4wlpe0 * BSIM4temp.Inv_W
                    + model.BSIM4plpe0 * BSIM4temp.Inv_LW;
                 BSIM4temp.pParam->BSIM4lpeb = model.BSIM4lpeb
                                    + model.BSIM4llpeb * BSIM4temp.Inv_L
                                    + model.BSIM4wlpeb * BSIM4temp.Inv_W
                                    + model.BSIM4plpeb * BSIM4temp.Inv_LW;
                 BSIM4temp.pParam->BSIM4dvtp0 = model.BSIM4dvtp0
                                     + model.BSIM4ldvtp0 * BSIM4temp.Inv_L
                                     + model.BSIM4wdvtp0 * BSIM4temp.Inv_W
                                     + model.BSIM4pdvtp0 * BSIM4temp.Inv_LW;
                 BSIM4temp.pParam->BSIM4dvtp1 = model.BSIM4dvtp1
                                     + model.BSIM4ldvtp1 * BSIM4temp.Inv_L
                                     + model.BSIM4wdvtp1 * BSIM4temp.Inv_W
                                     + model.BSIM4pdvtp1 * BSIM4temp.Inv_LW;
                 BSIM4temp.pParam->BSIM4dvtp2 = model.BSIM4dvtp2        /* v4.7  */
                                     + model.BSIM4ldvtp2 * BSIM4temp.Inv_L
                                     + model.BSIM4wdvtp2 * BSIM4temp.Inv_W
                                     + model.BSIM4pdvtp2 * BSIM4temp.Inv_LW;
                 BSIM4temp.pParam->BSIM4dvtp3 = model.BSIM4dvtp3        /* v4.7  */
                                     + model.BSIM4ldvtp3 * BSIM4temp.Inv_L
                                     + model.BSIM4wdvtp3 * BSIM4temp.Inv_W
                                     + model.BSIM4pdvtp3 * BSIM4temp.Inv_LW;
                 BSIM4temp.pParam->BSIM4dvtp4 = model.BSIM4dvtp4        /* v4.7  */
                                     + model.BSIM4ldvtp4 * BSIM4temp.Inv_L
                                     + model.BSIM4wdvtp4 * BSIM4temp.Inv_W
                                     + model.BSIM4pdvtp4 * BSIM4temp.Inv_LW;
                 BSIM4temp.pParam->BSIM4dvtp5 = model.BSIM4dvtp5        /* v4.7  */
                                     + model.BSIM4ldvtp5 * BSIM4temp.Inv_L
                                     + model.BSIM4wdvtp5 * BSIM4temp.Inv_W
                                     + model.BSIM4pdvtp5 * BSIM4temp.Inv_LW;
         BSIM4temp.pParam->BSIM4dvt0 = model.BSIM4dvt0
                    + model.BSIM4ldvt0 * BSIM4temp.Inv_L
                    + model.BSIM4wdvt0 * BSIM4temp.Inv_W
                    + model.BSIM4pdvt0 * BSIM4temp.Inv_LW;
         BSIM4temp.pParam->BSIM4dvt1 = model.BSIM4dvt1
                    + model.BSIM4ldvt1 * BSIM4temp.Inv_L
                    + model.BSIM4wdvt1 * BSIM4temp.Inv_W
                    + model.BSIM4pdvt1 * BSIM4temp.Inv_LW;
         BSIM4temp.pParam->BSIM4dvt2 = model.BSIM4dvt2
                    + model.BSIM4ldvt2 * BSIM4temp.Inv_L
                    + model.BSIM4wdvt2 * BSIM4temp.Inv_W
                    + model.BSIM4pdvt2 * BSIM4temp.Inv_LW;
         BSIM4temp.pParam->BSIM4dvt0w = model.BSIM4dvt0w
                    + model.BSIM4ldvt0w * BSIM4temp.Inv_L
                    + model.BSIM4wdvt0w * BSIM4temp.Inv_W
                    + model.BSIM4pdvt0w * BSIM4temp.Inv_LW;
         BSIM4temp.pParam->BSIM4dvt1w = model.BSIM4dvt1w
                    + model.BSIM4ldvt1w * BSIM4temp.Inv_L
                    + model.BSIM4wdvt1w * BSIM4temp.Inv_W
                    + model.BSIM4pdvt1w * BSIM4temp.Inv_LW;
         BSIM4temp.pParam->BSIM4dvt2w = model.BSIM4dvt2w
                    + model.BSIM4ldvt2w * BSIM4temp.Inv_L
                    + model.BSIM4wdvt2w * BSIM4temp.Inv_W
                    + model.BSIM4pdvt2w * BSIM4temp.Inv_LW;
         BSIM4temp.pParam->BSIM4drout = model.BSIM4drout
                     + model.BSIM4ldrout * BSIM4temp.Inv_L
                     + model.BSIM4wdrout * BSIM4temp.Inv_W
                     + model.BSIM4pdrout * BSIM4temp.Inv_LW;
         BSIM4temp.pParam->BSIM4dsub = model.BSIM4dsub
                    + model.BSIM4ldsub * BSIM4temp.Inv_L
                    + model.BSIM4wdsub * BSIM4temp.Inv_W
                    + model.BSIM4pdsub * BSIM4temp.Inv_LW;
         BSIM4temp.pParam->BSIM4vth0 = model.BSIM4vth0
                    + model.BSIM4lvth0 * BSIM4temp.Inv_L
                    + model.BSIM4wvth0 * BSIM4temp.Inv_W
                    + model.BSIM4pvth0 * BSIM4temp.Inv_LW;
         BSIM4temp.pParam->BSIM4ua = model.BSIM4ua
                  + model.BSIM4lua * BSIM4temp.Inv_L
                  + model.BSIM4wua * BSIM4temp.Inv_W
                  + model.BSIM4pua * BSIM4temp.Inv_LW;
         BSIM4temp.pParam->BSIM4ua1 = model.BSIM4ua1
                   + model.BSIM4lua1 * BSIM4temp.Inv_L
                   + model.BSIM4wua1 * BSIM4temp.Inv_W
                   + model.BSIM4pua1 * BSIM4temp.Inv_LW;
         BSIM4temp.pParam->BSIM4ub = model.BSIM4ub
                  + model.BSIM4lub * BSIM4temp.Inv_L
                  + model.BSIM4wub * BSIM4temp.Inv_W
                  + model.BSIM4pub * BSIM4temp.Inv_LW;
         BSIM4temp.pParam->BSIM4ub1 = model.BSIM4ub1
                   + model.BSIM4lub1 * BSIM4temp.Inv_L
                   + model.BSIM4wub1 * BSIM4temp.Inv_W
                   + model.BSIM4pub1 * BSIM4temp.Inv_LW;
         BSIM4temp.pParam->BSIM4uc = model.BSIM4uc
                  + model.BSIM4luc * BSIM4temp.Inv_L
                  + model.BSIM4wuc * BSIM4temp.Inv_W
                  + model.BSIM4puc * BSIM4temp.Inv_LW;
         BSIM4temp.pParam->BSIM4uc1 = model.BSIM4uc1
                   + model.BSIM4luc1 * BSIM4temp.Inv_L
                   + model.BSIM4wuc1 * BSIM4temp.Inv_W
                   + model.BSIM4puc1 * BSIM4temp.Inv_LW;
         BSIM4temp.pParam->BSIM4ud = model.BSIM4ud
                  + model.BSIM4lud * BSIM4temp.Inv_L
                  + model.BSIM4wud * BSIM4temp.Inv_W
                  + model.BSIM4pud * BSIM4temp.Inv_LW;
         BSIM4temp.pParam->BSIM4ud1 = model.BSIM4ud1
                  + model.BSIM4lud1 * BSIM4temp.Inv_L
                  + model.BSIM4wud1 * BSIM4temp.Inv_W
                  + model.BSIM4pud1 * BSIM4temp.Inv_LW;
         BSIM4temp.pParam->BSIM4up = model.BSIM4up
                  + model.BSIM4lup * BSIM4temp.Inv_L
                  + model.BSIM4wup * BSIM4temp.Inv_W
                  + model.BSIM4pup * BSIM4temp.Inv_LW;
         BSIM4temp.pParam->BSIM4lp = model.BSIM4lp
                  + model.BSIM4llp * BSIM4temp.Inv_L
                  + model.BSIM4wlp * BSIM4temp.Inv_W
                  + model.BSIM4plp * BSIM4temp.Inv_LW;
                 BSIM4temp.pParam->BSIM4eu = model.BSIM4eu
                                  + model.BSIM4leu * BSIM4temp.Inv_L
                                  + model.BSIM4weu * BSIM4temp.Inv_W
                                  + model.BSIM4peu * BSIM4temp.Inv_LW;
         BSIM4temp.pParam->BSIM4u0 = model.BSIM4u0
                  + model.BSIM4lu0 * BSIM4temp.Inv_L
                  + model.BSIM4wu0 * BSIM4temp.Inv_W
                  + model.BSIM4pu0 * BSIM4temp.Inv_LW;
         BSIM4temp.pParam->BSIM4ute = model.BSIM4ute
                   + model.BSIM4lute * BSIM4temp.Inv_L
                   + model.BSIM4wute * BSIM4temp.Inv_W
                   + model.BSIM4pute * BSIM4temp.Inv_LW;
        /*high k mobility*/
        BSIM4temp.pParam->BSIM4ucs = model.BSIM4ucs
                  + model.BSIM4lucs * BSIM4temp.Inv_L
                  + model.BSIM4wucs * BSIM4temp.Inv_W
                  + model.BSIM4pucs * BSIM4temp.Inv_LW;
         BSIM4temp.pParam->BSIM4ucste = model.BSIM4ucste
                   + model.BSIM4lucste * BSIM4temp.Inv_L
                   + model.BSIM4wucste * BSIM4temp.Inv_W
                   + model.BSIM4pucste * BSIM4temp.Inv_LW;

         BSIM4temp.pParam->BSIM4voff = model.BSIM4voff
                    + model.BSIM4lvoff * BSIM4temp.Inv_L
                    + model.BSIM4wvoff * BSIM4temp.Inv_W
                    + model.BSIM4pvoff * BSIM4temp.Inv_LW;
         BSIM4temp.pParam->BSIM4tvoff = model.BSIM4tvoff
                    + model.BSIM4ltvoff * BSIM4temp.Inv_L
                    + model.BSIM4wtvoff * BSIM4temp.Inv_W
                    + model.BSIM4ptvoff * BSIM4temp.Inv_LW;
                 BSIM4temp.pParam->BSIM4minv = model.BSIM4minv
                                    + model.BSIM4lminv * BSIM4temp.Inv_L
                                    + model.BSIM4wminv * BSIM4temp.Inv_W
                                    + model.BSIM4pminv * BSIM4temp.Inv_LW;
                 BSIM4temp.pParam->BSIM4minvcv = model.BSIM4minvcv
                                    + model.BSIM4lminvcv * BSIM4temp.Inv_L
                                    + model.BSIM4wminvcv * BSIM4temp.Inv_W
                                    + model.BSIM4pminvcv * BSIM4temp.Inv_LW;
         BSIM4temp.pParam->BSIM4fprout = model.BSIM4fprout
                                     + model.BSIM4lfprout * BSIM4temp.Inv_L
                                     + model.BSIM4wfprout * BSIM4temp.Inv_W
                                     + model.BSIM4pfprout * BSIM4temp.Inv_LW;
                 BSIM4temp.pParam->BSIM4pdits = model.BSIM4pdits
                                     + model.BSIM4lpdits * BSIM4temp.Inv_L
                                     + model.BSIM4wpdits * BSIM4temp.Inv_W
                                     + model.BSIM4ppdits * BSIM4temp.Inv_LW;
                 BSIM4temp.pParam->BSIM4pditsd = model.BSIM4pditsd
                                      + model.BSIM4lpditsd * BSIM4temp.Inv_L
                                      + model.BSIM4wpditsd * BSIM4temp.Inv_W
                                      + model.BSIM4ppditsd * BSIM4temp.Inv_LW;
         BSIM4temp.pParam->BSIM4delta = model.BSIM4delta
                     + model.BSIM4ldelta * BSIM4temp.Inv_L
                     + model.BSIM4wdelta * BSIM4temp.Inv_W
                     + model.BSIM4pdelta * BSIM4temp.Inv_LW;
         BSIM4temp.pParam->BSIM4rdsw = model.BSIM4rdsw
                    + model.BSIM4lrdsw * BSIM4temp.Inv_L
                    + model.BSIM4wrdsw * BSIM4temp.Inv_W
                    + model.BSIM4prdsw * BSIM4temp.Inv_LW;
                 BSIM4temp.pParam->BSIM4rdw = model.BSIM4rdw
                                    + model.BSIM4lrdw * BSIM4temp.Inv_L
                                    + model.BSIM4wrdw * BSIM4temp.Inv_W
                                    + model.BSIM4prdw * BSIM4temp.Inv_LW;
                 BSIM4temp.pParam->BSIM4rsw = model.BSIM4rsw
                                    + model.BSIM4lrsw * BSIM4temp.Inv_L
                                    + model.BSIM4wrsw * BSIM4temp.Inv_W
                                    + model.BSIM4prsw * BSIM4temp.Inv_LW;
         BSIM4temp.pParam->BSIM4prwg = model.BSIM4prwg
                    + model.BSIM4lprwg * BSIM4temp.Inv_L
                    + model.BSIM4wprwg * BSIM4temp.Inv_W
                    + model.BSIM4pprwg * BSIM4temp.Inv_LW;
         BSIM4temp.pParam->BSIM4prwb = model.BSIM4prwb
                    + model.BSIM4lprwb * BSIM4temp.Inv_L
                    + model.BSIM4wprwb * BSIM4temp.Inv_W
                    + model.BSIM4pprwb * BSIM4temp.Inv_LW;
         BSIM4temp.pParam->BSIM4prt = model.BSIM4prt
                    + model.BSIM4lprt * BSIM4temp.Inv_L
                    + model.BSIM4wprt * BSIM4temp.Inv_W
                    + model.BSIM4pprt * BSIM4temp.Inv_LW;
         BSIM4temp.pParam->BSIM4eta0 = model.BSIM4eta0
                    + model.BSIM4leta0 * BSIM4temp.Inv_L
                    + model.BSIM4weta0 * BSIM4temp.Inv_W
                    + model.BSIM4peta0 * BSIM4temp.Inv_LW;
         BSIM4temp.pParam->BSIM4teta0 = model.BSIM4teta0        /* v4.7  */
                    + model.BSIM4lteta0 * BSIM4temp.Inv_L
                    + model.BSIM4wteta0 * BSIM4temp.Inv_W
                    + model.BSIM4pteta0 * BSIM4temp.Inv_LW;
         BSIM4temp.pParam->BSIM4tvoffcv = model.BSIM4tvoffcv    /* v4.8.0  */
                    + model.BSIM4ltvoffcv * BSIM4temp.Inv_L
                    + model.BSIM4wtvoffcv * BSIM4temp.Inv_W
                    + model.BSIM4ptvoffcv * BSIM4temp.Inv_LW;
         BSIM4temp.pParam->BSIM4etab = model.BSIM4etab
                    + model.BSIM4letab * BSIM4temp.Inv_L
                    + model.BSIM4wetab * BSIM4temp.Inv_W
                    + model.BSIM4petab * BSIM4temp.Inv_LW;
         BSIM4temp.pParam->BSIM4pclm = model.BSIM4pclm
                    + model.BSIM4lpclm * BSIM4temp.Inv_L
                    + model.BSIM4wpclm * BSIM4temp.Inv_W
                    + model.BSIM4ppclm * BSIM4temp.Inv_LW;
         BSIM4temp.pParam->BSIM4pdibl1 = model.BSIM4pdibl1
                      + model.BSIM4lpdibl1 * BSIM4temp.Inv_L
                      + model.BSIM4wpdibl1 * BSIM4temp.Inv_W
                      + model.BSIM4ppdibl1 * BSIM4temp.Inv_LW;
         BSIM4temp.pParam->BSIM4pdibl2 = model.BSIM4pdibl2
                      + model.BSIM4lpdibl2 * BSIM4temp.Inv_L
                      + model.BSIM4wpdibl2 * BSIM4temp.Inv_W
                      + model.BSIM4ppdibl2 * BSIM4temp.Inv_LW;
         BSIM4temp.pParam->BSIM4pdiblb = model.BSIM4pdiblb
                      + model.BSIM4lpdiblb * BSIM4temp.Inv_L
                      + model.BSIM4wpdiblb * BSIM4temp.Inv_W
                      + model.BSIM4ppdiblb * BSIM4temp.Inv_LW;
         BSIM4temp.pParam->BSIM4pscbe1 = model.BSIM4pscbe1
                      + model.BSIM4lpscbe1 * BSIM4temp.Inv_L
                      + model.BSIM4wpscbe1 * BSIM4temp.Inv_W
                      + model.BSIM4ppscbe1 * BSIM4temp.Inv_LW;
         BSIM4temp.pParam->BSIM4pscbe2 = model.BSIM4pscbe2
                      + model.BSIM4lpscbe2 * BSIM4temp.Inv_L
                      + model.BSIM4wpscbe2 * BSIM4temp.Inv_W
                      + model.BSIM4ppscbe2 * BSIM4temp.Inv_LW;
         BSIM4temp.pParam->BSIM4pvag = model.BSIM4pvag
                    + model.BSIM4lpvag * BSIM4temp.Inv_L
                    + model.BSIM4wpvag * BSIM4temp.Inv_W
                    + model.BSIM4ppvag * BSIM4temp.Inv_LW;
         BSIM4temp.pParam->BSIM4wr = model.BSIM4wr
                  + model.BSIM4lwr * BSIM4temp.Inv_L
                  + model.BSIM4wwr * BSIM4temp.Inv_W
                  + model.BSIM4pwr * BSIM4temp.Inv_LW;
         BSIM4temp.pParam->BSIM4dwg = model.BSIM4dwg
                   + model.BSIM4ldwg * BSIM4temp.Inv_L
                   + model.BSIM4wdwg * BSIM4temp.Inv_W
                   + model.BSIM4pdwg * BSIM4temp.Inv_LW;
         BSIM4temp.pParam->BSIM4dwb = model.BSIM4dwb
                   + model.BSIM4ldwb * BSIM4temp.Inv_L
                   + model.BSIM4wdwb * BSIM4temp.Inv_W
                   + model.BSIM4pdwb * BSIM4temp.Inv_LW;
         BSIM4temp.pParam->BSIM4b0 = model.BSIM4b0
                  + model.BSIM4lb0 * BSIM4temp.Inv_L
                  + model.BSIM4wb0 * BSIM4temp.Inv_W
                  + model.BSIM4pb0 * BSIM4temp.Inv_LW;
         BSIM4temp.pParam->BSIM4b1 = model.BSIM4b1
                  + model.BSIM4lb1 * BSIM4temp.Inv_L
                  + model.BSIM4wb1 * BSIM4temp.Inv_W
                  + model.BSIM4pb1 * BSIM4temp.Inv_LW;
         BSIM4temp.pParam->BSIM4alpha0 = model.BSIM4alpha0
                      + model.BSIM4lalpha0 * BSIM4temp.Inv_L
                      + model.BSIM4walpha0 * BSIM4temp.Inv_W
                      + model.BSIM4palpha0 * BSIM4temp.Inv_LW;
                 BSIM4temp.pParam->BSIM4alpha1 = model.BSIM4alpha1
                                      + model.BSIM4lalpha1 * BSIM4temp.Inv_L
                                      + model.BSIM4walpha1 * BSIM4temp.Inv_W
                                      + model.BSIM4palpha1 * BSIM4temp.Inv_LW;
         BSIM4temp.pParam->BSIM4beta0 = model.BSIM4beta0
                     + model.BSIM4lbeta0 * BSIM4temp.Inv_L
                     + model.BSIM4wbeta0 * BSIM4temp.Inv_W
                     + model.BSIM4pbeta0 * BSIM4temp.Inv_LW;
                 BSIM4temp.pParam->BSIM4agidl = model.BSIM4agidl
                                     + model.BSIM4lagidl * BSIM4temp.Inv_L
                                     + model.BSIM4wagidl * BSIM4temp.Inv_W
                                     + model.BSIM4pagidl * BSIM4temp.Inv_LW;
                 BSIM4temp.pParam->BSIM4bgidl = model.BSIM4bgidl
                                     + model.BSIM4lbgidl * BSIM4temp.Inv_L
                                     + model.BSIM4wbgidl * BSIM4temp.Inv_W
                                     + model.BSIM4pbgidl * BSIM4temp.Inv_LW;
                 BSIM4temp.pParam->BSIM4cgidl = model.BSIM4cgidl
                                     + model.BSIM4lcgidl * BSIM4temp.Inv_L
                                     + model.BSIM4wcgidl * BSIM4temp.Inv_W
                                     + model.BSIM4pcgidl * BSIM4temp.Inv_LW;
                 BSIM4temp.pParam->BSIM4egidl = model.BSIM4egidl
                                     + model.BSIM4legidl * BSIM4temp.Inv_L
                                     + model.BSIM4wegidl * BSIM4temp.Inv_W
                                     + model.BSIM4pegidl * BSIM4temp.Inv_LW;
         BSIM4temp.pParam->BSIM4rgidl = model.BSIM4rgidl        /* v4.7 New GIDL/GISL */
                                     + model.BSIM4lrgidl * BSIM4temp.Inv_L
                                     + model.BSIM4wrgidl * BSIM4temp.Inv_W
                                     + model.BSIM4prgidl * BSIM4temp.Inv_LW;
                 BSIM4temp.pParam->BSIM4kgidl = model.BSIM4kgidl        /* v4.7 New GIDL/GISL */
                                     + model.BSIM4lkgidl * BSIM4temp.Inv_L
                                     + model.BSIM4wkgidl * BSIM4temp.Inv_W
                                     + model.BSIM4pkgidl * BSIM4temp.Inv_LW;
                 BSIM4temp.pParam->BSIM4fgidl = model.BSIM4fgidl        /* v4.7 New GIDL/GISL */
                                     + model.BSIM4lfgidl * BSIM4temp.Inv_L
                                     + model.BSIM4wfgidl * BSIM4temp.Inv_W
                                     + model.BSIM4pfgidl * BSIM4temp.Inv_LW;
                 BSIM4temp.pParam->BSIM4agisl = model.BSIM4agisl
                                     + model.BSIM4lagisl * BSIM4temp.Inv_L
                                     + model.BSIM4wagisl * BSIM4temp.Inv_W
                                     + model.BSIM4pagisl * BSIM4temp.Inv_LW;
                 BSIM4temp.pParam->BSIM4bgisl = model.BSIM4bgisl
                                     + model.BSIM4lbgisl * BSIM4temp.Inv_L
                                     + model.BSIM4wbgisl * BSIM4temp.Inv_W
                                     + model.BSIM4pbgisl * BSIM4temp.Inv_LW;
                 BSIM4temp.pParam->BSIM4cgisl = model.BSIM4cgisl
                                     + model.BSIM4lcgisl * BSIM4temp.Inv_L
                                     + model.BSIM4wcgisl * BSIM4temp.Inv_W
                                     + model.BSIM4pcgisl * BSIM4temp.Inv_LW;
                 BSIM4temp.pParam->BSIM4egisl = model.BSIM4egisl
                                     + model.BSIM4legisl * BSIM4temp.Inv_L
                                     + model.BSIM4wegisl * BSIM4temp.Inv_W
                                     + model.BSIM4pegisl * BSIM4temp.Inv_LW;
         BSIM4temp.pParam->BSIM4rgisl = model.BSIM4rgisl        /* v4.7 New GIDL/GISL */
                                     + model.BSIM4lrgisl * BSIM4temp.Inv_L
                                     + model.BSIM4wrgisl * BSIM4temp.Inv_W
                                     + model.BSIM4prgisl * BSIM4temp.Inv_LW;
                 BSIM4temp.pParam->BSIM4kgisl = model.BSIM4kgisl        /* v4.7 New GIDL/GISL */
                                     + model.BSIM4lkgisl * BSIM4temp.Inv_L
                                     + model.BSIM4wkgisl * BSIM4temp.Inv_W
                                     + model.BSIM4pkgisl * BSIM4temp.Inv_LW;
                 BSIM4temp.pParam->BSIM4fgisl = model.BSIM4fgisl        /* v4.7 New GIDL/GISL */
                                     + model.BSIM4lfgisl * BSIM4temp.Inv_L
                                     + model.BSIM4wfgisl * BSIM4temp.Inv_W
                                     + model.BSIM4pfgisl * BSIM4temp.Inv_LW;
                 BSIM4temp.pParam->BSIM4aigc = model.BSIM4aigc
                                     + model.BSIM4laigc * BSIM4temp.Inv_L
                                     + model.BSIM4waigc * BSIM4temp.Inv_W
                                     + model.BSIM4paigc * BSIM4temp.Inv_LW;
                 BSIM4temp.pParam->BSIM4bigc = model.BSIM4bigc
                                     + model.BSIM4lbigc * BSIM4temp.Inv_L
                                     + model.BSIM4wbigc * BSIM4temp.Inv_W
                                     + model.BSIM4pbigc * BSIM4temp.Inv_LW;
                 BSIM4temp.pParam->BSIM4cigc = model.BSIM4cigc
                                     + model.BSIM4lcigc * BSIM4temp.Inv_L
                                     + model.BSIM4wcigc * BSIM4temp.Inv_W
                                     + model.BSIM4pcigc * BSIM4temp.Inv_LW;
                 BSIM4temp.pParam->BSIM4aigsd = model.BSIM4aigsd
                                     + model.BSIM4laigsd * BSIM4temp.Inv_L
                                     + model.BSIM4waigsd * BSIM4temp.Inv_W
                                     + model.BSIM4paigsd * BSIM4temp.Inv_LW;
                 BSIM4temp.pParam->BSIM4bigsd = model.BSIM4bigsd
                                     + model.BSIM4lbigsd * BSIM4temp.Inv_L
                                     + model.BSIM4wbigsd * BSIM4temp.Inv_W
                                     + model.BSIM4pbigsd * BSIM4temp.Inv_LW;
                 BSIM4temp.pParam->BSIM4cigsd = model.BSIM4cigsd
                                     + model.BSIM4lcigsd * BSIM4temp.Inv_L
                                     + model.BSIM4wcigsd * BSIM4temp.Inv_W
                                     + model.BSIM4pcigsd * BSIM4temp.Inv_LW;
                 BSIM4temp.pParam->BSIM4aigs = model.BSIM4aigs
                                     + model.BSIM4laigs * BSIM4temp.Inv_L
                                     + model.BSIM4waigs * BSIM4temp.Inv_W
                                     + model.BSIM4paigs * BSIM4temp.Inv_LW;
                 BSIM4temp.pParam->BSIM4bigs = model.BSIM4bigs
                                     + model.BSIM4lbigs * BSIM4temp.Inv_L
                                     + model.BSIM4wbigs * BSIM4temp.Inv_W
                                     + model.BSIM4pbigs * BSIM4temp.Inv_LW;
                 BSIM4temp.pParam->BSIM4cigs = model.BSIM4cigs
                                     + model.BSIM4lcigs * BSIM4temp.Inv_L
                                     + model.BSIM4wcigs * BSIM4temp.Inv_W
                                     + model.BSIM4pcigs * BSIM4temp.Inv_LW;
                 BSIM4temp.pParam->BSIM4aigd = model.BSIM4aigd
                                     + model.BSIM4laigd * BSIM4temp.Inv_L
                                     + model.BSIM4waigd * BSIM4temp.Inv_W
                                     + model.BSIM4paigd * BSIM4temp.Inv_LW;
                 BSIM4temp.pParam->BSIM4bigd = model.BSIM4bigd
                                     + model.BSIM4lbigd * BSIM4temp.Inv_L
                                     + model.BSIM4wbigd * BSIM4temp.Inv_W
                                     + model.BSIM4pbigd * BSIM4temp.Inv_LW;
                 BSIM4temp.pParam->BSIM4cigd = model.BSIM4cigd
                                     + model.BSIM4lcigd * BSIM4temp.Inv_L
                                     + model.BSIM4wcigd * BSIM4temp.Inv_W
                                     + model.BSIM4pcigd * BSIM4temp.Inv_LW;
                 BSIM4temp.pParam->BSIM4aigbacc = model.BSIM4aigbacc
                                       + model.BSIM4laigbacc * BSIM4temp.Inv_L
                                       + model.BSIM4waigbacc * BSIM4temp.Inv_W
                                       + model.BSIM4paigbacc * BSIM4temp.Inv_LW;
                 BSIM4temp.pParam->BSIM4bigbacc = model.BSIM4bigbacc
                                       + model.BSIM4lbigbacc * BSIM4temp.Inv_L
                                       + model.BSIM4wbigbacc * BSIM4temp.Inv_W
                                       + model.BSIM4pbigbacc * BSIM4temp.Inv_LW;
                 BSIM4temp.pParam->BSIM4cigbacc = model.BSIM4cigbacc
                                       + model.BSIM4lcigbacc * BSIM4temp.Inv_L
                                       + model.BSIM4wcigbacc * BSIM4temp.Inv_W
                                       + model.BSIM4pcigbacc * BSIM4temp.Inv_LW;
                 BSIM4temp.pParam->BSIM4aigbinv = model.BSIM4aigbinv
                                       + model.BSIM4laigbinv * BSIM4temp.Inv_L
                                       + model.BSIM4waigbinv * BSIM4temp.Inv_W
                                       + model.BSIM4paigbinv * BSIM4temp.Inv_LW;
                 BSIM4temp.pParam->BSIM4bigbinv = model.BSIM4bigbinv
                                       + model.BSIM4lbigbinv * BSIM4temp.Inv_L
                                       + model.BSIM4wbigbinv * BSIM4temp.Inv_W
                                       + model.BSIM4pbigbinv * BSIM4temp.Inv_LW;
                 BSIM4temp.pParam->BSIM4cigbinv = model.BSIM4cigbinv
                                       + model.BSIM4lcigbinv * BSIM4temp.Inv_L
                                       + model.BSIM4wcigbinv * BSIM4temp.Inv_W
                                       + model.BSIM4pcigbinv * BSIM4temp.Inv_LW;
                 BSIM4temp.pParam->BSIM4nigc = model.BSIM4nigc
                                       + model.BSIM4lnigc * BSIM4temp.Inv_L
                                       + model.BSIM4wnigc * BSIM4temp.Inv_W
                                       + model.BSIM4pnigc * BSIM4temp.Inv_LW;
                 BSIM4temp.pParam->BSIM4nigbacc = model.BSIM4nigbacc
                                       + model.BSIM4lnigbacc * BSIM4temp.Inv_L
                                       + model.BSIM4wnigbacc * BSIM4temp.Inv_W
                                       + model.BSIM4pnigbacc * BSIM4temp.Inv_LW;
                 BSIM4temp.pParam->BSIM4nigbinv = model.BSIM4nigbinv
                                       + model.BSIM4lnigbinv * BSIM4temp.Inv_L
                                       + model.BSIM4wnigbinv * BSIM4temp.Inv_W
                                       + model.BSIM4pnigbinv * BSIM4temp.Inv_LW;
                 BSIM4temp.pParam->BSIM4ntox = model.BSIM4ntox
                                    + model.BSIM4lntox * BSIM4temp.Inv_L
                                    + model.BSIM4wntox * BSIM4temp.Inv_W
                                    + model.BSIM4pntox * BSIM4temp.Inv_LW;
                 BSIM4temp.pParam->BSIM4eigbinv = model.BSIM4eigbinv
                                       + model.BSIM4leigbinv * BSIM4temp.Inv_L
                                       + model.BSIM4weigbinv * BSIM4temp.Inv_W
                                       + model.BSIM4peigbinv * BSIM4temp.Inv_LW;
                 BSIM4temp.pParam->BSIM4pigcd = model.BSIM4pigcd
                                     + model.BSIM4lpigcd * BSIM4temp.Inv_L
                                     + model.BSIM4wpigcd * BSIM4temp.Inv_W
                                     + model.BSIM4ppigcd * BSIM4temp.Inv_LW;
                 BSIM4temp.pParam->BSIM4poxedge = model.BSIM4poxedge
                                       + model.BSIM4lpoxedge * BSIM4temp.Inv_L
                                       + model.BSIM4wpoxedge * BSIM4temp.Inv_W
                                       + model.BSIM4ppoxedge * BSIM4temp.Inv_LW;
                 BSIM4temp.pParam->BSIM4xrcrg1 = model.BSIM4xrcrg1
                                      + model.BSIM4lxrcrg1 * BSIM4temp.Inv_L
                                      + model.BSIM4wxrcrg1 * BSIM4temp.Inv_W
                                      + model.BSIM4pxrcrg1 * BSIM4temp.Inv_LW;
                 BSIM4temp.pParam->BSIM4xrcrg2 = model.BSIM4xrcrg2
                                      + model.BSIM4lxrcrg2 * BSIM4temp.Inv_L
                                      + model.BSIM4wxrcrg2 * BSIM4temp.Inv_W
                                      + model.BSIM4pxrcrg2 * BSIM4temp.Inv_LW;
                 BSIM4temp.pParam->BSIM4lambda = model.BSIM4lambda
                                      + model.BSIM4llambda * BSIM4temp.Inv_L
                                      + model.BSIM4wlambda * BSIM4temp.Inv_W
                                      + model.BSIM4plambda * BSIM4temp.Inv_LW;
                 BSIM4temp.pParam->BSIM4vtl = model.BSIM4vtl
                                      + model.BSIM4lvtl * BSIM4temp.Inv_L
                                      + model.BSIM4wvtl * BSIM4temp.Inv_W
                                      + model.BSIM4pvtl * BSIM4temp.Inv_LW;
                 BSIM4temp.pParam->BSIM4xn = model.BSIM4xn
                                      + model.BSIM4lxn * BSIM4temp.Inv_L
                                      + model.BSIM4wxn * BSIM4temp.Inv_W
                                      + model.BSIM4pxn * BSIM4temp.Inv_LW;
                 BSIM4temp.pParam->BSIM4vfbsdoff = model.BSIM4vfbsdoff
                                      + model.BSIM4lvfbsdoff * BSIM4temp.Inv_L
                                      + model.BSIM4wvfbsdoff * BSIM4temp.Inv_W
                                      + model.BSIM4pvfbsdoff * BSIM4temp.Inv_LW;
                 BSIM4temp.pParam->BSIM4tvfbsdoff = model.BSIM4tvfbsdoff
                                      + model.BSIM4ltvfbsdoff * BSIM4temp.Inv_L
                                      + model.BSIM4wtvfbsdoff * BSIM4temp.Inv_W
                                      + model.BSIM4ptvfbsdoff * BSIM4temp.Inv_LW;

         BSIM4temp.pParam->BSIM4cgsl = model.BSIM4cgsl
                    + model.BSIM4lcgsl * BSIM4temp.Inv_L
                    + model.BSIM4wcgsl * BSIM4temp.Inv_W
                    + model.BSIM4pcgsl * BSIM4temp.Inv_LW;
         BSIM4temp.pParam->BSIM4cgdl = model.BSIM4cgdl
                    + model.BSIM4lcgdl * BSIM4temp.Inv_L
                    + model.BSIM4wcgdl * BSIM4temp.Inv_W
                    + model.BSIM4pcgdl * BSIM4temp.Inv_LW;
         BSIM4temp.pParam->BSIM4ckappas = model.BSIM4ckappas
                       + model.BSIM4lckappas * BSIM4temp.Inv_L
                       + model.BSIM4wckappas * BSIM4temp.Inv_W
                       + model.BSIM4pckappas * BSIM4temp.Inv_LW;
                 BSIM4temp.pParam->BSIM4ckappad = model.BSIM4ckappad
                                       + model.BSIM4lckappad * BSIM4temp.Inv_L
                                       + model.BSIM4wckappad * BSIM4temp.Inv_W
                                       + model.BSIM4pckappad * BSIM4temp.Inv_LW;
         BSIM4temp.pParam->BSIM4cf = model.BSIM4cf
                  + model.BSIM4lcf * BSIM4temp.Inv_L
                  + model.BSIM4wcf * BSIM4temp.Inv_W
                  + model.BSIM4pcf * BSIM4temp.Inv_LW;
         BSIM4temp.pParam->BSIM4clc = model.BSIM4clc
                   + model.BSIM4lclc * BSIM4temp.Inv_L
                   + model.BSIM4wclc * BSIM4temp.Inv_W
                   + model.BSIM4pclc * BSIM4temp.Inv_LW;
         BSIM4temp.pParam->BSIM4cle = model.BSIM4cle
                   + model.BSIM4lcle * BSIM4temp.Inv_L
                   + model.BSIM4wcle * BSIM4temp.Inv_W
                   + model.BSIM4pcle * BSIM4temp.Inv_LW;
         BSIM4temp.pParam->BSIM4vfbcv = model.BSIM4vfbcv
                     + model.BSIM4lvfbcv * BSIM4temp.Inv_L
                     + model.BSIM4wvfbcv * BSIM4temp.Inv_W
                     + model.BSIM4pvfbcv * BSIM4temp.Inv_LW;
                 BSIM4temp.pParam->BSIM4acde = model.BSIM4acde
                                    + model.BSIM4lacde * BSIM4temp.Inv_L
                                    + model.BSIM4wacde * BSIM4temp.Inv_W
                                    + model.BSIM4pacde * BSIM4temp.Inv_LW;
                 BSIM4temp.pParam->BSIM4moin = model.BSIM4moin
                                    + model.BSIM4lmoin * BSIM4temp.Inv_L
                                    + model.BSIM4wmoin * BSIM4temp.Inv_W
                                    + model.BSIM4pmoin * BSIM4temp.Inv_LW;
                 BSIM4temp.pParam->BSIM4noff = model.BSIM4noff
                                    + model.BSIM4lnoff * BSIM4temp.Inv_L
                                    + model.BSIM4wnoff * BSIM4temp.Inv_W
                                    + model.BSIM4pnoff * BSIM4temp.Inv_LW;
                 BSIM4temp.pParam->BSIM4voffcv = model.BSIM4voffcv
                                      + model.BSIM4lvoffcv * BSIM4temp.Inv_L
                                      + model.BSIM4wvoffcv * BSIM4temp.Inv_W
                                      + model.BSIM4pvoffcv * BSIM4temp.Inv_LW;
                 BSIM4temp.pParam->BSIM4kvth0we = model.BSIM4kvth0we
                                      + model.BSIM4lkvth0we * BSIM4temp.Inv_L
                                      + model.BSIM4wkvth0we * BSIM4temp.Inv_W
                                      + model.BSIM4pkvth0we * BSIM4temp.Inv_LW;
                 BSIM4temp.pParam->BSIM4k2we = model.BSIM4k2we
                                      + model.BSIM4lk2we * BSIM4temp.Inv_L
                                      + model.BSIM4wk2we * BSIM4temp.Inv_W
                                      + model.BSIM4pk2we * BSIM4temp.Inv_LW;
                 BSIM4temp.pParam->BSIM4ku0we = model.BSIM4ku0we
                                      + model.BSIM4lku0we * BSIM4temp.Inv_L
                                      + model.BSIM4wku0we * BSIM4temp.Inv_W
                                      + model.BSIM4pku0we * BSIM4temp.Inv_LW;
                 BSIM4temp.pParam->BSIM4abulkCVfactor = 1.0 + std::pow((BSIM4temp.pParam->BSIM4clc
                         /BSIM4temp.pParam->BSIM4leffCV),
                        BSIM4temp.pParam->BSIM4cle);

         BSIM4temp.T0 = (BSIM4temp.TRatio - 1.0);

        BSIM4temp.PowWeffWr = std::pow(BSIM4temp.pParam->BSIM4weffCJ * 1.0e6,BSIM4temp.pParam->BSIM4wr) * BSIM4instance.BSIM4nf;

              BSIM4temp.T1 = BSIM4temp.T2 = BSIM4temp.T3 = BSIM4temp.T4 = 0.0;
         BSIM4temp.pParam->BSIM4ucs =BSIM4temp.pParam->BSIM4ucs * std::pow(BSIM4temp.TRatio,BSIM4temp.pParam->BSIM4ucste);
              if(model.BSIM4tempMod == 0)
          {
                 BSIM4temp.pParam->BSIM4ua =BSIM4temp.pParam->BSIM4ua +BSIM4temp.pParam->BSIM4ua1 *  BSIM4temp.T0;
                 BSIM4temp.pParam->BSIM4ub =BSIM4temp.pParam->BSIM4ub +BSIM4temp.pParam->BSIM4ub1 *  BSIM4temp.T0;
                 BSIM4temp.pParam->BSIM4uc =BSIM4temp.pParam->BSIM4uc +BSIM4temp.pParam->BSIM4uc1 *  BSIM4temp.T0;
                 BSIM4temp.pParam->BSIM4ud =BSIM4temp.pParam->BSIM4ud +BSIM4temp.pParam->BSIM4ud1 *  BSIM4temp.T0;
                     BSIM4temp.pParam->BSIM4vsattemp =BSIM4temp.pParam->BSIM4vsat -BSIM4temp.pParam->BSIM4at *  BSIM4temp.T0;
              BSIM4temp.T10 =BSIM4temp.pParam->BSIM4prt *  BSIM4temp.T0;
              if(model.BSIM4rdsMod)
              {
              /* External Rd(V) */
              BSIM4temp.T1 =BSIM4temp.pParam->BSIM4rdw + BSIM4temp.T10;
                      BSIM4temp.T2 = model.BSIM4rdwmin + BSIM4temp.T10;
              /* External Rs(V) */
              BSIM4temp.T3 =BSIM4temp.pParam->BSIM4rsw + BSIM4temp.T10;
                      BSIM4temp.T4 = model.BSIM4rswmin + BSIM4temp.T10;
                      }
              /* Internal Rds(V) in IV */
                 BSIM4temp.pParam->BSIM4rds0 = (BSIM4temp.pParam->BSIM4rdsw + BSIM4temp.T10)
                        * BSIM4instance.BSIM4nf / BSIM4temp.PowWeffWr;
             BSIM4temp.pParam->BSIM4rdswmin = (model.BSIM4rdswmin + BSIM4temp.T10)
                        * BSIM4instance.BSIM4nf / BSIM4temp.PowWeffWr;
                  }
          else
          {
              if (model.BSIM4tempMod == 3)
              {
                 BSIM4temp.pParam->BSIM4ua =BSIM4temp.pParam->BSIM4ua * std::pow(BSIM4temp.TRatio,BSIM4temp.pParam->BSIM4ua1) ;
                 BSIM4temp.pParam->BSIM4ub =BSIM4temp.pParam->BSIM4ub * std::pow(BSIM4temp.TRatio,BSIM4temp.pParam->BSIM4ub1);
                 BSIM4temp.pParam->BSIM4uc =BSIM4temp.pParam->BSIM4uc * std::pow(BSIM4temp.TRatio,BSIM4temp.pParam->BSIM4uc1);
                 BSIM4temp.pParam->BSIM4ud =BSIM4temp.pParam->BSIM4ud * std::pow(BSIM4temp.TRatio,BSIM4temp.pParam->BSIM4ud1);
              }
              else
              {
                  /* tempMod = 1, 2 */
                 BSIM4temp.pParam->BSIM4ua =BSIM4temp.pParam->BSIM4ua * (1.0 +BSIM4temp.pParam->BSIM4ua1 * BSIM4temp.delTemp) ;
                 BSIM4temp.pParam->BSIM4ub =BSIM4temp.pParam->BSIM4ub * (1.0 +BSIM4temp.pParam->BSIM4ub1 * BSIM4temp.delTemp);
                 BSIM4temp.pParam->BSIM4uc =BSIM4temp.pParam->BSIM4uc * (1.0 +BSIM4temp.pParam->BSIM4uc1 * BSIM4temp.delTemp);
                 BSIM4temp.pParam->BSIM4ud =BSIM4temp.pParam->BSIM4ud * (1.0 +BSIM4temp.pParam->BSIM4ud1 * BSIM4temp.delTemp);
              }
             BSIM4temp.pParam->BSIM4vsattemp =BSIM4temp.pParam->BSIM4vsat * (1.0 -BSIM4temp.pParam->BSIM4at * BSIM4temp.delTemp);
              BSIM4temp.T10 = 1.0 +BSIM4temp.pParam->BSIM4prt * BSIM4temp.delTemp;
              if(model.BSIM4rdsMod)
              {
              /* External Rd(V) */
              BSIM4temp.T1 =BSIM4temp.pParam->BSIM4rdw * BSIM4temp.T10;
                      BSIM4temp.T2 = model.BSIM4rdwmin * BSIM4temp.T10;

              /* External Rs(V) */
              BSIM4temp.T3 =BSIM4temp.pParam->BSIM4rsw * BSIM4temp.T10;
                      BSIM4temp.T4 = model.BSIM4rswmin * BSIM4temp.T10;
                      }
              /* Internal Rds(V) in IV */
                 BSIM4temp.pParam->BSIM4rds0 =BSIM4temp.pParam->BSIM4rdsw * BSIM4temp.T10 * BSIM4instance.BSIM4nf / BSIM4temp.PowWeffWr;
             BSIM4temp.pParam->BSIM4rdswmin = model.BSIM4rdswmin * BSIM4temp.T10 * BSIM4instance.BSIM4nf / BSIM4temp.PowWeffWr;
                  }

          if (BSIM4temp.T1 < 0.0)
          {   BSIM4temp.T1 = 0.0;
              printf("Warning: Rdw at current temperature is negative; set to 0.\n");
          }
          if (BSIM4temp.T2 < 0.0)
                  {   BSIM4temp.T2 = 0.0;
                      printf("Warning: Rdwmin at current temperature is negative; set to 0.\n");
                  }
         BSIM4temp.pParam->BSIM4rd0 = BSIM4temp.T1 / BSIM4temp.PowWeffWr;
                 BSIM4temp.pParam->BSIM4rdwmin = BSIM4temp.T2 / BSIM4temp.PowWeffWr;
                  if (BSIM4temp.T3 < 0.0)
                  {   BSIM4temp.T3 = 0.0;
                      printf("Warning: Rsw at current temperature is negative; set to 0.\n");
                  }
                  if (BSIM4temp.T4 < 0.0)
                  {   BSIM4temp.T4 = 0.0;
                      printf("Warning: Rswmin at current temperature is negative; set to 0.\n");
                  }
                 BSIM4temp.pParam->BSIM4rs0 = BSIM4temp.T3 / BSIM4temp.PowWeffWr;
                 BSIM4temp.pParam->BSIM4rswmin = BSIM4temp.T4 / BSIM4temp.PowWeffWr;

                  if (BSIM4temp.pParam->BSIM4u0 > 1.0)
                     BSIM4temp.pParam->BSIM4u0 =BSIM4temp.pParam->BSIM4u0 / 1.0e4;

                  /* mobility channel length dependence */
                 BSIM4temp.T5 = 1.0 -BSIM4temp.pParam->BSIM4up * std::exp( -BSIM4temp.pParam->BSIM4leff /BSIM4temp.pParam->BSIM4lp);
                 BSIM4temp.pParam->BSIM4u0temp =BSIM4temp.pParam->BSIM4u0 * BSIM4temp.T5
                      * std::pow(BSIM4temp.TRatio,BSIM4temp.pParam->BSIM4ute);
                  if (BSIM4temp.pParam->BSIM4eu < 0.0)
                  {  BSIM4temp.pParam->BSIM4eu = 0.0;
              printf("Warning: eu has been negative; reset to 0.0.\n");
          }
                  if (BSIM4temp.pParam->BSIM4ucs < 0.0)
                  {  BSIM4temp.pParam->BSIM4ucs = 0.0;
              printf("Warning: ucs has been negative; reset to 0.0.\n");
          }

            BSIM4temp.pParam->BSIM4vfbsdoff =BSIM4temp.pParam->BSIM4vfbsdoff * (1.0 +BSIM4temp.pParam->BSIM4tvfbsdoff * BSIM4temp.delTemp);
            BSIM4temp.pParam->BSIM4voff =BSIM4temp.pParam->BSIM4voff * (1.0 +BSIM4temp.pParam->BSIM4tvoff * BSIM4temp.delTemp);

        BSIM4temp.pParam->BSIM4nfactor =BSIM4temp.pParam->BSIM4nfactor +BSIM4temp.pParam->BSIM4tnfactor * BSIM4temp.delTemp / BSIM4temp.Tnom;  /* v4.7 temp dep of leakage currents */
            BSIM4temp.pParam->BSIM4voffcv =BSIM4temp.pParam->BSIM4voffcv * (1.0 +BSIM4temp.pParam->BSIM4tvoffcv * BSIM4temp.delTemp);   /*    v4.7 temp dep of leakage currents */
            BSIM4temp.pParam->BSIM4eta0 =BSIM4temp.pParam->BSIM4eta0 +BSIM4temp.pParam->BSIM4teta0 * BSIM4temp.delTemp / BSIM4temp.Tnom;   /*   v4.7 temp dep of leakage currents */

                /* Source End Velocity Limit  */
                  if((model.BSIM4vtlGiven) && (model.BSIM4vtl > 0.0) )
                  {
                     if(model.BSIM4lc < 0.0)BSIM4temp.pParam->BSIM4lc = 0.0;
                     else  BSIM4temp.pParam->BSIM4lc = model.BSIM4lc ;
                      BSIM4temp.T0 =BSIM4temp.pParam->BSIM4leff / (BSIM4temp.pParam->BSIM4xn *BSIM4temp.pParam->BSIM4leff +BSIM4temp.pParam->BSIM4lc);
                    BSIM4temp.pParam->BSIM4tfactor = (1.0 -  BSIM4temp.T0) / (1.0 +  BSIM4temp.T0 );
                  }

                 BSIM4temp.pParam->BSIM4cgdo = (model.BSIM4cgdo +BSIM4temp.pParam->BSIM4cf)
                    *BSIM4temp.pParam->BSIM4weffCV;
                 BSIM4temp.pParam->BSIM4cgso = (model.BSIM4cgso +BSIM4temp.pParam->BSIM4cf)
                    *BSIM4temp.pParam->BSIM4weffCV;
                 BSIM4temp.pParam->BSIM4cgbo = model.BSIM4cgbo *BSIM4temp.pParam->BSIM4leffCV * BSIM4instance.BSIM4nf;

                  if (!model.BSIM4ndepGiven && model.BSIM4gamma1Given)
                  {    BSIM4temp.T0 =BSIM4temp.pParam->BSIM4gamma1 * model.BSIM4coxe;
                     BSIM4temp.pParam->BSIM4ndep = 3.01248e22 *  BSIM4temp.T0 *  BSIM4temp.T0;
                  }

         BSIM4temp.pParam->BSIM4phi = BSIM4temp.Vtm0 * std::log(BSIM4temp.pParam->BSIM4ndep /  BSIM4temp.ni)
                   +BSIM4temp.pParam->BSIM4phin + 0.4;

             BSIM4temp.pParam->BSIM4sqrtPhi = std::sqrt(BSIM4temp.pParam->BSIM4phi);
             BSIM4temp.pParam->BSIM4phis3 =BSIM4temp.pParam->BSIM4sqrtPhi *BSIM4temp.pParam->BSIM4phi;

                 BSIM4temp.pParam->BSIM4Xdep0 = std::sqrt(2.0 * BSIM4temp.epssub / (Charge_q
                     *BSIM4temp.pParam->BSIM4ndep * 1.0e6))
                                     *BSIM4temp.pParam->BSIM4sqrtPhi;
                 BSIM4temp.pParam->BSIM4sqrtXdep0 = std::sqrt(BSIM4temp.pParam->BSIM4Xdep0);

          if(model.BSIM4mtrlMod == 0)
           BSIM4temp.pParam->BSIM4litl = std::sqrt(3.0 * 3.9 /   BSIM4temp.epsrox *BSIM4temp.pParam->BSIM4xj *  BSIM4temp.toxe);
          else
           BSIM4temp.pParam->BSIM4litl = std::sqrt(model.BSIM4epsrsub/ BSIM4temp.epsrox *BSIM4temp.pParam->BSIM4xj *  BSIM4temp.toxe);

                 BSIM4temp.pParam->BSIM4vbi = BSIM4temp.Vtm0 * std::log(BSIM4temp.pParam->BSIM4nsd
                       *BSIM4temp.pParam->BSIM4ndep / ( BSIM4temp.ni *  BSIM4temp.ni));

          if (model.BSIM4mtrlMod == 0)
          {
            if (BSIM4temp.pParam->BSIM4ngate > 0.0)
            {  BSIM4temp.pParam->BSIM4vfbsd = BSIM4temp.Vtm0 * std::log(BSIM4temp.pParam->BSIM4ngate
                                         /BSIM4temp.pParam->BSIM4nsd);
             }
            else
             BSIM4temp.pParam->BSIM4vfbsd = 0.0;
          }
          else
          {
             BSIM4temp.T0 = BSIM4temp.Vtm0 * std::log(BSIM4temp.pParam->BSIM4nsd/ BSIM4temp.ni);
            BSIM4temp.T1 = 0.5 *  BSIM4temp.Eg0;
            if( BSIM4temp.T0 > BSIM4temp.T1)
               BSIM4temp.T0 = BSIM4temp.T1;
            BSIM4temp.T2 = model.BSIM4easub + BSIM4temp.T1 - model.BSIM4type *  BSIM4temp.T0;
           BSIM4temp.pParam->BSIM4vfbsd = model.BSIM4phig - BSIM4temp.T2;
          }

                 BSIM4temp.pParam->BSIM4cdep0 = std::sqrt(Charge_q * BSIM4temp.epssub
                     *BSIM4temp.pParam->BSIM4ndep * 1.0e6 / 2.0
                     /BSIM4temp.pParam->BSIM4phi);

                 BSIM4temp.pParam->BSIM4ToxRatio = std::exp(BSIM4temp.pParam->BSIM4ntox
                    * std::log(model.BSIM4toxref /  BSIM4temp.toxe))
                    /  BSIM4temp.toxe /  BSIM4temp.toxe;
                 BSIM4temp.pParam->BSIM4ToxRatioEdge = std::exp(BSIM4temp.pParam->BSIM4ntox
                                            * std::log(model.BSIM4toxref
                                            / ( BSIM4temp.toxe *BSIM4temp.pParam->BSIM4poxedge)))
                                            /  BSIM4temp.toxe /  BSIM4temp.toxe
                                            /BSIM4temp.pParam->BSIM4poxedge /BSIM4temp.pParam->BSIM4poxedge;

                 BSIM4temp.pParam->BSIM4Aechvb = (model.BSIM4type == BSIM4_NMOS) ? 4.97232e-7 : 3.42537e-7;
                 BSIM4temp.pParam->BSIM4Bechvb = (model.BSIM4type == BSIM4_NMOS) ? 7.45669e11 : 1.16645e12;

                const std::string& version = model.BSIM4version;
                if (version != "4.8.1" && 
                     version.substr(0, 4) != "4.81" &&
                     version != "4.8.2" && 
                     version.substr(0, 4) != "4.82")
                {  /* check only for version <= 4.80 */
                  
                 BSIM4temp.pParam->BSIM4AechvbEdgeS =BSIM4temp.pParam->BSIM4Aechvb *BSIM4temp.pParam->BSIM4weff
                      * model.BSIM4dlcig *BSIM4temp.pParam->BSIM4ToxRatioEdge;
                 BSIM4temp.pParam->BSIM4AechvbEdgeD =BSIM4temp.pParam->BSIM4Aechvb *BSIM4temp.pParam->BSIM4weff
                      * model.BSIM4dlcigd *BSIM4temp.pParam->BSIM4ToxRatioEdge;

                  }
                  else
                  {
                      if (model.BSIM4dlcig < 0.0)
                  {
                          printf("Warning: dlcig = %g is negative. Set to zero.\n", model.BSIM4dlcig);
                          model.BSIM4dlcig = 0.0;
                  }
                 BSIM4temp.pParam->BSIM4AechvbEdgeS =BSIM4temp.pParam->BSIM4Aechvb *BSIM4temp.pParam->BSIM4weff
                      * model.BSIM4dlcig *BSIM4temp.pParam->BSIM4ToxRatioEdge;
                  if (model.BSIM4dlcigd < 0.0)
                  {
                          printf("Warning: dlcigd = %g is negative. Set to zero.\n", model.BSIM4dlcigd);
                          model.BSIM4dlcigd = 0.0;
                  }
                 BSIM4temp.pParam->BSIM4AechvbEdgeD =BSIM4temp.pParam->BSIM4Aechvb *BSIM4temp.pParam->BSIM4weff
                      * model.BSIM4dlcigd *BSIM4temp.pParam->BSIM4ToxRatioEdge;

                  }

                 BSIM4temp.pParam->BSIM4BechvbEdge =  -BSIM4temp.pParam->BSIM4Bechvb
                      *  BSIM4temp.toxe *BSIM4temp.pParam->BSIM4poxedge;
                 BSIM4temp.pParam->BSIM4Aechvb *=BSIM4temp.pParam->BSIM4weff *BSIM4temp.pParam->BSIM4leff
                       *BSIM4temp.pParam->BSIM4ToxRatio;
                 BSIM4temp.pParam->BSIM4Bechvb *= - BSIM4temp.toxe;




                 BSIM4temp.pParam->BSIM4mstar = 0.5 + atan(BSIM4temp.pParam->BSIM4minv) / PI;
                 BSIM4temp.pParam->BSIM4mstarcv = 0.5 + atan(BSIM4temp.pParam->BSIM4minvcv) / PI;
                 BSIM4temp.pParam->BSIM4voffcbn = BSIM4temp.pParam->BSIM4voff + model.BSIM4voffl /BSIM4temp.pParam->BSIM4leff;
                 BSIM4temp.pParam->BSIM4voffcbncv = BSIM4temp.pParam->BSIM4voffcv + model.BSIM4voffcvl /BSIM4temp.pParam->BSIM4leff;

                 BSIM4temp.pParam->BSIM4ldeb = std::sqrt(BSIM4temp.epssub * BSIM4temp.Vtm0 / (Charge_q
                                    *BSIM4temp.pParam->BSIM4ndep * 1.0e6)) / 3.0;
                 BSIM4temp.pParam->BSIM4acde *= std::pow((BSIM4temp.pParam->BSIM4ndep / 2.0e16), -0.25);


                  if (model.BSIM4k1Given || model.BSIM4k2Given)
              {   if (!model.BSIM4k1Given)
                  {   fprintf(stdout, "Warning: k1 should be specified with k2.\n");
                         BSIM4temp.pParam->BSIM4k1 = 0.53;
                      }
                      if (!model.BSIM4k2Given)
                  {   fprintf(stdout, "Warning: k2 should be specified with k1.\n");
                         BSIM4temp.pParam->BSIM4k2 = -0.0186;
                      }
                      if (model.BSIM4nsubGiven)
                          fprintf(stdout, "Warning: nsub is ignored because k1 or k2 is given.\n");
                      if (model.BSIM4xtGiven)
                          fprintf(stdout, "Warning: xt is ignored because k1 or k2 is given.\n");
                      if (model.BSIM4vbxGiven)
                          fprintf(stdout, "Warning: vbx is ignored because k1 or k2 is given.\n");
                      if (model.BSIM4gamma1Given)
                          fprintf(stdout, "Warning: gamma1 is ignored because k1 or k2 is given.\n");
                      if (model.BSIM4gamma2Given)
                          fprintf(stdout, "Warning: gamma2 is ignored because k1 or k2 is given.\n");
                  }
                  else
              {   if (!model.BSIM4vbxGiven)
                         BSIM4temp.pParam->BSIM4vbx =BSIM4temp.pParam->BSIM4phi - 7.7348e-4
                                           *BSIM4temp.pParam->BSIM4ndep
                       *BSIM4temp.pParam->BSIM4xt *BSIM4temp.pParam->BSIM4xt;
                  if (BSIM4temp.pParam->BSIM4vbx > 0.0)
                 BSIM4temp.pParam->BSIM4vbx =  -BSIM4temp.pParam->BSIM4vbx;
                  if (BSIM4temp.pParam->BSIM4vbm > 0.0)
                         BSIM4temp.pParam->BSIM4vbm =  -BSIM4temp.pParam->BSIM4vbm;

                      if (!model.BSIM4gamma1Given)
                         BSIM4temp.pParam->BSIM4gamma1 = 5.753e-12
                          * std::sqrt(BSIM4temp.pParam->BSIM4ndep)
                                              / model.BSIM4coxe;
                      if (!model.BSIM4gamma2Given)
                         BSIM4temp.pParam->BSIM4gamma2 = 5.753e-12
                          * std::sqrt(BSIM4temp.pParam->BSIM4nsub)
                                              / model.BSIM4coxe;

                       BSIM4temp.T0 =BSIM4temp.pParam->BSIM4gamma1 -BSIM4temp.pParam->BSIM4gamma2;
                      BSIM4temp.T1 = std::sqrt(BSIM4temp.pParam->BSIM4phi -BSIM4temp.pParam->BSIM4vbx)
             -BSIM4temp.pParam->BSIM4sqrtPhi;
                      BSIM4temp.T2 = std::sqrt(BSIM4temp.pParam->BSIM4phi * (BSIM4temp.pParam->BSIM4phi
             -BSIM4temp.pParam->BSIM4vbm)) -BSIM4temp.pParam->BSIM4phi;
                     BSIM4temp.pParam->BSIM4k2 =  BSIM4temp.T0 * BSIM4temp.T1 / (2.0 * BSIM4temp.T2 +BSIM4temp.pParam->BSIM4vbm);
                     BSIM4temp.pParam->BSIM4k1 =BSIM4temp.pParam->BSIM4gamma2 - 2.0
                      *BSIM4temp.pParam->BSIM4k2 * std::sqrt(BSIM4temp.pParam->BSIM4phi
                      -BSIM4temp.pParam->BSIM4vbm);
                  }

                  if (!model.BSIM4vfbGiven)
                  {
            if (model.BSIM4vth0Given)
                      {  BSIM4temp.pParam->BSIM4vfb = model.BSIM4type *BSIM4temp.pParam->BSIM4vth0
                                           -BSIM4temp.pParam->BSIM4phi -BSIM4temp.pParam->BSIM4k1
                                           *BSIM4temp.pParam->BSIM4sqrtPhi;
                      }
                      else
              {
            if ((model.BSIM4mtrlMod) && (model.BSIM4phigGiven) &&
                (model.BSIM4nsubGiven))
              {
                 BSIM4temp.T0 = BSIM4temp.Vtm0 * std::log(BSIM4temp.pParam->BSIM4nsub/ BSIM4temp.ni);
                BSIM4temp.T1 = 0.5 *  BSIM4temp.Eg0;
                if( BSIM4temp.T0 > BSIM4temp.T1)
                   BSIM4temp.T0 = BSIM4temp.T1;
                BSIM4temp.T2 = model.BSIM4easub + BSIM4temp.T1 + model.BSIM4type *  BSIM4temp.T0;
               BSIM4temp.pParam->BSIM4vfb = model.BSIM4phig - BSIM4temp.T2;
              }
            else
              {
               BSIM4temp.pParam->BSIM4vfb = -1.0;
              }
              }
                  }
                   if (!model.BSIM4vth0Given)
                  {  BSIM4temp.pParam->BSIM4vth0 = model.BSIM4type * (BSIM4temp.pParam->BSIM4vfb
                                        +BSIM4temp.pParam->BSIM4phi +BSIM4temp.pParam->BSIM4k1
                                        *BSIM4temp.pParam->BSIM4sqrtPhi);
                  }

                 BSIM4temp.pParam->BSIM4k1ox =BSIM4temp.pParam->BSIM4k1 *  BSIM4temp.toxe
                                    / model.BSIM4toxm;

                   BSIM4temp.tmp = std::sqrt(BSIM4temp.epssub / ( BSIM4temp.epsrox * EPS0) *  BSIM4temp.toxe *BSIM4temp.pParam->BSIM4Xdep0);
               BSIM4temp.T0 =BSIM4temp.pParam->BSIM4dsub *BSIM4temp.pParam->BSIM4leff /  BSIM4temp.tmp;
                  if ( BSIM4temp.T0 < EXP_THRESHOLD)
              {   BSIM4temp.T1 = std::exp( BSIM4temp.T0);
                      BSIM4temp.T2 = BSIM4temp.T1 - 1.0;
                      BSIM4temp.T3 = BSIM4temp.T2 * BSIM4temp.T2;
                      BSIM4temp.T4 = BSIM4temp.T3 + 2.0 * BSIM4temp.T1 * MIN_EXP;
                     BSIM4temp.pParam->BSIM4theta0vb0 = BSIM4temp.T1 / BSIM4temp.T4;
                  }
                  else
                     BSIM4temp.pParam->BSIM4theta0vb0 = 1.0 / (MAX_EXP - 2.0);

               BSIM4temp.T0 =BSIM4temp.pParam->BSIM4drout *BSIM4temp.pParam->BSIM4leff /  BSIM4temp.tmp;
              if ( BSIM4temp.T0 < EXP_THRESHOLD)
                  {   BSIM4temp.T1 = std::exp( BSIM4temp.T0);
                      BSIM4temp.T2 = BSIM4temp.T1 - 1.0;
                      BSIM4temp.T3 = BSIM4temp.T2 * BSIM4temp.T2;
                      BSIM4temp.T4 = BSIM4temp.T3 + 2.0 * BSIM4temp.T1 * MIN_EXP;
                      BSIM4temp.T5 = BSIM4temp.T1 / BSIM4temp.T4;
                  }
                  else
                      BSIM4temp.T5 = 1.0 / (MAX_EXP - 2.0); /* 3.0 * MIN_EXP omitted */
                 BSIM4temp.pParam->BSIM4thetaRout =BSIM4temp.pParam->BSIM4pdibl1 * BSIM4temp.T5
                                         +BSIM4temp.pParam->BSIM4pdibl2;

                   BSIM4temp.tmp = std::sqrt(BSIM4temp.pParam->BSIM4Xdep0);
                   BSIM4temp.tmp1 =BSIM4temp.pParam->BSIM4vbi -BSIM4temp.pParam->BSIM4phi;
                   BSIM4temp.tmp2 = model.BSIM4factor1 *  BSIM4temp.tmp;

                   BSIM4temp.T0 =BSIM4temp.pParam->BSIM4dvt1w *BSIM4temp.pParam->BSIM4weff
                     *BSIM4temp.pParam->BSIM4leff /  BSIM4temp.tmp2;
                  if ( BSIM4temp.T0 < EXP_THRESHOLD)
                  {   BSIM4temp.T1 = std::exp( BSIM4temp.T0);
                      BSIM4temp.T2 = BSIM4temp.T1 - 1.0;
                      BSIM4temp.T3 = BSIM4temp.T2 * BSIM4temp.T2;
                      BSIM4temp.T4 = BSIM4temp.T3 + 2.0 * BSIM4temp.T1 * MIN_EXP;
                       BSIM4temp.T8 = BSIM4temp.T1 / BSIM4temp.T4;
                  }
                  else
                       BSIM4temp.T8 = 1.0 / (MAX_EXP - 2.0);
                   BSIM4temp.T0 =BSIM4temp.pParam->BSIM4dvt0w *  BSIM4temp.T8;
                   BSIM4temp.T8 =  BSIM4temp.T0 *  BSIM4temp.tmp1;

                   BSIM4temp.T0 =BSIM4temp.pParam->BSIM4dvt1 *BSIM4temp.pParam->BSIM4leff /  BSIM4temp.tmp2;
                  if ( BSIM4temp.T0 < EXP_THRESHOLD)
                  {   BSIM4temp.T1 = std::exp( BSIM4temp.T0);
                      BSIM4temp.T2 = BSIM4temp.T1 - 1.0;
                      BSIM4temp.T3 = BSIM4temp.T2 * BSIM4temp.T2;
                      BSIM4temp.T4 = BSIM4temp.T3 + 2.0 * BSIM4temp.T1 * MIN_EXP;
                       BSIM4temp.T9 = BSIM4temp.T1 / BSIM4temp.T4;
                  }
                  else
                       BSIM4temp.T9 = 1.0 / (MAX_EXP - 2.0);
                   BSIM4temp.T9 =BSIM4temp.pParam->BSIM4dvt0 *  BSIM4temp.T9 *  BSIM4temp.tmp1;

                  BSIM4temp.T4 =  BSIM4temp.toxe *BSIM4temp.pParam->BSIM4phi
                     / (BSIM4temp.pParam->BSIM4weff +BSIM4temp.pParam->BSIM4w0);

                   BSIM4temp.T0 = std::sqrt(1.0 +BSIM4temp.pParam->BSIM4lpe0 /BSIM4temp.pParam->BSIM4leff);
                  if((model.BSIM4tempMod == 1) || (model.BSIM4tempMod == 0))
                    BSIM4temp.T3 = (BSIM4temp.pParam->BSIM4kt1 +BSIM4temp.pParam->BSIM4kt1l /BSIM4temp.pParam->BSIM4leff)
                            * (BSIM4temp.TRatio - 1.0);
                  if((model.BSIM4tempMod == 2)||(model.BSIM4tempMod == 3))
                        BSIM4temp.T3 = -BSIM4temp.pParam->BSIM4kt1 * (BSIM4temp.TRatio - 1.0);

                  BSIM4temp.T5 =BSIM4temp.pParam->BSIM4k1ox * ( BSIM4temp.T0 - 1.0) *BSIM4temp.pParam->BSIM4sqrtPhi
                     + BSIM4temp.T3;
                 BSIM4temp.pParam->BSIM4vfbzbfactor = -  BSIM4temp.T8 -  BSIM4temp.T9 +BSIM4temp.pParam->BSIM4k3 * BSIM4temp.T4 + BSIM4temp.T5
                       -BSIM4temp.pParam->BSIM4phi -BSIM4temp.pParam->BSIM4k1 *BSIM4temp.pParam->BSIM4sqrtPhi;

          /* stress effect */

              BSIM4temp.wlod = model.BSIM4wlod;
              if (model.BSIM4wlod < 0.0)
              {   fprintf(stderr, "Warning: WLOD = %g is less than 0. 0.0 is used\n",model.BSIM4wlod);
                      BSIM4temp.wlod = 0.0;
              }
                   BSIM4temp.T0 = std::pow(BSIM4temp.Lnew, model.BSIM4llodku0);
          BSIM4temp.W_tmp = BSIM4temp.Wnew + BSIM4temp.wlod;
                  BSIM4temp.T1 = std::pow(BSIM4temp.W_tmp, model.BSIM4wlodku0);
                   BSIM4temp.tmp1 = model.BSIM4lku0 /  BSIM4temp.T0 + model.BSIM4wku0 / BSIM4temp.T1
                         + model.BSIM4pku0 / ( BSIM4temp.T0 * BSIM4temp.T1);
                 BSIM4temp.pParam->BSIM4ku0 = 1.0 +  BSIM4temp.tmp1;

                   BSIM4temp.T0 = std::pow(BSIM4temp.Lnew, model.BSIM4llodvth);
                  BSIM4temp.T1 = std::pow(BSIM4temp.W_tmp, model.BSIM4wlodvth);
                   BSIM4temp.tmp1 = model.BSIM4lkvth0 /  BSIM4temp.T0 + model.BSIM4wkvth0 / BSIM4temp.T1
                       + model.BSIM4pkvth0 / ( BSIM4temp.T0 * BSIM4temp.T1);
                 BSIM4temp.pParam->BSIM4kvth0 = 1.0 +  BSIM4temp.tmp1;
         BSIM4temp.pParam->BSIM4kvth0 = std::sqrt(BSIM4temp.pParam->BSIM4kvth0 * BSIM4temp.pParam->BSIM4kvth0 + DELTA);

                   BSIM4temp.T0 = (BSIM4temp.TRatio - 1.0);
                 BSIM4temp.pParam->BSIM4ku0temp =BSIM4temp.pParam->BSIM4ku0 * (1.0 + model.BSIM4tku0 * BSIM4temp.T0) + DELTA;

                   BSIM4temp.Inv_saref = 1.0/(model.BSIM4saref + 0.5* BSIM4temp.Ldrn);
                   BSIM4temp.Inv_sbref = 1.0/(model.BSIM4sbref + 0.5* BSIM4temp.Ldrn);
         BSIM4temp.pParam->BSIM4inv_od_ref =  BSIM4temp.Inv_saref +  BSIM4temp.Inv_sbref;
         BSIM4temp.pParam->BSIM4rho_ref = model.BSIM4ku0 /BSIM4temp.pParam->BSIM4ku0temp *BSIM4temp.pParam->BSIM4inv_od_ref;

                  /*high k*/
                  /*Calculate VgsteffVth for mobMod=3*/
                  if(model.BSIM4mobMod==3)
                  { /*Calculate  BSIM4temp.n @ Vbs=Vds=0*/
                       BSIM4temp.lt1 = model.BSIM4factor1*BSIM4temp.pParam->BSIM4sqrtXdep0;
                       BSIM4temp.T0 =BSIM4temp.pParam->BSIM4dvt1 *BSIM4temp.pParam->BSIM4leff /  BSIM4temp.lt1;
                      if ( BSIM4temp.T0 < EXP_THRESHOLD)
                      {
                          BSIM4temp.T1 = std::exp( BSIM4temp.T0);
                          BSIM4temp.T2 = BSIM4temp.T1 - 1.0;
                          BSIM4temp.T3 = BSIM4temp.T2 * BSIM4temp.T2;
                          BSIM4temp.T4 = BSIM4temp.T3 + 2.0 * BSIM4temp.T1 * MIN_EXP;
                           BSIM4temp.Theta0 = BSIM4temp.T1 / BSIM4temp.T4;
                      }
                      else
                           BSIM4temp.Theta0 = 1.0 / (MAX_EXP - 2.0);

                       BSIM4temp.tmp1 = BSIM4temp.epssub /BSIM4temp.pParam->BSIM4Xdep0;
                       BSIM4temp.tmp2 =BSIM4temp.pParam->BSIM4nfactor *  BSIM4temp.tmp1;
                       BSIM4temp.tmp3 = ( BSIM4temp.tmp2 +BSIM4temp.pParam->BSIM4cdsc *  BSIM4temp.Theta0 +BSIM4temp.pParam->BSIM4cit) / model.BSIM4coxe;
                      if ( BSIM4temp.tmp3 >= -0.5)
                           BSIM4temp.n0 = 1.0 +  BSIM4temp.tmp3;
                      else
                      {
                           BSIM4temp.T0 = 1.0 / (3.0 + 8.0 *  BSIM4temp.tmp3);
                           BSIM4temp.n0 = (1.0 + 3.0 *  BSIM4temp.tmp3) *  BSIM4temp.T0;
                      }

                       BSIM4temp.T0 =  BSIM4temp.n0 * model.BSIM4vtm;
                      BSIM4temp.T1 =BSIM4temp.pParam->BSIM4voffcbn;
                      BSIM4temp.T2 = BSIM4temp.T1/ BSIM4temp.T0;
                      if (BSIM4temp.T2 < -EXP_THRESHOLD)
                      {   BSIM4temp.T3 = model.BSIM4coxe * MIN_EXP /BSIM4temp.pParam->BSIM4cdep0;
                          BSIM4temp.T4 =BSIM4temp.pParam->BSIM4mstar + BSIM4temp.T3 *  BSIM4temp.n0;
                      }
                      else if (BSIM4temp.T2 > EXP_THRESHOLD)
                      {   BSIM4temp.T3 = model.BSIM4coxe * MAX_EXP /BSIM4temp.pParam->BSIM4cdep0;
                          BSIM4temp.T4 =BSIM4temp.pParam->BSIM4mstar + BSIM4temp.T3 *  BSIM4temp.n0;
                      }
                      else
                      {  BSIM4temp.T3 = std::exp(BSIM4temp.T2)* model.BSIM4coxe /BSIM4temp.pParam->BSIM4cdep0;
                         BSIM4temp.T4 =BSIM4temp.pParam->BSIM4mstar + BSIM4temp.T3 *  BSIM4temp.n0;
                      }
                     BSIM4temp.pParam->BSIM4VgsteffVth =  BSIM4temp.T0 * std::log(2.0)/BSIM4temp.T4;
                  }

                  /* New DITS term added in 4.7 */
                   BSIM4temp.T0 =  -BSIM4temp.pParam->BSIM4dvtp3 * std::log(BSIM4temp.pParam->BSIM4leff);
                   BSIM4temp.T1 = DEXP( BSIM4temp.T0);
                 BSIM4temp.pParam->BSIM4dvtp2factor =BSIM4temp.pParam->BSIM4dvtp5 +BSIM4temp.pParam->BSIM4dvtp2 * BSIM4temp.T1;

              } /* End of SizeNotFound */

              /*  stress effect */
              if( (BSIM4instance.BSIM4sa > 0.0) && (BSIM4instance.BSIM4sb > 0.0) &&
                 ((BSIM4instance.BSIM4nf == 1.0) || ((BSIM4instance.BSIM4nf > 1.0) && (BSIM4instance.BSIM4sd > 0.0))) )
          {    BSIM4temp.Inv_sa = 0;
                   BSIM4temp.Inv_sb = 0;

                   BSIM4temp.kvsat = model.BSIM4kvsat;
          if (model.BSIM4kvsat < -1.0 )
              {   fprintf(stderr, "Warning: KVSAT = %g is too small; -1.0 is used.\n",model.BSIM4kvsat);
                   BSIM4temp.kvsat = -1.0;
                  }
                  if (model.BSIM4kvsat > 1.0)
                  {   fprintf(stderr, "Warning: KVSAT = %g is too big; 1.0 is used.\n",model.BSIM4kvsat);
                   BSIM4temp.kvsat = 1.0;
                  }

              for(i = 0; i < BSIM4instance.BSIM4nf; i++){
                     BSIM4temp.T0 = 1.0 / BSIM4instance.BSIM4nf / (BSIM4instance.BSIM4sa + 0.5* BSIM4temp.Ldrn + i * (BSIM4instance.BSIM4sd + BSIM4temp.Ldrn));
                        BSIM4temp.T1 = 1.0 / BSIM4instance.BSIM4nf / (BSIM4instance.BSIM4sb + 0.5* BSIM4temp.Ldrn + i * (BSIM4instance.BSIM4sd + BSIM4temp.Ldrn));
                     BSIM4temp.Inv_sa +=  BSIM4temp.T0;
                         BSIM4temp.Inv_sb += BSIM4temp.T1;
                  }
                   BSIM4temp.Inv_ODeff =  BSIM4temp.Inv_sa +  BSIM4temp.Inv_sb;
                   BSIM4temp.rho = model.BSIM4ku0 /BSIM4temp.pParam->BSIM4ku0temp *  BSIM4temp.Inv_ODeff;
                   BSIM4temp.T0 = (1.0 +  BSIM4temp.rho)/(1.0 +BSIM4temp.pParam->BSIM4rho_ref);
                  BSIM4instance.BSIM4u0temp =BSIM4temp.pParam->BSIM4u0temp *  BSIM4temp.T0;

                  BSIM4temp.T1 = (1.0 +  BSIM4temp.kvsat *  BSIM4temp.rho)/(1.0 +  BSIM4temp.kvsat *BSIM4temp.pParam->BSIM4rho_ref);
                  BSIM4instance.BSIM4vsattemp =BSIM4temp.pParam->BSIM4vsattemp * BSIM4temp.T1;

           BSIM4temp.OD_offset =  BSIM4temp.Inv_ODeff -BSIM4temp.pParam->BSIM4inv_od_ref;
           BSIM4temp.dvth0_lod = model.BSIM4kvth0 /BSIM4temp.pParam->BSIM4kvth0 *  BSIM4temp.OD_offset;
                   BSIM4temp.dk2_lod = model.BSIM4stk2 / std::pow(BSIM4temp.pParam->BSIM4kvth0, model.BSIM4lodk2) *
                                    BSIM4temp.OD_offset;
                   BSIM4temp.deta0_lod = model.BSIM4steta0 / std::pow(BSIM4temp.pParam->BSIM4kvth0, model.BSIM4lodeta0) *
                                      BSIM4temp.OD_offset;
          BSIM4instance.BSIM4vth0 =BSIM4temp.pParam->BSIM4vth0 +  BSIM4temp.dvth0_lod;

                  BSIM4instance.BSIM4eta0 =BSIM4temp.pParam->BSIM4eta0 +  BSIM4temp.deta0_lod;
          BSIM4instance.BSIM4k2 =BSIM4temp.pParam->BSIM4k2 +  BSIM4temp.dk2_lod;
           } else {
              BSIM4instance.BSIM4u0temp =BSIM4temp.pParam->BSIM4u0temp;
                      BSIM4instance.BSIM4vth0 =BSIM4temp.pParam->BSIM4vth0;
                      BSIM4instance.BSIM4vsattemp =BSIM4temp.pParam->BSIM4vsattemp;
                      BSIM4instance.BSIM4eta0 =BSIM4temp.pParam->BSIM4eta0;
                      BSIM4instance.BSIM4k2 =BSIM4temp.pParam->BSIM4k2;
              }

          /*  Well Proximity Effect  */
              if (model.BSIM4wpemod)
              { if( (!BSIM4instance.BSIM4scaGiven) && (!BSIM4instance.BSIM4scbGiven) && (!BSIM4instance.BSIM4sccGiven) )
        {   if((BSIM4instance.BSIM4scGiven) && (BSIM4instance.BSIM4sc > 0.0) )
                {   BSIM4temp.T1 = BSIM4instance.BSIM4sc +  BSIM4temp.Wdrn;
                    BSIM4temp.T2 = 1.0 / model.BSIM4scref;
            BSIM4instance.BSIM4sca = model.BSIM4scref * model.BSIM4scref
                    / (BSIM4instance.BSIM4sc * BSIM4temp.T1);
            BSIM4instance.BSIM4scb = ( (0.1 * BSIM4instance.BSIM4sc + 0.01 * model.BSIM4scref)
                    * std::exp(-10.0 * BSIM4instance.BSIM4sc * BSIM4temp.T2)
                    - (0.1 * BSIM4temp.T1 + 0.01 * model.BSIM4scref)
                    * std::exp(-10.0 * BSIM4temp.T1 * BSIM4temp.T2) ) /  BSIM4temp.Wdrn;
                        BSIM4instance.BSIM4scc = ( (0.05 * BSIM4instance.BSIM4sc + 0.0025 * model.BSIM4scref)
                                        * std::exp(-20.0 * BSIM4instance.BSIM4sc * BSIM4temp.T2)
                                        - (0.05 * BSIM4temp.T1 + 0.0025 * model.BSIM4scref)
                                        * std::exp(-20.0 * BSIM4temp.T1 * BSIM4temp.T2) ) /  BSIM4temp.Wdrn;
            } else {
                        fprintf(stderr, "Warning: No WPE as none of SCA, SCB, SCC, SC is given and/or SC not positive.\n");
            }
        }

               if (BSIM4instance.BSIM4sca < 0.0)
                {
                    printf("Warning: SCA = %g is negative. Set to 0.0.\n", BSIM4instance.BSIM4sca);
                    BSIM4instance.BSIM4sca = 0.0;
                }
                if (BSIM4instance.BSIM4scb < 0.0)
                {
                    printf("Warning: SCB = %g is negative. Set to 0.0.\n", BSIM4instance.BSIM4scb);
                    BSIM4instance.BSIM4scb = 0.0;
                }
                if (BSIM4instance.BSIM4scc < 0.0)
                {
                    printf("Warning: SCC = %g is negative. Set to 0.0.\n", BSIM4instance.BSIM4scc);
                    BSIM4instance.BSIM4scc = 0.0;
                }
                if (BSIM4instance.BSIM4sc < 0.0)
                {
                    printf("Warning: SC = %g is negative. Set to 0.0.\n", BSIM4instance.BSIM4sc);
                    BSIM4instance.BSIM4sc = 0.0;
                }
                /*4.6.2*/
         BSIM4temp.sceff = BSIM4instance.BSIM4sca + model.BSIM4web * BSIM4instance.BSIM4scb
                      + model.BSIM4wec * BSIM4instance.BSIM4scc;
                BSIM4instance.BSIM4vth0 +=BSIM4temp.pParam->BSIM4kvth0we *  BSIM4temp.sceff;
                BSIM4instance.BSIM4k2 += BSIM4temp.pParam->BSIM4k2we *  BSIM4temp.sceff;
        BSIM4temp.T3 =  1.0 +BSIM4temp.pParam->BSIM4ku0we *  BSIM4temp.sceff;
        if (BSIM4temp.T3 <= 0.0)
        {   BSIM4temp.T3 = 0.0;
                fprintf(stderr, "Warning: ku0we = %g is negatively too high. Negative mobility! \n",BSIM4temp.pParam->BSIM4ku0we);
        }
                BSIM4instance.BSIM4u0temp *= BSIM4temp.T3;
              }

        /* adding delvto  */
            BSIM4instance.BSIM4vth0 += BSIM4instance.BSIM4delvto;
            BSIM4instance.BSIM4vfb =BSIM4temp.pParam->BSIM4vfb + model.BSIM4type * BSIM4instance.BSIM4delvto;

        /* Instance variables calculation  */
            BSIM4temp.T3 = model.BSIM4type * BSIM4instance.BSIM4vth0
               - BSIM4instance.BSIM4vfb -BSIM4temp.pParam->BSIM4phi;
            BSIM4temp.T4 = BSIM4temp.T3 + BSIM4temp.T3;
            BSIM4temp.T5 = 2.5 * BSIM4temp.T3;
            BSIM4instance.BSIM4vtfbphi1 = (model.BSIM4type == BSIM4_NMOS) ? BSIM4temp.T4 : BSIM4temp.T5;
            if (BSIM4instance.BSIM4vtfbphi1 < 0.0)
                BSIM4instance.BSIM4vtfbphi1 = 0.0;

            BSIM4instance.BSIM4vtfbphi2 = 4.0 * BSIM4temp.T3;
            if (BSIM4instance.BSIM4vtfbphi2 < 0.0)
                BSIM4instance.BSIM4vtfbphi2 = 0.0;

            if (BSIM4instance.BSIM4k2 < 0.0)
            {    BSIM4temp.T0 = 0.5 *BSIM4temp.pParam->BSIM4k1 / BSIM4instance.BSIM4k2;
                BSIM4instance.BSIM4vbsc = 0.9 * (BSIM4temp.pParam->BSIM4phi -  BSIM4temp.T0 *  BSIM4temp.T0);
                if (BSIM4instance.BSIM4vbsc > -3.0)
                    BSIM4instance.BSIM4vbsc = -3.0;
                else if (BSIM4instance.BSIM4vbsc < -30.0)
                    BSIM4instance.BSIM4vbsc = -30.0;
            }
            else
                BSIM4instance.BSIM4vbsc = -30.0;
            if (BSIM4instance.BSIM4vbsc >BSIM4temp.pParam->BSIM4vbm)
                BSIM4instance.BSIM4vbsc =BSIM4temp.pParam->BSIM4vbm;
            BSIM4instance.BSIM4k2ox = BSIM4instance.BSIM4k2 *  BSIM4temp.toxe
                              / model.BSIM4toxm;

            BSIM4instance.BSIM4vfbzb =BSIM4temp.pParam->BSIM4vfbzbfactor
                +  model.BSIM4type * BSIM4instance.BSIM4vth0 ;

              BSIM4instance.BSIM4cgso =BSIM4temp.pParam->BSIM4cgso;
              BSIM4instance.BSIM4cgdo =BSIM4temp.pParam->BSIM4cgdo;

           BSIM4temp.lnl = std::log(BSIM4temp.pParam->BSIM4leff * 1.0e6);
           BSIM4temp.lnw = std::log(BSIM4temp.pParam->BSIM4weff * 1.0e6);
           BSIM4temp.lnnf = std::log(BSIM4instance.BSIM4nf);

           BSIM4temp.bodymode = 5;
          if( ( !model.BSIM4rbps0Given) ||
          ( !model.BSIM4rbpd0Given) )
         BSIM4temp.bodymode = 1;
          else
        if( (!model.BSIM4rbsbx0Given && !model.BSIM4rbsby0Given) ||
              (!model.BSIM4rbdbx0Given && !model.BSIM4rbdby0Given) )
           BSIM4temp.bodymode = 3;

          if(BSIM4instance.BSIM4rbodyMod == 2)
        {
          if ( BSIM4temp.bodymode == 5)
            {
              /* BSIM4temp.rbsbx =  std::exp( std::log(model.BSIM4rbsbx0) + model.BSIM4rbsdbxl *  BSIM4temp.lnl +
                    model.BSIM4rbsdbxw *  BSIM4temp.lnw + model.BSIM4rbsdbxnf *  BSIM4temp.lnnf );
               BSIM4temp.rbsby =  std::exp( std::log(model.BSIM4rbsby0) + model.BSIM4rbsdbyl *  BSIM4temp.lnl +
                    model.BSIM4rbsdbyw *  BSIM4temp.lnw + model.BSIM4rbsdbynf *  BSIM4temp.lnnf );
                     */
                       BSIM4temp.rbsbx =  model.BSIM4rbsbx0 * std::exp( model.BSIM4rbsdbxl *  BSIM4temp.lnl +
                    model.BSIM4rbsdbxw *  BSIM4temp.lnw + model.BSIM4rbsdbxnf *  BSIM4temp.lnnf );
               BSIM4temp.rbsby =  model.BSIM4rbsby0 * std::exp( model.BSIM4rbsdbyl *  BSIM4temp.lnl +
                    model.BSIM4rbsdbyw *  BSIM4temp.lnw + model.BSIM4rbsdbynf *  BSIM4temp.lnnf );
              BSIM4instance.BSIM4rbsb =  BSIM4temp.rbsbx *  BSIM4temp.rbsby / ( BSIM4temp.rbsbx +  BSIM4temp.rbsby);


              /* BSIM4temp.rbdbx =  std::exp( std::log(model.BSIM4rbdbx0) + model.BSIM4rbsdbxl *  BSIM4temp.lnl +
                    model.BSIM4rbsdbxw *  BSIM4temp.lnw + model.BSIM4rbsdbxnf *  BSIM4temp.lnnf );
               BSIM4temp.rbdby =  std::exp( std::log(model.BSIM4rbdby0) + model.BSIM4rbsdbyl *  BSIM4temp.lnl +
                    model.BSIM4rbsdbyw *  BSIM4temp.lnw + model.BSIM4rbsdbynf *  BSIM4temp.lnnf );
                      */

                       BSIM4temp.rbdbx =  model.BSIM4rbdbx0 * std::exp( model.BSIM4rbsdbxl *  BSIM4temp.lnl +
                    model.BSIM4rbsdbxw *  BSIM4temp.lnw + model.BSIM4rbsdbxnf *  BSIM4temp.lnnf );
               BSIM4temp.rbdby =  model.BSIM4rbdby0 * std::exp(  model.BSIM4rbsdbyl *  BSIM4temp.lnl +
                    model.BSIM4rbsdbyw *  BSIM4temp.lnw + model.BSIM4rbsdbynf *  BSIM4temp.lnnf );

              BSIM4instance.BSIM4rbdb =  BSIM4temp.rbdbx *  BSIM4temp.rbdby / ( BSIM4temp.rbdbx +  BSIM4temp.rbdby);
            }

          if (( BSIM4temp.bodymode == 3)|| ( BSIM4temp.bodymode == 5))
            {
              /*BSIM4instance.BSIM4rbps = std::exp( std::log(model.BSIM4rbps0) + model.BSIM4rbpsl *  BSIM4temp.lnl +
                         model.BSIM4rbpsw *  BSIM4temp.lnw + model.BSIM4rbpsnf *  BSIM4temp.lnnf );
              BSIM4instance.BSIM4rbpd = std::exp( std::log(model.BSIM4rbpd0) + model.BSIM4rbpdl *  BSIM4temp.lnl +
                                     model.BSIM4rbpdw *  BSIM4temp.lnw + model.BSIM4rbpdnf *  BSIM4temp.lnnf );
                     */
                     BSIM4instance.BSIM4rbps = model.BSIM4rbps0 * std::exp( model.BSIM4rbpsl *  BSIM4temp.lnl +
                         model.BSIM4rbpsw *  BSIM4temp.lnw + model.BSIM4rbpsnf *  BSIM4temp.lnnf );
             BSIM4instance.BSIM4rbpd = model.BSIM4rbpd0 * std::exp( model.BSIM4rbpdl *  BSIM4temp.lnl +
                                     model.BSIM4rbpdw *  BSIM4temp.lnw + model.BSIM4rbpdnf *  BSIM4temp.lnnf );

            }

          /* BSIM4temp.rbpbx =  std::exp( std::log(model.BSIM4rbpbx0) + model.BSIM4rbpbxl *  BSIM4temp.lnl +
                model.BSIM4rbpbxw *  BSIM4temp.lnw + model.BSIM4rbpbxnf *  BSIM4temp.lnnf );
           BSIM4temp.rbpby =  std::exp( std::log(model.BSIM4rbpby0) + model.BSIM4rbpbyl *  BSIM4temp.lnl +
                model.BSIM4rbpbyw *  BSIM4temp.lnw + model.BSIM4rbpbynf *  BSIM4temp.lnnf );
                  */
                   BSIM4temp.rbpbx =  model.BSIM4rbpbx0 * std::exp(  model.BSIM4rbpbxl *  BSIM4temp.lnl +
                model.BSIM4rbpbxw *  BSIM4temp.lnw + model.BSIM4rbpbxnf *  BSIM4temp.lnnf );
           BSIM4temp.rbpby =  model.BSIM4rbpby0 * std::exp(  model.BSIM4rbpbyl *  BSIM4temp.lnl +
                model.BSIM4rbpbyw *  BSIM4temp.lnw + model.BSIM4rbpbynf *  BSIM4temp.lnnf );

          BSIM4instance.BSIM4rbpb =  BSIM4temp.rbpbx* BSIM4temp.rbpby/( BSIM4temp.rbpbx +  BSIM4temp.rbpby);
        }


              if ((BSIM4instance.BSIM4rbodyMod == 1 ) || ((BSIM4instance.BSIM4rbodyMod == 2 ) && ( BSIM4temp.bodymode == 5)) )
              {   if (BSIM4instance.BSIM4rbdb < 1.0e-3)
                      BSIM4instance.BSIM4grbdb = 1.0e3; /* in mho */
                  else
                      BSIM4instance.BSIM4grbdb = model.BSIM4gbmin + 1.0 / BSIM4instance.BSIM4rbdb;
                  if (BSIM4instance.BSIM4rbpb < 1.0e-3)
                      BSIM4instance.BSIM4grbpb = 1.0e3;
                  else
                      BSIM4instance.BSIM4grbpb = model.BSIM4gbmin + 1.0 / BSIM4instance.BSIM4rbpb;
                  if (BSIM4instance.BSIM4rbps < 1.0e-3)
                      BSIM4instance.BSIM4grbps = 1.0e3;
                  else
                      BSIM4instance.BSIM4grbps = model.BSIM4gbmin + 1.0 / BSIM4instance.BSIM4rbps;
                  if (BSIM4instance.BSIM4rbsb < 1.0e-3)
                      BSIM4instance.BSIM4grbsb = 1.0e3;
                  else
                      BSIM4instance.BSIM4grbsb = model.BSIM4gbmin + 1.0 / BSIM4instance.BSIM4rbsb;
                  if (BSIM4instance.BSIM4rbpd < 1.0e-3)
                      BSIM4instance.BSIM4grbpd = 1.0e3;
                  else
                      BSIM4instance.BSIM4grbpd = model.BSIM4gbmin + 1.0 / BSIM4instance.BSIM4rbpd;

              }

          if((BSIM4instance.BSIM4rbodyMod == 2) && ( BSIM4temp.bodymode == 3))
              {
                      BSIM4instance.BSIM4grbdb = BSIM4instance.BSIM4grbsb = model.BSIM4gbmin;
                  if (BSIM4instance.BSIM4rbpb < 1.0e-3)
                      BSIM4instance.BSIM4grbpb = 1.0e3;
                  else
                      BSIM4instance.BSIM4grbpb = model.BSIM4gbmin + 1.0 / BSIM4instance.BSIM4rbpb;
                  if (BSIM4instance.BSIM4rbps < 1.0e-3)
                      BSIM4instance.BSIM4grbps = 1.0e3;
                  else
                      BSIM4instance.BSIM4grbps = model.BSIM4gbmin + 1.0 / BSIM4instance.BSIM4rbps;
                  if (BSIM4instance.BSIM4rbpd < 1.0e-3)
                      BSIM4instance.BSIM4grbpd = 1.0e3;
                  else
                      BSIM4instance.BSIM4grbpd = model.BSIM4gbmin + 1.0 / BSIM4instance.BSIM4rbpd;
              }

          if((BSIM4instance.BSIM4rbodyMod == 2) && ( BSIM4temp.bodymode == 1))
              {
                      BSIM4instance.BSIM4grbdb = BSIM4instance.BSIM4grbsb = model.BSIM4gbmin;
              BSIM4instance.BSIM4grbps = BSIM4instance.BSIM4grbpd = 1.0e3;
                  if (BSIM4instance.BSIM4rbpb < 1.0e-3)
                      BSIM4instance.BSIM4grbpb = 1.0e3;
                  else
                      BSIM4instance.BSIM4grbpb = model.BSIM4gbmin + 1.0 / BSIM4instance.BSIM4rbpb;
              }


              /*
               * Process geomertry dependent parasitics
           */

              BSIM4instance.BSIM4grgeltd = model.BSIM4rshg * (BSIM4instance.BSIM4xgw
                      +BSIM4temp.pParam->BSIM4weffCJ / 3.0 / BSIM4instance.BSIM4ngcon) /
                      (BSIM4instance.BSIM4ngcon * BSIM4instance.BSIM4nf *
                      (BSIM4temp.Lnew - model.BSIM4xgl));
              if (BSIM4instance.BSIM4grgeltd > 0.0)
                  BSIM4instance.BSIM4grgeltd = 1.0 / BSIM4instance.BSIM4grgeltd;
              else
              {   BSIM4instance.BSIM4grgeltd = 1.0e3; /* mho */
          if (BSIM4instance.BSIM4rgateMod != 0)
                  printf("Warning: The gate conductance reset to 1.0e3 mho.\n");
              }

           BSIM4temp.DMCGeff = model.BSIM4dmcg - model.BSIM4dmcgt;
               BSIM4temp.DMCIeff = model.BSIM4dmci;
               BSIM4temp.DMDGeff = model.BSIM4dmdg - model.BSIM4dmcgt;

/*        if (BSIM4instance.BSIM4sourcePerimeterGiven)
          {   if (model.BSIM4perMod == 0)
                  BSIM4instance.BSIM4Pseff = BSIM4instance.BSIM4sourcePerimeter;
          else
              BSIM4instance.BSIM4Pseff = BSIM4instance.BSIM4sourcePerimeter
                       -BSIM4temp.pParam->BSIM4weffCJ * BSIM4instance.BSIM4nf;
          }
          else
              BSIM4PAeffGeo(BSIM4instance.BSIM4nf, BSIM4instance.BSIM4geoMod, BSIM4instance.BSIM4min,
                                   BSIM4temp.pParam->BSIM4weffCJ,  BSIM4temp.DMCGeff,  BSIM4temp.DMCIeff,  BSIM4temp.DMDGeff,
                    &(BSIM4instance.BSIM4Pseff), & BSIM4temp.dumPd, & BSIM4temp.dumAs, & BSIM4temp.dumAd);
              if (BSIM4instance.BSIM4Pseff < 0.0) /4.6.2/
              BSIM4instance.BSIM4Pseff = 0.0; */

    /* New Diode Model v4.7*/
          if (BSIM4instance.BSIM4sourcePerimeterGiven)
          {   /* given */
                if (BSIM4instance.BSIM4sourcePerimeter == 0.0)
                    BSIM4instance.BSIM4Pseff = 0.0;
                else if (BSIM4instance.BSIM4sourcePerimeter < 0.0)
                {
                    printf("Warning: Source Perimeter is specified as negative, it is set to zero.\n");
                    BSIM4instance.BSIM4Pseff = 0.0;
                } else
                {
                    if (model.BSIM4perMod == 0)
                        BSIM4instance.BSIM4Pseff = BSIM4instance.BSIM4sourcePerimeter;
                    else
                        BSIM4instance.BSIM4Pseff = BSIM4instance.BSIM4sourcePerimeter
                            -BSIM4temp.pParam->BSIM4weffCJ * BSIM4instance.BSIM4nf;
                }
          } else /* not given */
              BSIM4PAeffGeo(BSIM4instance.BSIM4nf, BSIM4instance.BSIM4geoMod, BSIM4instance.BSIM4min,
                                   BSIM4temp.pParam->BSIM4weffCJ,  BSIM4temp.DMCGeff,  BSIM4temp.DMCIeff,  BSIM4temp.DMDGeff,
                    &(BSIM4instance.BSIM4Pseff), & BSIM4temp.dumPd, & BSIM4temp.dumAs, & BSIM4temp.dumAd);

          if (BSIM4instance.BSIM4Pseff < 0.0){ /* v4.7 final check */
              BSIM4instance.BSIM4Pseff = 0.0;
              printf("Warning: Pseff is negative, it is set to zero.\n");
      }
            /*  if (BSIM4instance.BSIM4drainPerimeterGiven)
              {   if (model.BSIM4perMod == 0)
                      BSIM4instance.BSIM4Pdeff = BSIM4instance.BSIM4drainPerimeter;
                  else
                      BSIM4instance.BSIM4Pdeff = BSIM4instance.BSIM4drainPerimeter
                       -BSIM4temp.pParam->BSIM4weffCJ * BSIM4instance.BSIM4nf;
              }
              else
                  BSIM4PAeffGeo(BSIM4instance.BSIM4nf, BSIM4instance.BSIM4geoMod, BSIM4instance.BSIM4min,
                                   BSIM4temp.pParam->BSIM4weffCJ,  BSIM4temp.DMCGeff,  BSIM4temp.DMCIeff,  BSIM4temp.DMDGeff,
                    & BSIM4temp.dumPs, &(BSIM4instance.BSIM4Pdeff), & BSIM4temp.dumAs, & BSIM4temp.dumAd);
               if (BSIM4instance.BSIM4Pdeff < 0.0) /4.6.2/
              BSIM4instance.BSIM4Pdeff = 0.0; */

 if (BSIM4instance.BSIM4drainPerimeterGiven)
          {   /* given */
              if (BSIM4instance.BSIM4drainPerimeter == 0.0)
                    BSIM4instance.BSIM4Pdeff = 0.0;
              else if (BSIM4instance.BSIM4drainPerimeter < 0.0)
              {
                    printf("Warning: Drain Perimeter is specified as negative, it is set to zero.\n");
                    BSIM4instance.BSIM4Pdeff = 0.0;
              } else
              {
                    if (model.BSIM4perMod == 0)
                        BSIM4instance.BSIM4Pdeff = BSIM4instance.BSIM4drainPerimeter;
                    else
                        BSIM4instance.BSIM4Pdeff = BSIM4instance.BSIM4drainPerimeter
                            -BSIM4temp.pParam->BSIM4weffCJ * BSIM4instance.BSIM4nf;
              }
          } else /* not given */
              BSIM4PAeffGeo(BSIM4instance.BSIM4nf, BSIM4instance.BSIM4geoMod, BSIM4instance.BSIM4min,
                   BSIM4temp.pParam->BSIM4weffCJ,  BSIM4temp.DMCGeff,  BSIM4temp.DMCIeff,  BSIM4temp.DMDGeff,
                    & BSIM4temp.dumPs, &(BSIM4instance.BSIM4Pdeff), & BSIM4temp.dumAs, & BSIM4temp.dumAd);

          if (BSIM4instance.BSIM4Pdeff < 0.0){
              BSIM4instance.BSIM4Pdeff = 0.0; /*New Diode v4.7*/
              printf("Warning: Pdeff is negative, it is set to zero.\n");
          }
              if (BSIM4instance.BSIM4sourceAreaGiven)
                  BSIM4instance.BSIM4Aseff = BSIM4instance.BSIM4sourceArea;
              else
                  BSIM4PAeffGeo(BSIM4instance.BSIM4nf, BSIM4instance.BSIM4geoMod, BSIM4instance.BSIM4min,
                                   BSIM4temp.pParam->BSIM4weffCJ,  BSIM4temp.DMCGeff,  BSIM4temp.DMCIeff,  BSIM4temp.DMDGeff,
                    & BSIM4temp.dumPs, & BSIM4temp.dumPd, &(BSIM4instance.BSIM4Aseff), & BSIM4temp.dumAd);
          if (BSIM4instance.BSIM4Aseff < 0.0){
          BSIM4instance.BSIM4Aseff = 0.0; /* v4.7 */
              printf("Warning: Aseff is negative, it is set to zero.\n");
              }
              if (BSIM4instance.BSIM4drainAreaGiven)
                  BSIM4instance.BSIM4Adeff = BSIM4instance.BSIM4drainArea;
              else
                  BSIM4PAeffGeo(BSIM4instance.BSIM4nf, BSIM4instance.BSIM4geoMod, BSIM4instance.BSIM4min,
                                   BSIM4temp.pParam->BSIM4weffCJ,  BSIM4temp.DMCGeff,  BSIM4temp.DMCIeff,  BSIM4temp.DMDGeff,
                    & BSIM4temp.dumPs, & BSIM4temp.dumPd, & BSIM4temp.dumAs, &(BSIM4instance.BSIM4Adeff));
          if (BSIM4instance.BSIM4Adeff < 0.0){
          BSIM4instance.BSIM4Adeff = 0.0; /* v4.7 */
              printf("Warning: Adeff is negative, it is set to zero.\n");
              }
          /* Processing S/D resistance and conductance below */
              if(BSIM4instance.BSIM4sNodePrime != BSIM4instance.BSIM4sNode)
              {
                 BSIM4instance.BSIM4sourceConductance = 0.0;
                 if(BSIM4instance.BSIM4sourceSquaresGiven)
                 {
                    BSIM4instance.BSIM4sourceConductance = model.BSIM4sheetResistance
                                               * BSIM4instance.BSIM4sourceSquares;
                 } else if (BSIM4instance.BSIM4rgeoMod > 0)
                 {
                    BSIM4RdseffGeo(BSIM4instance.BSIM4nf, BSIM4instance.BSIM4geoMod,
                      BSIM4instance.BSIM4rgeoMod, BSIM4instance.BSIM4min,
                     BSIM4temp.pParam->BSIM4weffCJ, model.BSIM4sheetResistance,
                   BSIM4temp.DMCGeff,  BSIM4temp.DMCIeff,  BSIM4temp.DMDGeff, 1, &(BSIM4instance.BSIM4sourceConductance));
                 } else
                 {
                    BSIM4instance.BSIM4sourceConductance = 0.0;
                 }

                 if (BSIM4instance.BSIM4sourceConductance > 0.0)
                     BSIM4instance.BSIM4sourceConductance = 1.0
                                            / BSIM4instance.BSIM4sourceConductance;
                 else
                 {
                     BSIM4instance.BSIM4sourceConductance = 1.0e3; /* mho */
                     printf ("Warning: Source conductance reset to 1.0e3 mho.\n");
                 }
              } else
              {
                  BSIM4instance.BSIM4sourceConductance = 0.0;
              }

              if(BSIM4instance.BSIM4dNodePrime != BSIM4instance.BSIM4dNode)
              {
                 BSIM4instance.BSIM4drainConductance = 0.0;
                 if(BSIM4instance.BSIM4drainSquaresGiven)
                 {
                    BSIM4instance.BSIM4drainConductance = model.BSIM4sheetResistance
                                              * BSIM4instance.BSIM4drainSquares;
                 } else if (BSIM4instance.BSIM4rgeoMod > 0)
                 {
                    BSIM4RdseffGeo(BSIM4instance.BSIM4nf, BSIM4instance.BSIM4geoMod,
                      BSIM4instance.BSIM4rgeoMod, BSIM4instance.BSIM4min,
                     BSIM4temp.pParam->BSIM4weffCJ, model.BSIM4sheetResistance,
                   BSIM4temp.DMCGeff,  BSIM4temp.DMCIeff,  BSIM4temp.DMDGeff, 0, &(BSIM4instance.BSIM4drainConductance));
                 } else
                 {
                    BSIM4instance.BSIM4drainConductance = 0.0;
                 }

                 if (BSIM4instance.BSIM4drainConductance > 0.0)
                     BSIM4instance.BSIM4drainConductance = 1.0
                                           / BSIM4instance.BSIM4drainConductance;
                 else
                 {
                     BSIM4instance.BSIM4drainConductance = 1.0e3; /* mho */
                     printf ("Warning: Drain conductance reset to 1.0e3 mho.\n");
                  }
              } else
              {
                  BSIM4instance.BSIM4drainConductance = 0.0;
              }

               /* End of Rsd processing */


               BSIM4temp.Nvtms = model.BSIM4vtm * model.BSIM4SjctEmissionCoeff;
              if ((BSIM4instance.BSIM4Aseff <= 0.0) && (BSIM4instance.BSIM4Pseff <= 0.0))
              {    BSIM4temp.SourceSatCurrent = 0.0; /* v4.7 */
          /*  BSIM4temp.SourceSatCurrent = 1.0e-14; */
              }
              else
              {    BSIM4temp.SourceSatCurrent = BSIM4instance.BSIM4Aseff * model.BSIM4SjctTempSatCurDensity
                   + BSIM4instance.BSIM4Pseff * model.BSIM4SjctSidewallTempSatCurDensity
                                   +BSIM4temp.pParam->BSIM4weffCJ * BSIM4instance.BSIM4nf
                                   * model.BSIM4SjctGateSidewallTempSatCurDensity;
              }
              if ( BSIM4temp.SourceSatCurrent > 0.0)
              {   switch(model.BSIM4dioMod)
                  {   case 0:
              if ((model.BSIM4bvs /  BSIM4temp.Nvtms) > EXP_THRESHOLD)
                  BSIM4instance.BSIM4XExpBVS = model.BSIM4xjbvs * MIN_EXP;
              else
                          BSIM4instance.BSIM4XExpBVS = model.BSIM4xjbvs * std::exp(-model.BSIM4bvs /  BSIM4temp.Nvtms);
                  break;
                      case 1:
                          BSIM4DioIjthVjmEval( BSIM4temp.Nvtms, model.BSIM4ijthsfwd,  BSIM4temp.SourceSatCurrent,
                                  0.0, &(BSIM4instance.BSIM4vjsmFwd));
                          BSIM4instance.BSIM4IVjsmFwd =  BSIM4temp.SourceSatCurrent * std::exp(BSIM4instance.BSIM4vjsmFwd /  BSIM4temp.Nvtms);
                          break;
                      case 2:
                          if ((model.BSIM4bvs /  BSIM4temp.Nvtms) > EXP_THRESHOLD)
                          {   BSIM4instance.BSIM4XExpBVS = model.BSIM4xjbvs * MIN_EXP;
                   BSIM4temp.tmp = MIN_EXP;
              }
                          else
              {   BSIM4instance.BSIM4XExpBVS = std::exp(-model.BSIM4bvs /  BSIM4temp.Nvtms);
                   BSIM4temp.tmp = BSIM4instance.BSIM4XExpBVS;
                      BSIM4instance.BSIM4XExpBVS *= model.BSIM4xjbvs;
              }

                          BSIM4DioIjthVjmEval( BSIM4temp.Nvtms, model.BSIM4ijthsfwd,  BSIM4temp.SourceSatCurrent,
                                          BSIM4instance.BSIM4XExpBVS, &(BSIM4instance.BSIM4vjsmFwd));
                   BSIM4temp.T0 = std::exp(BSIM4instance.BSIM4vjsmFwd /  BSIM4temp.Nvtms);
                          BSIM4instance.BSIM4IVjsmFwd =  BSIM4temp.SourceSatCurrent * ( BSIM4temp.T0 - BSIM4instance.BSIM4XExpBVS /  BSIM4temp.T0
                          + BSIM4instance.BSIM4XExpBVS - 1.0);
                  BSIM4instance.BSIM4SslpFwd =  BSIM4temp.SourceSatCurrent
                           * ( BSIM4temp.T0 + BSIM4instance.BSIM4XExpBVS /  BSIM4temp.T0) /  BSIM4temp.Nvtms;

              BSIM4temp.T2 = model.BSIM4ijthsrev /  BSIM4temp.SourceSatCurrent;
              if (BSIM4temp.T2 < 1.0)
              {   BSIM4temp.T2 = 10.0;
                  fprintf(stderr, "Warning: ijthsrev too small and set to 10 times IsbSat.\n");
              }
                          BSIM4instance.BSIM4vjsmRev = -model.BSIM4bvs
                         -  BSIM4temp.Nvtms * std::log((BSIM4temp.T2 - 1.0) / model.BSIM4xjbvs);
              BSIM4temp.T1 = model.BSIM4xjbvs * std::exp(-(model.BSIM4bvs
                 + BSIM4instance.BSIM4vjsmRev) /  BSIM4temp.Nvtms);
              BSIM4instance.BSIM4IVjsmRev =  BSIM4temp.SourceSatCurrent * (1.0 + BSIM4temp.T1);
                          BSIM4instance.BSIM4SslpRev = - BSIM4temp.SourceSatCurrent * BSIM4temp.T1 /  BSIM4temp.Nvtms;
                          break;
                  default:
                          printf("Specified dioMod = %d not matched\n", model.BSIM4dioMod);
                  }
              }

               BSIM4temp.Nvtmd = model.BSIM4vtm * model.BSIM4DjctEmissionCoeff;
          if ((BSIM4instance.BSIM4Adeff <= 0.0) && (BSIM4instance.BSIM4Pdeff <= 0.0))
              {  /*  BSIM4temp.DrainSatCurrent = 1.0e-14;  v4.7 */
               BSIM4temp.DrainSatCurrent = 0.0;
              }
              else
              {    BSIM4temp.DrainSatCurrent = BSIM4instance.BSIM4Adeff * model.BSIM4DjctTempSatCurDensity
                  + BSIM4instance.BSIM4Pdeff * model.BSIM4DjctSidewallTempSatCurDensity
                                  +BSIM4temp.pParam->BSIM4weffCJ * BSIM4instance.BSIM4nf
                                  * model.BSIM4DjctGateSidewallTempSatCurDensity;
              }
              if ( BSIM4temp.DrainSatCurrent > 0.0)
              {   switch(model.BSIM4dioMod)
                  {   case 0:
                          if ((model.BSIM4bvd /  BSIM4temp.Nvtmd) > EXP_THRESHOLD)
                          BSIM4instance.BSIM4XExpBVD = model.BSIM4xjbvd * MIN_EXP;
                          else
                          BSIM4instance.BSIM4XExpBVD = model.BSIM4xjbvd * std::exp(-model.BSIM4bvd /  BSIM4temp.Nvtmd);
                          break;
                      case 1:
                          BSIM4DioIjthVjmEval( BSIM4temp.Nvtmd, model.BSIM4ijthdfwd,  BSIM4temp.DrainSatCurrent,
                                              0.0, &(BSIM4instance.BSIM4vjdmFwd));
                          BSIM4instance.BSIM4IVjdmFwd =  BSIM4temp.DrainSatCurrent * std::exp(BSIM4instance.BSIM4vjdmFwd /  BSIM4temp.Nvtmd);
                          break;
                      case 2:
                          if ((model.BSIM4bvd /  BSIM4temp.Nvtmd) > EXP_THRESHOLD)
                          {   BSIM4instance.BSIM4XExpBVD = model.BSIM4xjbvd * MIN_EXP;
                               BSIM4temp.tmp = MIN_EXP;
                          }
                          else
                          {   BSIM4instance.BSIM4XExpBVD = std::exp(-model.BSIM4bvd /  BSIM4temp.Nvtmd);
                               BSIM4temp.tmp = BSIM4instance.BSIM4XExpBVD;
                              BSIM4instance.BSIM4XExpBVD *= model.BSIM4xjbvd;
                          }

                          BSIM4DioIjthVjmEval( BSIM4temp.Nvtmd, model.BSIM4ijthdfwd,  BSIM4temp.DrainSatCurrent,
                                              BSIM4instance.BSIM4XExpBVD, &(BSIM4instance.BSIM4vjdmFwd));
                           BSIM4temp.T0 = std::exp(BSIM4instance.BSIM4vjdmFwd /  BSIM4temp.Nvtmd);
                          BSIM4instance.BSIM4IVjdmFwd =  BSIM4temp.DrainSatCurrent * ( BSIM4temp.T0 - BSIM4instance.BSIM4XExpBVD /  BSIM4temp.T0
                                              + BSIM4instance.BSIM4XExpBVD - 1.0);
                          BSIM4instance.BSIM4DslpFwd =  BSIM4temp.DrainSatCurrent
                                               * ( BSIM4temp.T0 + BSIM4instance.BSIM4XExpBVD /  BSIM4temp.T0) /  BSIM4temp.Nvtmd;

                          BSIM4temp.T2 = model.BSIM4ijthdrev /  BSIM4temp.DrainSatCurrent;
                          if (BSIM4temp.T2 < 1.0)
                          {   BSIM4temp.T2 = 10.0;
                              fprintf(stderr, "Warning: ijthdrev too small and set to 10 times IdbSat.\n");
                          }
                          BSIM4instance.BSIM4vjdmRev = -model.BSIM4bvd
                                             -  BSIM4temp.Nvtmd * std::log((BSIM4temp.T2 - 1.0) / model.BSIM4xjbvd); /* bugfix */
                          BSIM4temp.T1 = model.BSIM4xjbvd * std::exp(-(model.BSIM4bvd
                             + BSIM4instance.BSIM4vjdmRev) /  BSIM4temp.Nvtmd);
                          BSIM4instance.BSIM4IVjdmRev =  BSIM4temp.DrainSatCurrent * (1.0 + BSIM4temp.T1);
                          BSIM4instance.BSIM4DslpRev = - BSIM4temp.DrainSatCurrent * BSIM4temp.T1 /  BSIM4temp.Nvtmd;
                          break;
                  default:
                          printf("Specified dioMod = %d not matched\n", model.BSIM4dioMod);
                  }
              }

        /* GEDL current reverse bias */
             BSIM4temp.T0 = (BSIM4temp.TRatio - 1.0);
                model.BSIM4njtsstemp = model.BSIM4njts * (1.0 + model.BSIM4tnjts *  BSIM4temp.T0);
                model.BSIM4njtsswstemp = model.BSIM4njtssw * (1.0 + model.BSIM4tnjtssw *  BSIM4temp.T0);
                model.BSIM4njtsswgstemp = model.BSIM4njtsswg * (1.0 + model.BSIM4tnjtsswg *  BSIM4temp.T0);
                model.BSIM4njtsdtemp = model.BSIM4njtsd * (1.0 + model.BSIM4tnjtsd *  BSIM4temp.T0);
                model.BSIM4njtsswdtemp = model.BSIM4njtsswd * (1.0 + model.BSIM4tnjtsswd *  BSIM4temp.T0);
                model.BSIM4njtsswgdtemp = model.BSIM4njtsswgd * (1.0 + model.BSIM4tnjtsswgd *  BSIM4temp.T0);
                 BSIM4temp.T7 =  BSIM4temp.Eg0 / model.BSIM4vtm *  BSIM4temp.T0;
                 BSIM4temp.T9 = model.BSIM4xtss *  BSIM4temp.T7;
                 BSIM4temp.T1 = DEXP( BSIM4temp.T9);
                 BSIM4temp.T9 = model.BSIM4xtsd *  BSIM4temp.T7;
                 BSIM4temp.T2 = DEXP( BSIM4temp.T9);
                 BSIM4temp.T9 = model.BSIM4xtssws *  BSIM4temp.T7;
                 BSIM4temp.T3 = DEXP( BSIM4temp.T9);
                 BSIM4temp.T9 = model.BSIM4xtsswd *  BSIM4temp.T7;
                 BSIM4temp.T4 = DEXP( BSIM4temp.T9);
                 BSIM4temp.T9 = model.BSIM4xtsswgs *  BSIM4temp.T7;
                 BSIM4temp.T5 = DEXP( BSIM4temp.T9);
                 BSIM4temp.T9 = model.BSIM4xtsswgd *  BSIM4temp.T7;
                 BSIM4temp.T6 = DEXP( BSIM4temp.T9);
                /*IBM TAT*/
                if(model.BSIM4jtweff < 0.0)
                  {   model.BSIM4jtweff = 0.0;
                      fprintf(stderr, "TAT width dependence effect is negative. Jtweff is clamped to zero.\n");
                  }
                 BSIM4temp.T11 = std::sqrt(model.BSIM4jtweff /BSIM4temp.pParam->BSIM4weffCJ) + 1.0;

        BSIM4temp.T10 =BSIM4temp.pParam->BSIM4weffCJ * BSIM4instance.BSIM4nf;
        BSIM4instance.BSIM4SjctTempRevSatCur = BSIM4temp.T1 * BSIM4instance.BSIM4Aseff * model.BSIM4jtss;
        BSIM4instance.BSIM4DjctTempRevSatCur = BSIM4temp.T2 * BSIM4instance.BSIM4Adeff * model.BSIM4jtsd;
        BSIM4instance.BSIM4SswTempRevSatCur = BSIM4temp.T3 * BSIM4instance.BSIM4Pseff * model.BSIM4jtssws;
        BSIM4instance.BSIM4DswTempRevSatCur = BSIM4temp.T4 * BSIM4instance.BSIM4Pdeff * model.BSIM4jtsswd;
        BSIM4instance.BSIM4SswgTempRevSatCur = BSIM4temp.T5 * BSIM4temp.T10 *  BSIM4temp.T11 * model.BSIM4jtsswgs;
        BSIM4instance.BSIM4DswgTempRevSatCur = BSIM4temp.T6 * BSIM4temp.T10 *  BSIM4temp.T11 * model.BSIM4jtsswgd;

                if(model.BSIM4mtrlMod != 0 && model.BSIM4mtrlCompatMod == 0)
        {
            /* Calculate TOXP from EOT */
            /* Calculate  BSIM4temp.Vgs_eff @ Vgs = VDD with Poly Depletion Effect */
                     BSIM4temp.Vtm0eot = KboQ * model.BSIM4tempeot;
             BSIM4temp.Vtmeot  =  BSIM4temp.Vtm0eot;
             BSIM4temp.vbieot =  BSIM4temp.Vtm0eot * std::log(BSIM4temp.pParam->BSIM4nsd
                       *BSIM4temp.pParam->BSIM4ndep / ( BSIM4temp.ni *  BSIM4temp.ni));
             BSIM4temp.phieot =  BSIM4temp.Vtm0eot * std::log(BSIM4temp.pParam->BSIM4ndep /  BSIM4temp.ni)
                   +BSIM4temp.pParam->BSIM4phin + 0.4;
             BSIM4temp.tmp2 = BSIM4instance.BSIM4vfb +  BSIM4temp.phieot;
             BSIM4temp.vddeot = model.BSIM4type * model.BSIM4vddeot;
             BSIM4temp.T0 = model.BSIM4epsrgate * EPS0;
            if ((BSIM4temp.pParam->BSIM4ngate > 1.0e18) && (BSIM4temp.pParam->BSIM4ngate < 1.0e25)
            && ( BSIM4temp.vddeot >  BSIM4temp.tmp2) && ( BSIM4temp.T0!=0))
              {
            BSIM4temp.T1 = 1.0e6 * CHARGE *  BSIM4temp.T0 *BSIM4temp.pParam->BSIM4ngate /
              (model.BSIM4coxe * model.BSIM4coxe);
             BSIM4temp.T8 =  BSIM4temp.vddeot -  BSIM4temp.tmp2;
            BSIM4temp.T4 = std::sqrt(1.0 + 2.0 *  BSIM4temp.T8 / BSIM4temp.T1);
            BSIM4temp.T2 = 2.0 *  BSIM4temp.T8 / (BSIM4temp.T4 + 1.0);
            BSIM4temp.T3 = 0.5 * BSIM4temp.T2 * BSIM4temp.T2 / BSIM4temp.T1;
             BSIM4temp.T7 = 1.12 - BSIM4temp.T3 - 0.05;
            BSIM4temp.T6 = std::sqrt( BSIM4temp.T7 *  BSIM4temp.T7 + 0.224);
            BSIM4temp.T5 = 1.12 - 0.5 * ( BSIM4temp.T7 + BSIM4temp.T6);
             BSIM4temp.Vgs_eff =  BSIM4temp.vddeot - BSIM4temp.T5;
              }
            else
               BSIM4temp.Vgs_eff =  BSIM4temp.vddeot;

            /* Calculate  BSIM4temp.Vth @ Vds=Vbs=0 */

             BSIM4temp.V0 =  BSIM4temp.vbieot -  BSIM4temp.phieot;
             BSIM4temp.lt1 = model.BSIM4factor1*BSIM4temp.pParam->BSIM4sqrtXdep0;
             BSIM4temp.ltw =  BSIM4temp.lt1;
             BSIM4temp.T0 =BSIM4temp.pParam->BSIM4dvt1 * model.BSIM4leffeot /  BSIM4temp.lt1;
            if ( BSIM4temp.T0 < EXP_THRESHOLD)
              {
            BSIM4temp.T1 = std::exp( BSIM4temp.T0);
            BSIM4temp.T2 = BSIM4temp.T1 - 1.0;
            BSIM4temp.T3 = BSIM4temp.T2 * BSIM4temp.T2;
            BSIM4temp.T4 = BSIM4temp.T3 + 2.0 * BSIM4temp.T1 * MIN_EXP;
             BSIM4temp.Theta0 = BSIM4temp.T1 / BSIM4temp.T4;
              }
            else
               BSIM4temp.Theta0 = 1.0 / (MAX_EXP - 2.0);
             BSIM4temp.Delt_vth =BSIM4temp.pParam->BSIM4dvt0 *  BSIM4temp.Theta0 *  BSIM4temp.V0;
             BSIM4temp.T0 =BSIM4temp.pParam->BSIM4dvt1w * model.BSIM4weffeot * model.BSIM4leffeot /  BSIM4temp.ltw;
            if ( BSIM4temp.T0 < EXP_THRESHOLD)
              {   BSIM4temp.T1 = std::exp( BSIM4temp.T0);
              BSIM4temp.T2 = BSIM4temp.T1 - 1.0;
              BSIM4temp.T3 = BSIM4temp.T2 * BSIM4temp.T2;
              BSIM4temp.T4 = BSIM4temp.T3 + 2.0 * BSIM4temp.T1 * MIN_EXP;
              BSIM4temp.T5 = BSIM4temp.T1 / BSIM4temp.T4;
              }
            else
              BSIM4temp.T5 = 1.0 / (MAX_EXP - 2.0); /* 3.0 * MIN_EXP omitted */
            BSIM4temp.T2 =BSIM4temp.pParam->BSIM4dvt0w * BSIM4temp.T5 *  BSIM4temp.V0;
             BSIM4temp.TempRatioeot =  model.BSIM4tempeot / model.BSIM4tnom - 1.0;
             BSIM4temp.T0 = std::sqrt(1.0 +BSIM4temp.pParam->BSIM4lpe0 / model.BSIM4leffeot);
            BSIM4temp.T1 =BSIM4temp.pParam->BSIM4k1ox * ( BSIM4temp.T0 - 1.0) * std::sqrt( BSIM4temp.phieot)
              + (BSIM4temp.pParam->BSIM4kt1 +BSIM4temp.pParam->BSIM4kt1l / model.BSIM4leffeot) *  BSIM4temp.TempRatioeot;
             BSIM4temp.Vth_NarrowW =  BSIM4temp.toxe *  BSIM4temp.phieot
                  / (model.BSIM4weffeot +BSIM4temp.pParam->BSIM4w0);
             BSIM4temp.Lpe_Vb = std::sqrt(1.0 +BSIM4temp.pParam->BSIM4lpeb / model.BSIM4leffeot);
             BSIM4temp.Vth = model.BSIM4type * BSIM4instance.BSIM4vth0 +
              (BSIM4temp.pParam->BSIM4k1ox -BSIM4temp.pParam->BSIM4k1)*sqrt( BSIM4temp.phieot)* BSIM4temp.Lpe_Vb
              -  BSIM4temp.Delt_vth - BSIM4temp.T2 +BSIM4temp.pParam->BSIM4k3 *  BSIM4temp.Vth_NarrowW + BSIM4temp.T1;

            /* Calculate  BSIM4temp.n */
             BSIM4temp.tmp1 = BSIM4temp.epssub /BSIM4temp.pParam->BSIM4Xdep0;
             BSIM4temp.tmp2 =BSIM4temp.pParam->BSIM4nfactor *  BSIM4temp.tmp1;
             BSIM4temp.tmp3 = ( BSIM4temp.tmp2 +BSIM4temp.pParam->BSIM4cdsc *  BSIM4temp.Theta0 +BSIM4temp.pParam->BSIM4cit) / model.BSIM4coxe;
            if ( BSIM4temp.tmp3 >= -0.5)
               BSIM4temp.n = 1.0 +  BSIM4temp.tmp3;
            else
              {
             BSIM4temp.T0 = 1.0 / (3.0 + 8.0 *  BSIM4temp.tmp3);
             BSIM4temp.n = (1.0 + 3.0 *  BSIM4temp.tmp3) *  BSIM4temp.T0;
              }


            /*  BSIM4temp.Vth correction for Pocket implant */
            if (BSIM4temp.pParam->BSIM4dvtp0 > 0.0)
              {
            BSIM4temp.T3 = model.BSIM4leffeot +BSIM4temp.pParam->BSIM4dvtp0 * 2.0;
            if (model.BSIM4tempMod < 2)
              BSIM4temp.T4 =  BSIM4temp.Vtmeot * std::log(model.BSIM4leffeot / BSIM4temp.T3);
            else
              BSIM4temp.T4 =  BSIM4temp.Vtm0eot * std::log(model.BSIM4leffeot / BSIM4temp.T3);
             BSIM4temp.Vth -=  BSIM4temp.n * BSIM4temp.T4;
              }
              BSIM4temp.Vgsteff =  BSIM4temp.Vgs_eff- BSIM4temp.Vth;
            /* calculating Toxp */
            BSIM4temp.T3 = model.BSIM4type * BSIM4instance.BSIM4vth0
               - BSIM4instance.BSIM4vfb -  BSIM4temp.phieot;
            BSIM4temp.T4 = BSIM4temp.T3 + BSIM4temp.T3;
            BSIM4temp.T5 = 2.5 * BSIM4temp.T3;

              BSIM4temp.vtfbphi2eot = 4.0 * BSIM4temp.T3;
            if (  BSIM4temp.vtfbphi2eot < 0.0)
                  BSIM4temp.vtfbphi2eot = 0.0;


              BSIM4temp.niter = 0;
              BSIM4temp.toxpf =  BSIM4temp.toxe;
            do
              {
              BSIM4temp.toxpi =   BSIM4temp.toxpf;
             BSIM4temp.tmp2 = 2.0e8 *   BSIM4temp.toxpf;
             BSIM4temp.T0 = (  BSIM4temp.Vgsteff +   BSIM4temp.vtfbphi2eot) /  BSIM4temp.tmp2;
            BSIM4temp.T1 = 1.0 + std::exp(model.BSIM4bdos * 0.7 * std::log( BSIM4temp.T0));
              BSIM4temp.Tcen = model.BSIM4ados * 1.9e-9 / BSIM4temp.T1;
              BSIM4temp.toxpf =  BSIM4temp.toxe -   BSIM4temp.epsrox/model.BSIM4epsrsub *   BSIM4temp.Tcen;
              BSIM4temp.niter++;
              } while ((  BSIM4temp.niter<=4)&&(std::abs(BSIM4temp.toxpf -  BSIM4temp.toxpi)>1e-12));
              BSIM4instance.BSIM4toxp =   BSIM4temp.toxpf;
              BSIM4instance.BSIM4coxp =   BSIM4temp.epsrox * EPS0 / BSIM4instance.BSIM4toxp;


          } else {
          BSIM4instance.BSIM4toxp = model.BSIM4toxp;
          BSIM4instance.BSIM4coxp = model.BSIM4coxp;
              }

            if (BSIM4checkInstance(model, BSIM4instance, *BSIM4instance.pParam))
              {   
                  std::cerr << "Fatal error(s) detected during BSIM4 parameter checking for " 
                  << BSIM4instance.BSIM4name << " in model " << model.BSIM4modName << std::endl;
                  exit(1);
              }
  /* End instance */
  return 0;
}


}// namespace bsim4
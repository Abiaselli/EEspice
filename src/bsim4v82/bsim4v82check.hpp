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
#include <sstream>
#include <vector>
#include <iostream>
#include <cmath>
#include "bsim4v82.hpp"
namespace bsim4{
bool BSIM4checkModel(BSIM4model &model, double CKTtemp){
    // 0 = no error, 1 = error
    // Skip check if "modcheck" is true
    if(model.modcheck) {
        return false;
    }

    // Start model checking
    bool Fatal_Flag = false;

    // Use a vector of strings to accumulate messages.
    std::vector<std::string> messages;
    {
        std::ostringstream oss;
        oss << "\nChecking parameters for BSIM 4.8 model " << model.BSIM4modName << "\n";
        messages.push_back(oss.str());
    }
    // Check the version string (only allowed versions start with "4.8" or "4.81")
    if (!( model.BSIM4version == "4.8.2" ||
        model.BSIM4version.substr(0, 4) == "4.82" ) )
    {
        std::string warn1 = "Warning: This model supports BSIM4 version 4.8.2\n";
        std::string warn2 = "You specified a wrong version number. Working now with BSIM4.8.2\n";
        std::cout << warn1 << warn2;
        messages.push_back(warn1);
        messages.push_back(warn2);
    }
    

    // Begin parameter checks; if a check fails, append a message and mark fatal if needed.
    auto appendMsg = [&messages](std::string msg) {
        messages.push_back(std::move(msg));
    };

    if (model.BSIM4toxe <= 0.0 )
    {   
        std::ostringstream oss;
        oss << "Fatal: Toxe = " << model.BSIM4toxe << " is not positive.\n";
        appendMsg(oss.str());
        Fatal_Flag = true;
    }
    if (model.BSIM4toxp <= 0.0 )
    {
        std::ostringstream oss;
        oss << "Fatal: Toxp = " << model.BSIM4toxp << " is not positive.\n";
        appendMsg(oss.str());
        Fatal_Flag = true;
    }
    if (model.BSIM4eot <= 0.0 )
    {
        std::ostringstream oss;
        oss << "Fatal: EOT = " << model.BSIM4eot << " is not positive.\n";
        appendMsg(oss.str());
        Fatal_Flag = true;
    }
    if (model.BSIM4epsrgate < 0.0 )
    {
        std::ostringstream oss;
        oss << "Fatal: Epsrgate = " << model.BSIM4epsrgate << " is not positive.\n";
        appendMsg(oss.str());
        Fatal_Flag = true;
    }
    if (model.BSIM4epsrsub < 0.0 )
    {
        std::ostringstream oss;
        oss << "Fatal: Epsrsub = " << model.BSIM4epsrsub << " is not positive.\n";
        appendMsg(oss.str());
        Fatal_Flag = true;
    }
    if (model.BSIM4easub < 0.0 )
    {
        std::ostringstream oss;
        oss << "Fatal: Easub = " << model.BSIM4easub << " is not positive.\n";
        appendMsg(oss.str());
        Fatal_Flag = true;
    }
    if (model.BSIM4ni0sub <= 0.0 )
    {
        std::ostringstream oss;
        oss << "Fatal: ni0sub = " << model.BSIM4ni0sub << " is not positive.\n";
        appendMsg(oss.str());
        Fatal_Flag = true;
    }

    if (model.BSIM4toxm <= 0.0 )
    {
        std::ostringstream oss;
        oss << "Fatal: Toxm = " << model.BSIM4toxm << " is not positive.\n";
        appendMsg(oss.str());
        Fatal_Flag = true;
    }
    if (model.BSIM4toxref <= 0.0 )
    {
        std::ostringstream oss;
        oss << "Fatal: Toxref = " << model.BSIM4toxref << " is not positive.\n";
        appendMsg(oss.str());
        Fatal_Flag = true;
    }


    
    
    if (model.BSIM4gbmin < 1.0e-20 )
    {
        std::ostringstream oss;
        oss << "Warning: Gbmin = " << model.BSIM4gbmin << " is too small.\n";
        appendMsg(oss.str());
    }

    /* Check saturation parameters */
    if (model.BSIM4pditsl < 0.0 )
    {
        std::ostringstream oss;
        oss << "Fatal: pditsl = " << model.BSIM4pditsl << " is negative.\n";
        appendMsg(oss.str());
        Fatal_Flag = true;
    }

    

    if (model.BSIM4vtss < 0.0 )
    {
        std::ostringstream oss;
        oss << "Fatal: Vtss = " << model.BSIM4vtss << " is negative.\n";
        appendMsg(oss.str());
        Fatal_Flag = true;
    }
    if (model.BSIM4vtsd < 0.0 )
    {
        std::ostringstream oss;
        oss << "Fatal: Vtsd = " << model.BSIM4vtsd << " is negative.\n";
        appendMsg(oss.str());
        Fatal_Flag = true;
    }
    if (model.BSIM4vtssws < 0.0 )
    {
        std::ostringstream oss;
        oss << "Fatal: Vtssws = " << model.BSIM4vtssws << " is negative.\n";
        appendMsg(oss.str());
        Fatal_Flag = true;
    }
    if (model.BSIM4vtsswd < 0.0 )
    {
        std::ostringstream oss;
        oss << "Fatal: Vtsswd = " << model.BSIM4vtsswd << " is negative.\n";
        appendMsg(oss.str());
        Fatal_Flag = true;
    }
    if (model.BSIM4vtsswgs < 0.0 )
    {
        std::ostringstream oss;
        oss << "Fatal: Vtsswgs = " << model.BSIM4vtsswgs << " is negative.\n";
        appendMsg(oss.str());
        Fatal_Flag = true;
    }
    if (model.BSIM4vtsswgd < 0.0 )
    {
        std::ostringstream oss;
        oss << "Fatal: Vtsswgd = " << model.BSIM4vtsswgd << " is negative.\n";
        appendMsg(oss.str());
        Fatal_Flag = true;
    }


    if (model.BSIM4paramChk == 1)
    {

        /* Check threshold voltage parameters */
        if (model.BSIM4toxe < 1.0e-10)
        {
            std::ostringstream oss;
            oss << "Warning: Toxe = " << model.BSIM4toxe << " is less than 1A. Recommended Toxe >= 5A\n";
            appendMsg(oss.str());
        }
        if (model.BSIM4toxp < 1.0e-10)
        {
            std::ostringstream oss;
            oss << "Warning: Toxp = " << model.BSIM4toxp << " is less than 1A. Recommended Toxp >= 5A\n";
            appendMsg(oss.str());
        }
        if (model.BSIM4toxm < 1.0e-10)
        {
            std::ostringstream oss;
            oss << "Warning: Toxm = " << model.BSIM4toxm << " is less than 1A. Recommended Toxm >= 5A\n";
            appendMsg(oss.str());
        }
        /* Check GIDL parameters */
        if (model.BSIM4gidlclamp >= 0.0)
        {
            std::ostringstream oss;
            oss << "Warning: gidlclamp = " << model.BSIM4gidlclamp << " is zero or positive.\n";
            appendMsg(oss.str());
        }
    }

    /* Check body resistance parameters */
    if (model.BSIM4rbps0 <= 0.0)
    {   std::ostringstream oss;
        oss << "Fatal: RBPS0 = " << model.BSIM4rbps0 << " is not positive.\n";
        appendMsg(oss.str());
        Fatal_Flag = true;
    }
    if (model.BSIM4rbpd0 <= 0.0)
    {   std::ostringstream oss;
        oss << "Fatal: RBPD0 = " << model.BSIM4rbpd0 << " is not positive.\n";
        appendMsg(oss.str());
        Fatal_Flag = true;
    }
    if (model.BSIM4rbpbx0 <= 0.0)
    {   std::ostringstream oss;
        oss << "Fatal: RBPBX0 = " << model.BSIM4rbpbx0 << " is not positive.\n";
        appendMsg(oss.str());
        Fatal_Flag = true;
    }
    if (model.BSIM4rbpby0 <= 0.0)
    {   std::ostringstream oss;
        oss << "Fatal: RBPBY0 = " << model.BSIM4rbpby0 << " is not positive.\n";
        appendMsg(oss.str());
        Fatal_Flag = true;
    }
    if (model.BSIM4rbdbx0 <= 0.0)
    {   std::ostringstream oss;
        oss << "Fatal: RBDBX0 = " << model.BSIM4rbdbx0 << " is not positive.\n";
        appendMsg(oss.str());
        Fatal_Flag = true;
    }
    if (model.BSIM4rbdby0 <= 0.0)
    {   std::ostringstream oss;
        oss << "Fatal: RBDBY0 = " << model.BSIM4rbdby0 << " is not positive.\n";
        appendMsg(oss.str());
        Fatal_Flag = true;
    }
    if (model.BSIM4rbsbx0 <= 0.0)
    {   std::ostringstream oss;
        oss << "Fatal: RBSBX0 = " << model.BSIM4rbsbx0 << " is not positive.\n";
        appendMsg(oss.str());
        Fatal_Flag = true;
    }
    if (model.BSIM4rbsby0 <= 0.0)
    {   std::ostringstream oss;
        oss << "Fatal: RBSBY0 = " << model.BSIM4rbsby0 << " is not positive.\n";
        appendMsg(oss.str());
        Fatal_Flag = true;
    }

    /* Check overlap capacitance parameters */
    if (model.BSIM4cgdo < 0.0)
    {
        std::ostringstream oss;
        oss << "Warning: cgdo = " << model.BSIM4cgdo << " is negative. Set to zero.\n";
        appendMsg(oss.str());
        model.BSIM4cgdo = 0.0;
    }
    if (model.BSIM4cgso < 0.0)
    {
        std::ostringstream oss;
        oss << "Warning: cgso = " << model.BSIM4cgso << " is negative. Set to zero.\n";
        appendMsg(oss.str());
        model.BSIM4cgso = 0.0;
    }
    if (model.BSIM4cgbo < 0.0)
    {
        std::ostringstream oss;
        oss << "Warning: cgbo = " << model.BSIM4cgbo << " is negative. Set to zero.\n";
        appendMsg(oss.str());
        model.BSIM4cgbo = 0.0;
    }
    if (model.BSIM4tnoiMod == 1){
        std::ostringstream oss;
        oss << "Warning: TNOIMOD=1 is not supported and may be removed from future version.\n";
        appendMsg(oss.str());
    }
    if (model.BSIM4idovvdsc <= 0.0)
    {
        std::ostringstream oss;
        oss << "Warning: idovvdsc = " << model.BSIM4idovvdsc << " is zero or negative.\n";
        appendMsg(oss.str());
    }

    // if ((strcmp(model.BSIM4version, "4.8.1")) && (strncmp(model.BSIM4version, "4.81", 4)) &&
    //     (strcmp(model.BSIM4version, "4.8.2")) && (strncmp(model.BSIM4version, "4.82", 4)))
    // { /* checking for version <= 4.8 */
    //     /* v4.7 */
    //     if (model.BSIM4tnoiMod == 1 || model.BSIM4tnoiMod == 2) {
    //         if (model.BSIM4tnoia < 0.0)
    //         {
    //             std::ostringstream oss;
    //             oss << "Warning: tnoia = " << model.BSIM4tnoia << " is negative. Set to zero.\n";
    //             appendMsg(oss.str());
    //             model.BSIM4tnoia = 0.0;
    //         }
    //         if (model.BSIM4tnoib < 0.0)
    //         {
    //             std::ostringstream oss;
    //             oss << "Warning: tnoib = " << model.BSIM4tnoib << " is negative. Set to zero.\n";
    //             appendMsg(oss.str());
    //             model.BSIM4tnoib = 0.0;
    //         }

    //         if (model.BSIM4rnoia < 0.0)
    //         {
    //             std::ostringstream oss;
    //             oss << "Warning: rnoia = " << model.BSIM4rnoia << " is negative. Set to zero.\n";
    //             appendMsg(oss.str());
    //             model.BSIM4rnoia = 0.0;
    //         }
    //         if (model.BSIM4rnoib < 0.0)
    //         {
    //             std::ostringstream oss;
    //             oss << "Warning: rnoib = " << model.BSIM4rnoib << " is negative. Set to zero.\n";
    //             appendMsg(oss.str());
    //             model.BSIM4rnoib = 0.0;
    //         }
    //     }

    //     /* v4.7 */
    //     if (model.BSIM4tnoiMod == 2) {
    //         if (model.BSIM4tnoic < 0.0) {
    //             std::ostringstream oss;
    //             oss << "Warning: tnoic = " << model.BSIM4tnoic << " is negative. Set to zero.\n";
    //             appendMsg(oss.str());
    //             model.BSIM4tnoic = 0.0;
    //         }
    //         if (model.BSIM4rnoic < 0.0) {
    //             std::ostringstream oss;
    //             oss << "Warning: rnoic = " << model.BSIM4rnoic << " is negative. Set to zero.\n";
    //             appendMsg(oss.str());
    //             model.BSIM4rnoic = 0.0;
    //         }
    //     }
    // }
    // else
    // {
    // version = 4.8
        if (model.BSIM4tnoiMod == 1){
            if (model.BSIM4tnoia < 0.0) {
                std::ostringstream oss;
                oss << "Warning: tnoia = " << model.BSIM4tnoia << " is negative. Set to zero.\n";
                appendMsg(oss.str());
                model.BSIM4tnoia = 0.0;
            }
            if (model.BSIM4tnoib < 0.0) {
                std::ostringstream oss;
                oss << "Warning: tnoib = " << model.BSIM4tnoib << " is negative. Set to zero.\n";
                appendMsg(oss.str());
                model.BSIM4tnoib = 0.0;
            }
            if (model.BSIM4rnoia < 0.0) {
                std::ostringstream oss;
                oss << "Warning: rnoia = " << model.BSIM4rnoia << " is negative. Set to zero.\n";
                appendMsg(oss.str());
                model.BSIM4rnoia = 0.0;
            }
            if (model.BSIM4rnoib < 0.0) {
                std::ostringstream oss;
                oss << "Warning: rnoib = " << model.BSIM4rnoib << " is negative. Set to zero.\n";
                appendMsg(oss.str());
                model.BSIM4rnoib = 0.0;
            }
        }
    // }

    /* Limits of Njs and Njd modified in BSIM4.7 */
    if (model.BSIM4SjctEmissionCoeff < 0.1) {
        std::ostringstream oss;
        oss << "Warning: Njs = " << model.BSIM4SjctEmissionCoeff << " is less than 0.1. Setting Njs to 0.1.\n";
        appendMsg(oss.str());
        model.BSIM4SjctEmissionCoeff = 0.1;
    }
    else if (model.BSIM4SjctEmissionCoeff < 0.7) {
        std::ostringstream oss;
        oss << "Warning: Njs = " << model.BSIM4SjctEmissionCoeff << " is less than 0.7.\n";
        appendMsg(oss.str());
    }
    if (model.BSIM4DjctEmissionCoeff < 0.1)
    {
        std::ostringstream oss;
        oss << "Warning: Njd = " << model.BSIM4DjctEmissionCoeff << " is less than 0.1. Setting Njd to 0.1.\n";
        appendMsg(oss.str());
        model.BSIM4DjctEmissionCoeff = 0.1;
    }
    else if (model.BSIM4DjctEmissionCoeff < 0.7) {
        std::ostringstream oss;
        oss << "Warning: Njd = " << model.BSIM4DjctEmissionCoeff << " is less than 0.7.\n";
        appendMsg(oss.str());
    }

    if (model.BSIM4njtsstemp < 0.0) {
        std::ostringstream oss;
        oss << "Warning: Njts = " << model.BSIM4njtsstemp << " is negative at temperature = " << CKTtemp << ".\n";
        appendMsg(oss.str());
    }
    if (model.BSIM4njtsswstemp < 0.0) {
        std::ostringstream oss;
        oss << "Warning: Njtssw = " << model.BSIM4njtsswstemp << " is negative at temperature = " << CKTtemp << ".\n";
        appendMsg(oss.str());
    }
    if (model.BSIM4njtsswgstemp < 0.0) {
        std::ostringstream oss;
        oss << "Warning: Njtsswg = " << model.BSIM4njtsswgstemp << " is negative at temperature = " << CKTtemp << ".\n";
        appendMsg(oss.str());
    }

    if (model.BSIM4njtsdGiven && model.BSIM4njtsdtemp < 0.0)
    {
        std::ostringstream oss;
        oss << "Warning: Njtsd = " << model.BSIM4njtsdtemp << " is negative at temperature = " << CKTtemp << ".\n";
        appendMsg(oss.str());
    }
    if (model.BSIM4njtsswdGiven && model.BSIM4njtsswdtemp < 0.0)
    {
        std::ostringstream oss;
        oss << "Warning: Njtsswd = " << model.BSIM4njtsswdtemp << " is negative at temperature = " << CKTtemp << ".\n";
        appendMsg(oss.str());
    }
    if (model.BSIM4njtsswgdGiven && model.BSIM4njtsswgdtemp < 0.0)
    {
        std::ostringstream oss;
        oss << "Warning: Njtsswgd = " << model.BSIM4njtsswgdtemp << " is negative at temperature = " << CKTtemp << ".\n";
        appendMsg(oss.str());
    }

    if (model.BSIM4ntnoi < 0.0)
    {
        std::ostringstream oss;
        oss << "Warning: ntnoi = " << model.BSIM4ntnoi << " is negative. Set to zero.\n";
        appendMsg(oss.str());
        model.BSIM4ntnoi = 0.0;
    }

    /* diode model */
    if (model.BSIM4SbulkJctBotGradingCoeff >= 0.99)
    {
        std::ostringstream oss;
        oss << "Warning: MJS = " << model.BSIM4SbulkJctBotGradingCoeff << " is too big. Set to 0.99.\n";
        appendMsg(oss.str());
        model.BSIM4SbulkJctBotGradingCoeff = 0.99;
    }
    if (model.BSIM4SbulkJctSideGradingCoeff >= 0.99)
    {
        std::ostringstream oss;
        oss << "Warning: MJSWS = " << model.BSIM4SbulkJctSideGradingCoeff << " is too big. Set to 0.99.\n";
        appendMsg(oss.str());
        model.BSIM4SbulkJctSideGradingCoeff = 0.99;
    }
    if (model.BSIM4SbulkJctGateSideGradingCoeff >= 0.99)
    {
        std::ostringstream oss;
        oss << "Warning: MJSWGS = " << model.BSIM4SbulkJctGateSideGradingCoeff << " is too big. Set to 0.99.\n";
        appendMsg(oss.str());
        model.BSIM4SbulkJctGateSideGradingCoeff = 0.99;
    }

    if (model.BSIM4DbulkJctBotGradingCoeff >= 0.99)
    {
        std::ostringstream oss;
        oss << "Warning: MJD = " << model.BSIM4DbulkJctBotGradingCoeff << " is too big. Set to 0.99.\n";
        appendMsg(oss.str());
        model.BSIM4DbulkJctBotGradingCoeff = 0.99;
    }
    if (model.BSIM4DbulkJctSideGradingCoeff >= 0.99)
    {
        std::ostringstream oss;
        oss << "Warning: MJSWD = " << model.BSIM4DbulkJctSideGradingCoeff << " is too big. Set to 0.99.\n";
        appendMsg(oss.str());
        model.BSIM4DbulkJctSideGradingCoeff = 0.99;
    }
    if (model.BSIM4DbulkJctGateSideGradingCoeff >= 0.99)
    {
        std::ostringstream oss;
        oss << "Warning: MJSWGD = " << model.BSIM4DbulkJctGateSideGradingCoeff << " is too big. Set to 0.99.\n";
        appendMsg(oss.str());
        model.BSIM4DbulkJctGateSideGradingCoeff = 0.99;
    }
    if (model.BSIM4wpemod == 1)
    {
        if (model.BSIM4scref <= 0.0)
        {
            std::ostringstream oss;
            oss << "Warning: SCREF = " << model.BSIM4scref << " is not positive. Set to 1e-6.\n";
            appendMsg(oss.str());
            model.BSIM4scref = 1e-6;
        }
    }

    if (Fatal_Flag){
        for (auto &msg : messages){
            std::cout << msg << std::endl;
        }
    }

    model.modcheck = true;
    return Fatal_Flag;
}

bool BSIM4checkInstance(const BSIM4model &model, BSIM4V82 &instance, bsim4SizeDependParam &pParam){
    // 0 = no error, 1 = error

    // Start instance checking
    bool Fatal_Flag = false;

    // Use a vector of strings to accumulate messages.
    std::vector<std::string> messages;
    {
        std::ostringstream oss;
        oss << "\nChecking parameters for BSIM 4.8 instance " << instance.BSIM4name << "\n";
        messages.push_back(oss.str());
    }

    // Begin parameter checks; if a check fails, append a message and mark fatal if needed.
    auto appendMsg = [&messages](std::string msg) {
        messages.push_back(std::move(msg));
    };
    // Check for instance model pointer
    if (instance.BSIM4modPtr == nullptr) {
        appendMsg("Error: BSIM4modPtr is nullptr.");
        Fatal_Flag = true;
    }
    // Check for instance bsim4SizeDependParam pointer
    if (instance.pParam == nullptr) {
        appendMsg("Error: pParam is nullptr.");
        Fatal_Flag = true;
    }

    // Check for instance parameters and instance's bsim4SizeDependParam
    if ((instance.BSIM4rgateMod == 2) || (instance.BSIM4rgateMod == 3))
    {   if ((instance.BSIM4trnqsMod == 1) || (instance.BSIM4acnqsMod == 1)) {
            appendMsg("Warning: You've selected both Rg and charge deficit NQS; select one only.\n");
        }
    }
    if (pParam.BSIM4lpe0 < -pParam.BSIM4leff)
    {
        std::ostringstream oss;
        oss << "Fatal: Lpe0 = " << pParam.BSIM4lpe0 << " is less than -Leff.\n";
        appendMsg(oss.str());
        Fatal_Flag = true;
    }
    if (model.BSIM4lintnoi > pParam.BSIM4leff / 2.0)
    {
        std::ostringstream oss;
        oss << "Fatal: Lintnoi = " << model.BSIM4lintnoi << " is too large - Leff for noise is negative.\n";
        appendMsg(oss.str());
        Fatal_Flag = true;
    }
    if (pParam.BSIM4lpeb < -pParam.BSIM4leff)
    {
        std::ostringstream oss;
        oss << "Fatal: Lpeb = " << pParam.BSIM4lpeb << " is less than -Leff.\n";
        appendMsg(oss.str());
        Fatal_Flag = true;
    }
    if (pParam.BSIM4ndep <= 0.0)
    {
        std::ostringstream oss;
        oss << "Fatal: Ndep = " << pParam.BSIM4ndep << " is not positive.\n";
        appendMsg(oss.str());
        Fatal_Flag = true;
    }
    if (pParam.BSIM4phi <= 0.0)
    {
        std::ostringstream oss;
        oss << "Fatal: Phi = " << pParam.BSIM4phi << " is not positive. Please check Phin and Ndep\n" << std::endl;
        oss << "Phin = " << pParam.BSIM4phin << "  Ndep = " << pParam.BSIM4ndep << "\n";
        appendMsg(oss.str());
        Fatal_Flag = true;
    }
    if (pParam.BSIM4nsub <= 0.0)
    {
        std::ostringstream oss;
        oss << "Fatal: Nsub = " << pParam.BSIM4nsub << " is not positive.\n";
        appendMsg(oss.str());
        Fatal_Flag = true;
    }
    if (pParam.BSIM4ngate < 0.0)
    {
        std::ostringstream oss;
        oss << "Fatal: Ngate = " << pParam.BSIM4ngate << " is not positive.\n";
        appendMsg(oss.str());
        Fatal_Flag = true;
    }
    if (pParam.BSIM4ngate > 1.e25)
    {
        std::ostringstream oss;
        oss << "Fatal: Ngate = " << pParam.BSIM4ngate << " is too high.\n";
        appendMsg(oss.str());
        Fatal_Flag = true;
    }
    if (pParam.BSIM4xj <= 0.0)
    {
        std::ostringstream oss;
        oss << "Fatal: Xj = " << pParam.BSIM4xj << " is not positive.\n";
        appendMsg(oss.str());
        Fatal_Flag = true;
    }
    if (pParam.BSIM4dvt1 < 0.0)
    {
        std::ostringstream oss;
        oss << "Fatal: Dvt1 = " << pParam.BSIM4dvt1 << " is negative.\n";
        appendMsg(oss.str());
        Fatal_Flag = true;
    }
    if (pParam.BSIM4dvt1w < 0.0)
    {
        std::ostringstream oss;
        oss << "Fatal: Dvt1w = " << pParam.BSIM4dvt1w << " is negative.\n";
        appendMsg(oss.str());
        Fatal_Flag = true;
    }
    if (pParam.BSIM4w0 == -pParam.BSIM4weff)
    {
        std::ostringstream oss;
        oss << "Fatal: (W0 + Weff) = 0 causing divided-by-zero.\n";
        appendMsg(oss.str());
        Fatal_Flag = true;
    }
    if (pParam.BSIM4dsub < 0.0)
    {
        std::ostringstream oss;
        oss << "Fatal: Dsub = " << pParam.BSIM4dsub << " is negative.\n";
        appendMsg(oss.str());
        Fatal_Flag = true;
    }
    if (pParam.BSIM4b1 == -pParam.BSIM4weff)
    {
        std::ostringstream oss;
        oss << "Fatal: (B1 + Weff) = 0 causing divided-by-zero.\n";
        appendMsg(oss.str());
        Fatal_Flag = true;
    }
    if (instance.BSIM4u0temp <= 0.0)
    {
        std::ostringstream oss;
        oss << "Fatal: u0 at current temperature = " << instance.BSIM4u0temp << " is not positive.\n";
        appendMsg(oss.str());
        Fatal_Flag = true;
    }
    if (pParam.BSIM4delta < 0.0)
    {
        std::ostringstream oss;
        oss << "Fatal: Delta = " << pParam.BSIM4delta << " is less than zero.\n";
        appendMsg(oss.str());
        Fatal_Flag = true;
    }
    if (instance.BSIM4vsattemp <= 0.0)
    {
        std::ostringstream oss;
        oss << "Fatal: Vsat at current temperature = " << instance.BSIM4vsattemp << " is not positive.\n";
        appendMsg(oss.str());
        Fatal_Flag = true;
    }
    if (pParam.BSIM4pclm <= 0.0)
    {
        std::ostringstream oss;
        oss << "Fatal: Pclm = " << pParam.BSIM4pclm << " is not positive.\n";
        appendMsg(oss.str());
        Fatal_Flag = true;
    }
    if (pParam.BSIM4drout < 0.0)
    {
        std::ostringstream oss;
        oss << "Fatal: Drout = " << pParam.BSIM4drout << " is negative.\n";
        appendMsg(oss.str());
        Fatal_Flag = true;
    }
    // if (instance.BSIM4m <= 0.0)
    // {
    //     std::ostringstream oss;
    //     oss << "Fatal: multiplier = " << instance.BSIM4m << " is not positive.\n";
    //     appendMsg(oss.str());
    //     Fatal_Flag = true;
    // }
    if (instance.BSIM4nf < 1.0)
    {
        std::ostringstream oss;
        oss << "Fatal: Number of finger = " << instance.BSIM4nf << " is smaller than one.\n";
        appendMsg(oss.str());
        Fatal_Flag = true;
    }
    if ((instance.BSIM4sa > 0.0) && (instance.BSIM4sb > 0.0) &&
        ((instance.BSIM4nf == 1.0) || ((instance.BSIM4nf > 1.0) && (instance.BSIM4sd > 0.0))))
    {
        if (model.BSIM4saref <= 0.0 )
        {;
            std::ostringstream oss;
            oss << "Fatal: SAref = " << model.BSIM4saref << " is not positive.\n";
            appendMsg(oss.str());
            Fatal_Flag = true;
        }
        if (model.BSIM4sbref <= 0.0 )
        {
            std::ostringstream oss;
            oss << "Fatal: SBref = " << model.BSIM4sbref << " is not positive.\n";
            appendMsg(oss.str());
            Fatal_Flag = true;
        }
    }
    if ((instance.BSIM4l + model.BSIM4xl) <= model.BSIM4xgl)
    {
        std::ostringstream oss;
        oss << "Fatal: The parameter xgl must be smaller than Ldrawn+XL.\n";
        appendMsg(oss.str());
        Fatal_Flag = true;
    }
    if (instance.BSIM4ngcon < 1.0)
    {
        appendMsg("Fatal: The parameter ngcon cannot be smaller than one.\n");
        Fatal_Flag = true;
    }
    if ((instance.BSIM4ngcon != 1.0) && (instance.BSIM4ngcon != 2.0))
    {
        instance.BSIM4ngcon = 1.0;
        appendMsg("Warning: Ngcon must be equal to one or two; reset to 1.0.\n");
    }
    /* Check saturation parameters */
    if (pParam.BSIM4fprout < 0.0)
    {
        std::ostringstream oss;
        oss << "Fatal: fprout = " << pParam.BSIM4fprout << " is negative.\n";
        appendMsg(oss.str());
        Fatal_Flag = true;
    }
    if (pParam.BSIM4pdits < 0.0)
    {
        std::ostringstream oss;
        oss << "Fatal: pdits = " << pParam.BSIM4pdits << " is negative.\n";
        appendMsg(oss.str());
        Fatal_Flag = true;
    }
    /* Check gate current parameters */
    if (model.BSIM4igbMod) {
        if (pParam.BSIM4nigbinv <= 0.0)
        {
            std::ostringstream oss;
            oss << "Fatal: nigbinv = " << pParam.BSIM4nigbinv << " is non-positive.\n";
            appendMsg(oss.str());
            Fatal_Flag = true;
        }
        if (pParam.BSIM4nigbacc <= 0.0)
        {
            std::ostringstream oss;
            oss << "Fatal: nigbacc = " << pParam.BSIM4nigbacc << " is non-positive.\n";
            appendMsg(oss.str());
            Fatal_Flag = true;
        }
    }
    if (model.BSIM4igcMod) {
        if (pParam.BSIM4nigc <= 0.0)
        {
            std::ostringstream oss;
            oss << "Fatal: nigc = " << pParam.BSIM4nigc << " is non-positive.\n";
            appendMsg(oss.str());
            Fatal_Flag = true;
        }
        if (pParam.BSIM4poxedge <= 0.0)
        {
            std::ostringstream oss;
            oss << "Fatal: poxedge = " << pParam.BSIM4poxedge << " is non-positive.\n";
            appendMsg(oss.str());
            Fatal_Flag = true;
        }
        if (pParam.BSIM4pigcd <= 0.0)
        {
            std::ostringstream oss;
            oss << "Fatal: pigcd = " << pParam.BSIM4pigcd << " is non-positive.\n";
            appendMsg(oss.str());
            Fatal_Flag = true;
        }
    }
    /* Check capacitance parameters */
    if (pParam.BSIM4clc < 0.0)
    {
        std::ostringstream oss;
        oss << "Fatal: Clc = " << pParam.BSIM4clc << " is negative.\n";
        appendMsg(oss.str());
        Fatal_Flag = true;
    }
    /* Check overlap capacitance parameters */
    if (pParam.BSIM4ckappas < 0.02)
    {
        std::ostringstream oss;
        oss << "Warning: ckappas = " << pParam.BSIM4ckappas << " is too small.\n";
        appendMsg(oss.str());
        pParam.BSIM4ckappas = 0.02;
    }
    if (pParam.BSIM4ckappad < 0.02)
    {
        std::ostringstream oss;
        oss << "Warning: ckappad = " << pParam.BSIM4ckappad << " is too small.\n";
        appendMsg(oss.str());
        pParam.BSIM4ckappad = 0.02;
    }


    if (model.BSIM4paramChk == 1)
    {
        if (pParam.BSIM4ndep <= 1.0e12)
        {
            std::ostringstream oss;
            oss << "Warning: Ndep = " << pParam.BSIM4ndep << " may be too small.\n";
            appendMsg(oss.str());
        }
        else if (pParam.BSIM4ndep >= 1.0e21)
        {
            std::ostringstream oss;
            oss << "Warning: Ndep = " << pParam.BSIM4ndep << " may be too large.\n";
            appendMsg(oss.str());
        }

        if (pParam.BSIM4nsub <= 1.0e14)
        {
            std::ostringstream oss;
            oss << "Warning: Nsub = " << pParam.BSIM4nsub << " may be too small.\n";
            appendMsg(oss.str());
        }
        else if (pParam.BSIM4nsub >= 1.0e21)
        {
            std::ostringstream oss;
            oss << "Warning: Nsub = " << pParam.BSIM4nsub << " may be too large.\n";
            appendMsg(oss.str());
        }

        if ((pParam.BSIM4ngate > 0.0) &&
            (pParam.BSIM4ngate <= 1.e18))
        {
            std::ostringstream oss;
            oss << "Warning: Ngate = " << pParam.BSIM4ngate << " is less than 1.E18cm^-3.\n";
            appendMsg(oss.str());
        }
        if (pParam.BSIM4dvt0 < 0.0)
        {
            std::ostringstream oss;
            oss << "Warning: Dvt0 = " << pParam.BSIM4dvt0 << " is negative.\n";
            appendMsg(oss.str());
        }
        if (std::abs(1.0e-8 / (pParam.BSIM4w0 + pParam.BSIM4weff)) > 10.0)
        {
            std::string msg = "Warning: (W0 + Weff) may be too small.\n";
            appendMsg(msg);
        }
        /* Check subthreshold parameters */
        if (pParam.BSIM4nfactor < 0.0)
        {
            std::ostringstream oss;
            oss << "Warning: Nfactor = " << pParam.BSIM4nfactor << " is negative.\n";
            appendMsg(oss.str());
        }
        if (pParam.BSIM4cdsc < 0.0)
        {
            std::ostringstream oss;
            oss << "Warning: Cdsc = " << pParam.BSIM4cdsc << " is negative.\n";
            appendMsg(oss.str());
        }
        if (pParam.BSIM4cdscd < 0.0)
        {
            std::ostringstream oss;
            oss << "Warning: Cdscd = " << pParam.BSIM4cdscd << " is negative.\n";
            appendMsg(oss.str());
        }
        /* Check DIBL parameters */
        if (instance.BSIM4eta0 < 0.0)
        {
            std::ostringstream oss;
            oss << "Warning: Eta0 = " << instance.BSIM4eta0 << " is negative.\n";
            appendMsg(oss.str());
        }
        /* Check GIDL parameters */
        if (model.BSIM4gidlclamp >= 0.0)
        {
            std::ostringstream oss;
            oss << "Warning: gidlclamp = " << model.BSIM4gidlclamp << " is zero or positive.\n";
            appendMsg(oss.str());
        }
        /* Check Abulk parameters */
        if (std::abs(1.0e-8 / (pParam.BSIM4b1 + pParam.BSIM4weff)) > 10.0)
        {
            std::string msg = "Warning: (B1 + Weff) may be too small.\n";
            appendMsg(msg);
        }

        /* Check Saturation parameters */
        if (pParam.BSIM4a2 < 0.01)
        {
            std::ostringstream oss;
            oss << "Warning: A2 = " << pParam.BSIM4a2 << " is too small. Set to 0.01.\n";
            appendMsg(oss.str());
            pParam.BSIM4a2 = 0.01;
        }
        else if (pParam.BSIM4a2 > 1.0)
        {
            std::ostringstream oss;
            oss << "Warning: A2 = " << pParam.BSIM4a2 << " is larger than 1. A2 is set to 1 and A1 is set to 0.\n";
            appendMsg(oss.str());
            pParam.BSIM4a2 = 1.0;
            pParam.BSIM4a1 = 0.0;
        }

        if (pParam.BSIM4prwg < 0.0)
        {
            std::ostringstream oss;
            oss << "Warning: Prwg = " << pParam.BSIM4prwg << " is negative. Set to zero.\n";
            appendMsg(oss.str());
            pParam.BSIM4prwg = 0.0;
        }

        if (pParam.BSIM4rdsw < 0.0)
        {
            std::ostringstream oss;
            oss << "Warning: Rdsw = " << pParam.BSIM4rdsw << " is negative. Set to zero.\n";
            appendMsg(oss.str());
            pParam.BSIM4rdsw = 0.0;
            pParam.BSIM4rds0 = 0.0;
        }

        if (pParam.BSIM4rds0 < 0.0)
        {
            std::ostringstream oss;
            oss << "Warning: Rds at current temperature = " << pParam.BSIM4rds0 << " is negative. Set to zero.\n";
            appendMsg(oss.str());
            pParam.BSIM4rds0 = 0.0;
        }

        if (pParam.BSIM4rdswmin < 0.0)
        {
            std::ostringstream oss;
            oss << "Warning: Rdswmin at current temperature = " << pParam.BSIM4rdswmin << " is negative. Set to zero.\n";
            appendMsg(oss.str());
            pParam.BSIM4rdswmin = 0.0;
        }

        if (pParam.BSIM4pscbe2 <= 0.0)
        {
            std::ostringstream oss;
            oss << "Warning: Pscbe2 = " << pParam.BSIM4pscbe2 << " is not positive.\n";
            appendMsg(oss.str());
        }

        if (pParam.BSIM4vsattemp < 1.0e3)
        {
            std::ostringstream oss;
            oss << "Warning: Vsat at current temperature = " << pParam.BSIM4vsattemp << " may be too small.\n";
            appendMsg(oss.str());
        }

        if ((model.BSIM4lambdaGiven) && (pParam.BSIM4lambda > 0.0))
        {
            if (pParam.BSIM4lambda > 1.0e-9)
            {
                std::ostringstream oss;
                oss << "Warning: Lambda = " << pParam.BSIM4lambda << " is larger than 1.0e-9.\n";
                appendMsg(oss.str());
            }
        }

        if ((model.BSIM4vtlGiven) && (pParam.BSIM4vtl > 0.0))
        {
            if (pParam.BSIM4vtl < 6.0e4)
            {
                std::ostringstream oss;
                oss << "Warning: Thermal velocity vtl = " << pParam.BSIM4vtl << " ay be too small.\n";
                appendMsg(oss.str());
            }

            if (pParam.BSIM4xn < 3.0)
            {
                std::ostringstream oss;
                oss << "Warning: back scattering coeff xn = " << pParam.BSIM4xn << " is too small. Reset to 3.0 \n";
                appendMsg(oss.str());
                pParam.BSIM4xn = 3.0;
            }

            if (model.BSIM4lc < 0.0 )
            {
                std::ostringstream oss;
                oss << "Warning: back scattering coeff lc = " << model.BSIM4lc << " is too small. Reset to 0.0 \n";
                appendMsg(oss.str());
                pParam.BSIM4lc = 0.0;
            }
        }
        if (pParam.BSIM4pdibl1 < 0.0)
        {
            std::ostringstream oss;
            oss << "Warning: Pdibl1 = " << pParam.BSIM4pdibl1 << " is negative.\n";
            appendMsg(oss.str());
        }
    }

    if (pParam.BSIM4pdibl2 < 0.0)
    {
        std::ostringstream oss;
        oss << "Warning: Pdibl2 = " << pParam.BSIM4pdibl2 << " is negative.\n";
        appendMsg(oss.str());
    }

    /* Check stress effect parameters */
    if ((instance.BSIM4sa > 0.0) && (instance.BSIM4sb > 0.0) &&
        ((instance.BSIM4nf == 1.0) || ((instance.BSIM4nf > 1.0) && (instance.BSIM4sd > 0.0))) 
        )
    {
        if (model.BSIM4lodk2 <= 0.0)
        {
            std::ostringstream oss;
            oss << "Warning: LODK2 = " << model.BSIM4lodk2 << " is not positive.\n";
            appendMsg(oss.str());
        }
        if (model.BSIM4lodeta0 <= 0.0)
        {
            std::ostringstream oss;
            oss << "Warning: LODETA0 = " << model.BSIM4lodeta0 << " is not positive.\n";
            appendMsg(oss.str());
        }
    }

    /* Check gate resistance parameters */
    if (instance.BSIM4rgateMod == 1)
    {   if (model.BSIM4rshg <= 0.0 ){
            std::string msg = "Warning: rshg should be positive for rgateMod = 1.\n";
            appendMsg(msg);
        }
    }
    else if (instance.BSIM4rgateMod == 2)
    {   if (model.BSIM4rshg <= 0.0 ){
            std::string msg = "Warning: rshg <= 0.0 for rgateMod = 2.\n";
            appendMsg(msg);
        }
        else if (pParam.BSIM4xrcrg1 <= 0.0){
            std::string msg = "Warning: xrcrg1 <= 0.0 for rgateMod = 2.\n";
            appendMsg(msg);
        }
    }
    if (instance.BSIM4rgateMod == 3)
    {   if (model.BSIM4rshg <= 0.0 ){
            std::string msg = "Warning: rshg should be positive for rgateMod = 3.\n";
            appendMsg(msg);
        }
        else if (pParam.BSIM4xrcrg1 <= 0.0){
            std::string msg = "Warning: xrcrg1 should be positive for rgateMod = 3.\n";
            appendMsg(msg);
        }
    }
    /* Check capacitance parameters */
    if (pParam.BSIM4noff < 0.1)
    {
        std::ostringstream oss;
        oss << "Warning: Noff = " << pParam.BSIM4noff << " is too small.\n";
        appendMsg(oss.str());
    }

    if (pParam.BSIM4voffcv < -0.5)
    {
        std::ostringstream oss;
        oss << "Warning: Voffcv = " << pParam.BSIM4voffcv << " is too small.\n";
        appendMsg(oss.str());
    }

    if (pParam.BSIM4moin < 5.0)
    {
        std::ostringstream oss;
        oss << "Warning: Moin = " << pParam.BSIM4moin << " is too small.\n";
        appendMsg(oss.str());
    }
    if (pParam.BSIM4moin > 25.0)
    {
        std::ostringstream oss;
        oss << "Warning: Moin = " << pParam.BSIM4moin << " is too large.\n";
        appendMsg(oss.str());
    }
    if (model.BSIM4capMod == 2) {
        if (pParam.BSIM4acde < 0.1)
        {
            std::ostringstream oss;
            oss << "Warning: Acde = " << pParam.BSIM4acde << " is too small.\n";
            appendMsg(oss.str());
        }
        if (pParam.BSIM4acde > 1.6)
        {
            std::ostringstream oss;
            oss << "Warning: Acde = " << pParam.BSIM4acde << " is too large.\n";
            appendMsg(oss.str());
        }
    }

    if (model.BSIM4wpemod == 1){
        if (instance.BSIM4sca < 0.0)
        {
            std::ostringstream oss;
            oss << "Warning: SCA = " << instance.BSIM4sca << " is negative. Set to 0.0.\n";
            appendMsg(oss.str());
            instance.BSIM4sca = 0.0;
        }
        if (instance.BSIM4scb < 0.0)
        {
            std::ostringstream oss;
            oss << "Warning: SCB = " << instance.BSIM4scb << " is negative. Set to 0.0.\n";
            appendMsg(oss.str());
            instance.BSIM4scb = 0.0;
        }
        if (instance.BSIM4scc < 0.0)
        {
            std::ostringstream oss;
            oss << "Warning: SCC = " << instance.BSIM4scc << " is negative. Set to 0.0.\n";
            appendMsg(oss.str());
            instance.BSIM4scc = 0.0;
        }
        if (instance.BSIM4sc < 0.0)
        {
            std::ostringstream oss;
            oss << "Warning: SC = " << instance.BSIM4sc << " is negative. Set to 0.0.\n";
            appendMsg(oss.str());
            instance.BSIM4sc = 0.0;
        }
    }

    if (Fatal_Flag){
        for (auto &msg : messages){
            std::cout << msg << std::endl;
        }
    }

    return Fatal_Flag;
}
} // namespace bsim4
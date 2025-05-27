#pragma once
#include <cmath>
#include <iostream>
#include "AC_calcs.hpp"

// This function is used to calculate the AC sweep values and ACfreqDelta
void sweepCalculate(ACSweepSpec &sweepSpec){
   
    if(sweepSpec.fstart <= 0){
        std::cerr << "ERROR: AC startfreq <= 0" << std::endl;
        exit(1);
    }

    switch (sweepSpec.interval)
    {
    case ACSweepSpec::DEC:
        // Calculate for DEC
        if(sweepSpec.fstop / 10.0 < sweepSpec.fstart){
            /* start-stop frequencies less than a decade apart */
            if (sweepSpec.fstop == sweepSpec.fstart){
                sweepSpec.ACfreqDelta = 1;
            }
            else{
                sweepSpec.ACfreqDelta = std::exp(std::log(10.0) / sweepSpec.numpts);
            }
        }
        else{
            double num_steps = std::floor(std::fabs(std::log10(sweepSpec.fstop / sweepSpec.fstart)) * sweepSpec.numpts);
            sweepSpec.ACfreqDelta = std::exp((std::log(sweepSpec.fstop / sweepSpec.fstart)) / num_steps);
        }
        break;
    case ACSweepSpec::OCT:
        // Calculate for OCT
        sweepSpec.ACfreqDelta = std::exp(std::log(2.0) / sweepSpec.numpts);
        break;
    case ACSweepSpec::LIN:
        // Calculate for LIN
        if (sweepSpec.numpts - 1 > 1){
            sweepSpec.ACfreqDelta = (sweepSpec.fstop - sweepSpec.fstart) / (sweepSpec.numpts - 1);
        }
        else{
            sweepSpec.ACfreqDelta = 0;
        }
        break;
    default:
        std::cerr << "Invalid AC sweep interval type." << std::endl;
        exit(1);
        break;
    }
}
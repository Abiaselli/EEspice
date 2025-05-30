#pragma once
#include <cmath>
#include <iostream>
#include "AC_calcs.hpp"
#include "sim_variables.hpp"
#include "CKT.hpp"
#include "map.hpp"
#include "circuit_parser.hpp"

// This function is used to calculate the ACfreqDelta
double ACfreqDeltaCalculate(const ACSweepSpec &sweepSpec){
   
    if(sweepSpec.fstart <= 0){
        std::cerr << "ERROR: AC startfreq <= 0" << std::endl;
        exit(1);
    }

    double ACfreqDelta = 0.0;

    switch (sweepSpec.interval)
    {
    case ACSweepSpec::DEC:
        // Calculate for DEC
        if(sweepSpec.fstop / 10.0 < sweepSpec.fstart){
            /* start-stop frequencies less than a decade apart */
            if (sweepSpec.fstop == sweepSpec.fstart){
                ACfreqDelta = 1;
            }
            else{
                ACfreqDelta = std::exp(std::log(10.0) / sweepSpec.numpts);
            }
        }
        else{
            double num_steps = std::floor(std::fabs(std::log10(sweepSpec.fstop / sweepSpec.fstart)) * sweepSpec.numpts);
            ACfreqDelta = std::exp((std::log(sweepSpec.fstop / sweepSpec.fstart)) / num_steps);
        }
        break;
    case ACSweepSpec::OCT:
        // Calculate for OCT
        ACfreqDelta = std::exp(std::log(2.0) / sweepSpec.numpts);
        break;
    case ACSweepSpec::LIN:
        // Calculate for LIN
        if (sweepSpec.numpts - 1 > 1){
            ACfreqDelta = (sweepSpec.fstop - sweepSpec.fstart) / (sweepSpec.numpts - 1);
        }
        else{
            ACfreqDelta = 0;
        }
        break;
    default:
        std::cerr << "Invalid AC sweep interval type in ACfreqDeltaCalculate." << std::endl;
        exit(1);
    }
}

// This function calculates the frequency tolerance based on the sweep specification
// tolerence parameter for finding final frequency
double freqTolCalculate(const ACSweepSpec &sweepSpec){
    double freqTol = 0.0;
    switch (sweepSpec.interval)
    {
        case ACSweepSpec::DEC:
        case ACSweepSpec::OCT:
            // For DEC and OCT
            freqTol = sweepSpec.ACfreqDelta * sweepSpec.fstop * RELTOL;
            break;
        case ACSweepSpec::LIN:
            // For LIN,
            freqTol = sweepSpec.ACfreqDelta * RELTOL;
            break;
        default:
            std::cerr << "Invalid AC sweep interval type in freqTolCalculate." << std::endl;
            exit(1);
    }
    return freqTol;
}

// This function generates the sweep values for the AC sweep
std::vector<double> generateSweepValues(const ACSweepSpec &sweepSpec){
    std::vector<double> sweep_values;
    double freq = sweepSpec.fstart;
    sweep_values.push_back(freq);
    /* main loop through all scheduled frequencies */
    while (freq <= sweepSpec.fstop + sweepSpec.freqTol) {
        
        switch (sweepSpec.interval)
        {
            case ACSweepSpec::DEC:
            case ACSweepSpec::OCT:
                // For DEC and OCT
                freq *= sweepSpec.ACfreqDelta;
                if(sweepSpec.ACfreqDelta == 1){
                    // If ACfreqDelta is 1, we need to break to avoid infinite loop
                    return sweep_values;
                }
                break;
            case ACSweepSpec::LIN:
                // For LIN,
                freq += sweepSpec.ACfreqDelta;
                if(sweepSpec.ACfreqDelta == 0){
                    // If ACfreqDelta is 0, we need to break to avoid infinite loop
                    return sweep_values;
                }
                break;
            default:
                std::cerr << "Invalid AC sweep interval type in freqUpdate." << std::endl;
                exit(1);
        }

        sweep_values.push_back(freq);
    }
    return sweep_values;
}

ACsimulator ACsetup(CircuitParser &parser, const CKTcircuit &ckt){
    ACsimulator acsim;
    // Check if the AC sweep simulation is non-linear
    if(!ckt.CKTelements.nmos.empty() || !ckt.CKTelements.pmos.empty() || !ckt.CKTelements.diodes.empty()){
        acsim.non_linear = true;
    }
    else{
        acsim.non_linear = false;
    }
    // Set the AC sweep parameters
    acsim.acsweep = std::move(parser.acSweep_parser);
    acsim.acsweep.ACfreqDelta = ACfreqDeltaCalculate(acsim.acsweep);
    acsim.acsweep.freqTol = freqTolCalculate(acsim.acsweep);
    acsim.acsweep.sweep_values = generateSweepValues(acsim.acsweep);
    if(acsim.acsweep.sweep_values.empty()){
        std::cerr << "Error: No valid AC sweep values generated." << std::endl;
        exit(1);
    }
    return acsim;
}

std::vector<AC> AC_ops(CKTcircuit &ckt, ACsimulator &acSim, const Modelmap &modmap){

    
}

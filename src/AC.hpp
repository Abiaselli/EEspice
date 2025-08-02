#pragma once
#include <cmath>
#include <iostream>
#include "AC_calcs.hpp"
#include "sim_variables.hpp"
#include "CKT.hpp"
#include "map.hpp"
#include "circuit_parser.hpp"
#include "Newton.hpp"
#include "OP.hpp"
#include "simulation_exceptions.hpp"
#include "bsim4v82/bsim4v82acstamp.hpp"

namespace AC{

// This function is used to calculate the ACfreqDelta
double ACfreqDeltaCalculate(const ACSweepSpec &sweepSpec){
   
    if(sweepSpec.fstart <= 0){
        throw SetupException("ERROR: AC startfreq <= 0", "INVALID_AC_START_FREQUENCY");
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
        throw SetupException("Invalid AC sweep interval type in ACfreqDeltaCalculate.", "INVALID_AC_SWEEP_TYPE");
    }
    return ACfreqDelta;
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

ACsimulator ACsetup(const CircuitParser &parser, const CKTcircuit &ckt){
    ACsimulator acsim;
    // Check if the AC sweep simulation is non-linear
    acsim.non_linear = CKTisNonLinear(ckt.CKTelements);
    // Set the AC sweep parameters
    // acsim.acsweep = std::move(parser.acSweep_parser);
    acsim.acsweep = parser.acSweep_parser;
    acsim.acsweep.ACfreqDelta = ACfreqDeltaCalculate(acsim.acsweep);
    acsim.acsweep.freqTol = freqTolCalculate(acsim.acsweep);
    acsim.acsweep.sweep_values = generateSweepValues(acsim.acsweep);
    if(acsim.acsweep.sweep_values.empty()){
        std::cerr << "Error: No valid AC sweep values generated." << std::endl;
        exit(1);
    }
    // Reserve space for the AC results
    acsim.vec_ac.reserve(acsim.acsweep.sweep_values.size());

    return acsim;
}

// Update the parameters for non-linear devices and save in cktstate based on OP results
void SaveOP(CKTcircuit &ckt, const arma::vec &pre_NR_solution){
    // Make sure the cktstate is updated outside
    // Only BSIM4V82 is supported for now
    arma::mat init_LHS = ckt.cktdematrix->get_init_LHS();
    arma::vec init_RHS = ckt.cktdematrix->get_init_RHS();
    for (auto &nmos : ckt.CKTelements.nmos){
        if (nmos.modelType == MosfetModelType::BSIM4V82){
            const bsim4::BSIM4model &b4model = *nmos.bsim4v82Instance.BSIM4modPtr;
            bsim4::BSIM4V82 &b4instance = nmos.bsim4v82Instance;
            bsim4::BSIM4load(ckt, b4model, b4instance, ckt.spiceCompatible, pre_NR_solution, ckt.CKTtemp, ckt.CKTgmin, init_LHS, init_RHS);
        }
    }
    for (auto &pmos : ckt.CKTelements.pmos){
        if (pmos.modelType == MosfetModelType::BSIM4V82){
            const bsim4::BSIM4model &b4model = *pmos.bsim4v82Instance.BSIM4modPtr;
            bsim4::BSIM4V82 &b4instance = pmos.bsim4v82Instance;
            bsim4::BSIM4load(ckt, b4model, b4instance, ckt.spiceCompatible, pre_NR_solution, ckt.CKTtemp, ckt.CKTgmin, init_LHS, init_RHS);
        }
    }
}

void DynamicNonLinear(arma::cx_dmat &LHS, arma::cx_dvec &RHS, const CKTcircuit &ckt, double omega){
    for (auto &nmos : ckt.CKTelements.nmos){
        if (nmos.modelType == MosfetModelType::LEVEL1){

        }
        else if (nmos.modelType == MosfetModelType::BSIM4V82){
            const bsim4::BSIM4model &b4model = *nmos.bsim4v82Instance.BSIM4modPtr;
            const bsim4::BSIM4V82 &b4instance = nmos.bsim4v82Instance;
            bsim4::BSIM4acLoad(ckt, b4model, b4instance, ckt.spiceCompatible, omega, LHS, RHS);
        }
        else{
            std::cerr << "Error: NMOS model type is not supported!" << std::endl;
            exit(1);
        }
    }
    for (auto &pmos : ckt.CKTelements.pmos){
        if (pmos.modelType == MosfetModelType::LEVEL1){

        }
        else if (pmos.modelType == MosfetModelType::BSIM4V82){
            const bsim4::BSIM4model &b4model = *pmos.bsim4v82Instance.BSIM4modPtr;
            const bsim4::BSIM4V82 &b4instance = pmos.bsim4v82Instance;
            bsim4::BSIM4acLoad(ckt, b4model, b4instance, ckt.spiceCompatible, omega, LHS, RHS);
        }
        else{
            std::cerr << "Error: PMOS model type is not supported!" << std::endl;
            exit(1);
        }
    }
}

// AC simulation function
std::vector<AC> AC_ops(CKTcircuit &ckt, ACsimulator &acSim, const Modelmap &modmap){
    // 1. OP analysis
    ckt.spiceCompatible.setFlagsOP();
    arma::vec op_solution = OperatingPointAnalysis(ckt, modmap, acSim.non_linear);
    printOperatingPoint(op_solution, ckt);

    // 2. Update the small signal parameters for non-linear devices
    //    For BSIM4, we update the charge and capacitance into cktstate
    ckt.spiceCompatible.setFlagsSmallSig();
    if(acSim.non_linear){
        SaveOP(ckt, op_solution);
    }

    // 3. AC analysis
    for(const auto &freq : acSim.acsweep.sweep_values){
        ckt.spiceCompatible.setFlagsAC();
        AC ac;
        ac.freq = freq;
        ac.omega = 2 * M_PI * freq;

        // Update the LHS and RHS matrices
        ac.LHS = ckt.cktdematrix->get_init_cxLHS();
        ac.RHS = ckt.cktdematrix->get_init_cxRHS();
        DynamicNonLinear(ac.LHS, ac.RHS, ckt, ac.omega);

        // solve the MNA matrix
        ac.solution = arma::solve(ac.LHS, ac.RHS);

        // Store the AC result
        acSim.vec_ac.push_back(ac);
        
    }
    // 4. Return the AC results
    return acSim.vec_ac;
}

}// namespace AC

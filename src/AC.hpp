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

namespace ac{

// This function is used to calculate the ACfreqDelta
double ACfreqDeltaCalculate(const ACSweepSpec &sweepSpec){
   
    if(sweepSpec.fstart <= 0){
        throw SetupException("ERROR: AC startfreq <= 0", "ACfreqDeltaCalculate");
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
        throw SetupException("Invalid AC sweep interval type in ACfreqDeltaCalculate.", "ACfreqDeltaCalculate");
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
            throw SetupException("Invalid AC sweep interval type.", "freqTolCalculate");
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
                throw SetupException("Invalid AC sweep interval type.", "generateSweepValues");
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
        throw SetupException("No valid AC sweep values generated.", "ACsetup");
    }
    // Reserve space for the AC results
    acsim.vec_ac.reserve(acsim.acsweep.sweep_values.size());

    return acsim;
}

// Update the parameters for non-linear devices and save in cktstate based on OP results
void SaveOP(CKTcircuit &ckt, const arma::vec &pre_NR_solution){
    // Make sure the cktstate is updated outside
    // Only BSIM4V82 is supported for now
    HybridMatrix init_LHS = ckt.cktdematrix->get_init_LHS();
    arma::vec init_RHS = ckt.cktdematrix->get_init_RHS();

    for (auto &bsim4 : ckt.CKTelements.bsim4){
        const bsim4::BSIM4model &b4model = *bsim4.bsim4v82Instance.BSIM4modPtr;
        bsim4::BSIM4V82 &b4instance = bsim4.bsim4v82Instance;
        bsim4::BSIM4load(ckt, b4model, b4instance, ckt.spiceCompatible, pre_NR_solution, ckt.CKTtemp, ckt.CKTgmin, init_LHS, init_RHS);
    }
}

void DynamicNonLinear(arma::cx_dmat &LHS, arma::cx_dvec &RHS, const CKTcircuit &ckt, double omega, SimulationTime &simTime){
    ScopedTimer loadTimer(simTime.matrix_load_time);
    for (auto &nmos : ckt.CKTelements.nmos){
        throw SetupException("NMOS model type LEVEL1 is not supported in AC analysis.", "AC::DynamicNonLinear");
    }
    for (auto &pmos : ckt.CKTelements.pmos){
        throw SetupException("PMOS model type LEVEL1 is not supported in AC analysis.", "AC::DynamicNonLinear");
    }
    for (auto &bsim4 : ckt.CKTelements.bsim4){
        const bsim4::BSIM4model &b4model = *bsim4.bsim4v82Instance.BSIM4modPtr;
        const bsim4::BSIM4V82 &b4instance = bsim4.bsim4v82Instance;
        bsim4::BSIM4acLoad(ckt, b4model, b4instance, ckt.spiceCompatible, omega, LHS, RHS);
    }
    for (auto &cap : ckt.CKTelements.capacitors){
        C_ACassigner(cap.nodePos, cap.nodeNeg, cap.value, omega, LHS);
    }
}

// AC simulation function
std::vector<ACResult> AC_ops(CKTcircuit &ckt, ACsimulator &acSim, const Modelmap &modmap){
    ScopedTimer analysisTimer(ckt.sim_stats.simTime.analysis_time); // Time the analysis
    // 1. OP analysis
    ckt.spiceCompatible.setFlagsOP();
    arma::vec op_solution = OperatingPointAnalysis(ckt, modmap, acSim.non_linear);
    if (!batchMode){
        printOperatingPointWithNames(op_solution, ckt.map);
    }

    // 2. Update the small signal parameters for non-linear devices
    //    For BSIM4, we update the charge and capacitance into cktstate
    ckt.spiceCompatible.setFlagsSmallSig();
    if(acSim.non_linear){
        SaveOP(ckt, op_solution);
    }

    // 3. AC analysis
    // Initial complex matrices for whole AC simulation
    ACMat acMat;

    // AC simulation loop
    for(const auto &freq : acSim.acsweep.sweep_values){
        ckt.spiceCompatible.setFlagsAC();
        ACResult ac;
        ac.freq = freq;
        ac.omega = 2 * M_PI * freq;

        // Update the LHS and RHS matrices
        acMat.LHS = ckt.cktdematrix->get_init_cxLHS();
        acMat.RHS = ckt.cktdematrix->get_init_cxRHS();
        DynamicNonLinear(acMat.LHS, acMat.RHS, ckt, ac.omega, ckt.sim_stats.simTime);

        // solve the MNA matrix
        ac.solution = arma::solve(acMat.LHS, acMat.RHS);

        // Store the AC result
        acSim.vec_ac.push_back(ac);
        
    }

    ckt.sim_stats.num_data_points = static_cast<int>(acSim.vec_ac.size());
    
    // 4. Return the AC results
    return acSim.vec_ac;
}

std::vector<ACResultPolar> convertToPolar(const std::vector<ACResult> &vec_ac){
    std::vector<ACResultPolar> polarResults;
    polarResults.reserve(vec_ac.size());
    constexpr double coff = 180.0 / M_PI;

    // Loop through the results for each frequency
    for (const auto &ac : vec_ac) {
        polarResults.emplace_back( ac.freq, 
                                    arma::abs(ac.solution),        // Magnitude
                                    arma::arg(ac.solution) * coff  // Phase in degrees
        );
    }
    return polarResults;
}

}// namespace ac

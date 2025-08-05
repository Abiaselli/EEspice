#pragma once
#include <vector>
#include <armadillo>

namespace ac {

struct ACSweepSpec{
    int interval{};                     // decade (DEC), octave (OCT) or linearly (LIN)
    double numpts{};                    // the number of frequency points used per interval
    double fstart{};                    // start frequency
    double fstop{};                     // stop frequency
    std::vector<double> sweep_values;   // All sweep values from fstart to fstop

    double ACfreqDelta{};               // Multiplication factor for an AC analysis
    double freqTol{};                   // tolerence parameter for finding final frequency

    enum ACSweepType : int {
        DEC = 1,    // Decade
        OCT,        // Octave
        LIN         // Linear
    };
};

struct ACMat{
    arma::cx_dmat LHS;
    arma::cx_dvec RHS;  
};

struct ACResult{
    double freq{};                // Frequency of the AC point
    double omega{};               // Angular frequency (2 * pi * freq)          
    arma::cx_dvec solution;        
};

// A new struct to store the results in polar form (Magnitude and Phase)
struct ACResultPolar {
    double freq{};
    arma::vec magnitudes;
    arma::vec phases_rad; // Phase in radians
};

enum class ACResultType{
    MagPhase, // Magnitude and Phase
    RealImag // Real and Imaginary parts
};

struct ACsimulator {
    ACSweepSpec acsweep;                 // sweep device data
    std::vector<ACResult> vec_ac;        // All AC results
    bool non_linear = false;             // Non-linear flag
    ACResultType type = ACResultType::MagPhase;
};

} // namespace ac

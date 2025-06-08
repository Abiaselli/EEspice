#pragma once
#include <vector>
#include <armadillo>

namespace AC {

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

struct AC{
    double freq{};                // Frequency of the AC point
    double omega{};               // Angular frequency (2 * pi * freq)
    arma::cx_dmat LHS;
    arma::cx_dvec RHS;            
    arma::cx_dvec solution;        

};

struct ACsimulator {
    ACSweepSpec acsweep;                 // sweep device data
    std::vector<AC> vec_ac;              // All AC results
    bool non_linear = false;             // Non-linear flag
};

} // namespace AC

#pragma once
#include <vector>
#include <armadillo>

struct ACSweepSpec{
    int interval{};                     // decade (DEC), octave (OCT) or linearly (LIN)
    double numpts{};                    // the number of frequency points used per interval
    double fstart{};                    // start frequency
    double fstop{};                     // stop frequency
    std::vector<double> sweep_values;   // All sweep values from fstart to fstop

    double ACfreqDelta{};               //Multiplication factor for an AC analysis

    enum ACSweepType : int {
        DEC = 1,    // Decade
        OCT,        // Octave
        LIN         // Linear
    };
};

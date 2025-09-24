#pragma once
#include "XB_timer.hpp"


// causes accounting and run time statistics
struct SimulationStatistics{
    SimulationTime simTime;
    int total_NR_iteration{}; // Total Newton Raphson iterations in the entire simulation
    int MNA_Matrix_size{};    // Size of the MNA matrix
};
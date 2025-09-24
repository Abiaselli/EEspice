#pragma once
#include "XB_timer.hpp"
#include <iostream>
#include <iomanip>
#include <chrono>
#include <string>


// Helper function to format durations
template<typename Duration>
std::string formatDuration(Duration d) {
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(d).count();
    auto us = std::chrono::duration_cast<std::chrono::microseconds>(d).count() % 1000;
    return std::to_string(ms) + "." + std::to_string(us) + " ms";
}

// causes accounting and run time statistics
struct SimulationStatistics{
    SimulationTime simTime;
    int total_NR_iteration{}; // Total Newton Raphson iterations in the entire simulation
    int MNA_Matrix_size{};    // Size of the MNA matrix

    void printStatistics() const;
};

inline void SimulationStatistics::printStatistics() const {
    std::cout << "\n========================================\n";
    std::cout << "Simulation Statistics\n";
    std::cout << "========================================\n";

    // Timing information
    std::cout << "Timing Breakdown:\n";
    std::cout << "  Total time:         " << formatDuration(simTime.total_time.total()) << "\n";
    std::cout << "  Parse time:         " << formatDuration(simTime.parse_time.total()) << "\n";
    std::cout << "  Setup time:         " << formatDuration(simTime.setup_time.total()) << "\n";
    std::cout << "  Analysis time:      " << formatDuration(simTime.analysis_time.total()) << "\n";
    std::cout << "  Matrix loading:     " << formatDuration(simTime.matrix_load_time.total()) << "\n";
    std::cout << "  Solver time:        " << formatDuration(simTime.solve_time.total()) << "\n";
    std::cout << "  Newton method:      " << formatDuration(simTime.newton_time.total()) << "\n";

    // Iteration and matrix information
    std::cout << "\nSimulation Details:\n";
    std::cout << "  Total NR iterations: " << total_NR_iteration << "\n";
    std::cout << "  MNA matrix size:     " << MNA_Matrix_size << "\n";
    std::cout << "========================================\n";
}
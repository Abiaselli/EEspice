#pragma once
#include "XB_timer.hpp"
#include <iostream>
#include <iomanip>
#include <chrono>
#include <string>
#include <format>


// Helper function to format durations with automatic unit selection
template<typename Duration>
std::string formatDuration(Duration d) {
    auto total_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(d).count();
    
    // Use nanoseconds for durations < 1 µs
    if (total_ns < 1'000) {
        return std::format("{:.3f} ns", static_cast<double>(total_ns));
    }
    // Use milliseconds for durations < 1 s
    if (total_ns < 1'000'000'000) {
        auto ms = total_ns / 1'000'000.0;
        return std::format("{:.3f} ms", ms);
    }
    // Use seconds otherwise
    else {
        auto s = total_ns / 1'000'000'000.0;
        return std::format("{:.3f} s", s);
    }
}

// causes accounting and run time statistics
struct SimulationStatistics{
    SimulationTime simTime;
    int total_NR_iteration{}; // Total Newton Raphson iterations in the entire simulation
    int MNA_Matrix_size{};    // Size of the MNA matrix
    int num_data_points{};   // Number of data points in the simulation (AC, DC sweep, Transient)

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
    std::cout << "  No. Data points:     " << num_data_points << "\n";
    std::cout << "========================================\n";
}
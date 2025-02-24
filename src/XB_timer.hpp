#pragma once

#include <chrono>

struct XB_Timer;
struct SimulationTime;

struct XB_Timer
{
public:
    // Returns the current time in milliseconds
    double current_ms() const
    {
        return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start_time).count();
    }
    // Starts the timer
    void start()
    {
        start_time = std::chrono::steady_clock::now();
    }
    // Stops the timer and calculates the elapsed time
    void stop()
    {
        elapsed_time = std::chrono::steady_clock::now() - start_time;
    }
    // Returns the elapsed time in milliseconds
    double ms() const
    {
        return std::chrono::duration<double, std::milli>(elapsed_time).count();
    }
    // Adds the elapsed time to the total time
    void total()
    {
        total_time += ms();
    }
    // Returns the total time in milliseconds
    double total_ms() const
    {
        return total_time;
    }

private:
    std::chrono::time_point<std::chrono::steady_clock> start_time = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_time = std::chrono::duration<double>::zero();
    double total_time{};
};

struct SimulationTime{
    XB_Timer total_time;
    XB_Timer parse_time;
    XB_Timer setup_time;
    XB_Timer analysis_time;
    XB_Timer solve_time;
    XB_Timer newton_time;
};

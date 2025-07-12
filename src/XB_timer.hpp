#pragma once

#include <chrono>

struct XB_Timer;
struct SimulationTime;

struct XB_Timer
{
public:
    // Starts or restarts the timer.
    void start()
    {
        start_time = std::chrono::steady_clock::now();
        is_running = true;
    }

    // Stops the timer, calculates the last elapsed duration, and adds it to the total.
    void stop()
    {
        if (is_running)
        {
            last_elapsed_time = std::chrono::steady_clock::now() - start_time;
            total_duration += last_elapsed_time;
            is_running = false;
        }
    }

    // Returns the duration of the last interval (from start() to stop()) in milliseconds.
    double last_elapsed_ms() const
    {
        return std::chrono::duration<double, std::milli>(last_elapsed_time).count();
    }

    // Returns the total accumulated time in milliseconds.
    double total_ms() const
    {
        return std::chrono::duration<double, std::milli>(total_duration).count();
    }

    // Resets the total accumulated time to zero.
    void reset()
    {
        total_duration = std::chrono::duration<double>::zero();
        last_elapsed_time = std::chrono::duration<double>::zero();
    }

private:
    std::chrono::time_point<std::chrono::steady_clock> start_time;
    std::chrono::duration<double> last_elapsed_time{};
    std::chrono::duration<double> total_duration{};
    bool is_running = false;
};

struct SimulationTime{
    XB_Timer total_time;
    XB_Timer parse_time;
    XB_Timer setup_time;
    XB_Timer analysis_time;
    XB_Timer solve_time;
    XB_Timer newton_time;
};

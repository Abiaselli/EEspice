#pragma once

#include <chrono>

struct XB_Timer;

struct XB_Timer
{
public:
    double current_ms() const
    {
        return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start_time).count();
    }

    void start()
    {
        start_time = std::chrono::steady_clock::now();
    }

    void stop()
    {
        elapsed_time = std::chrono::steady_clock::now() - start_time;
    }

    double ms() const
    {
        return std::chrono::duration<double, std::milli>(elapsed_time).count();
    }

    void total()
    {
        total_time += ms();
    }

    double total_ms() const
    {
        return total_time;
    }

private:
    std::chrono::time_point<std::chrono::steady_clock> start_time = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_time = std::chrono::duration<double>::zero();
    double total_time{};
};

#pragma once

#include <chrono>


// -Timer overhead benchmark: ~84.9 ns per start/stop cycle

struct XB_Timer
{
public:
    using Clock = std::chrono::steady_clock;
    using Duration = std::chrono::nanoseconds;

    // Starts or restarts the timer.
    void start()
    {
        start_time = Clock::now();
        is_running = true;
    }

    // Stops the timer, calculates the last elapsed duration, and adds it to the total.
    void stop()
    {
        if (is_running)
        {
            last_elapsed_time = Clock::now() - start_time;
            total_duration += last_elapsed_time;
            is_running = false;
        }
    }

    // Resets the total accumulated time to zero.
    void reset()
    {
        total_duration = Duration::zero();
        last_elapsed_time = Duration::zero();
        is_running = false;
    }

    // --- Getters ---
    [[nodiscard]] Duration last_elapsed() const { return last_elapsed_time; }
    [[nodiscard]] Duration total() const { return total_duration; }
    [[nodiscard]] Duration current_elapsed() const
    {
        return is_running ? (Clock::now() - start_time) : Duration::zero();
    }

private:
    std::chrono::time_point<Clock> start_time;
    Duration last_elapsed_time{};
    Duration total_duration{};
    bool is_running = false;
};

// The RAII helper class
class ScopedTimer
{
public:
    // 1. Constructor: When a ScopedTimer is created, it takes a reference
    //    to an XB_Timer and immediately calls start() on it.
    explicit ScopedTimer(XB_Timer& timer) 
        : timer_ref(timer)
    {
        timer_ref.start();
    }

    // 2. Destructor: When the ScopedTimer goes out of scope (e.g., at the
    //    end of a function or block), this is automatically called, stopping the timer.
    ~ScopedTimer()
    {
        timer_ref.stop();
    }

    // 3. Disable Copying and Moving: A scoped timer's job is unique to its
    //    scope. Copying it would lead to confusion (e.g., who stops the timer?).
    //    Deleting these prevents bugs.
    ScopedTimer(const ScopedTimer&) = delete;
    ScopedTimer& operator=(const ScopedTimer&) = delete;

private:
    XB_Timer& timer_ref; // Holds a reference to the actual timer object
};

struct SimulationTime{
    XB_Timer total_time;
    XB_Timer parse_time;
    XB_Timer setup_time;        // ckt setup time
    XB_Timer analysis_time;     // dc, tran, ac analysis time
    XB_Timer matrix_load_time;  // matrix loading time (Device stamping time)
    XB_Timer solve_time;        // solver time
    XB_Timer newton_time;       // Newton's method time
    XB_Timer bsim4_time;        // BSIM4 device evaluation time
    // Profiling timers for non-Newton code
    XB_Timer get_cap_state_time;     // get_cap_state() time
    XB_Timer lte_calc_time;          // single_LTE_ngspice() time
    XB_Timer update_device_time;     // updateDeviceState() time
    XB_Timer history_update_time;    // history_trans_update() time
    XB_Timer nicomcof_time;          // NIcomCof() time
};

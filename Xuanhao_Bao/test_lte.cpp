// #include <iostream>
// #include <chrono>
// #include <thread>
// #include "BS_thread_pool.hpp"
// #include "BS_thread_pool_utils.hpp"
#include "Transient_code_parser.hpp"

XB_Timer timer;
void measureTimeWithSteadyClock() {
    timer.start();
    // Simulate some work with a sleep
    std::this_thread::sleep_for(std::chrono::milliseconds(100));
    timer.stop();
    // auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Steady clock duration: " << timer.ms() << " milliseconds\n";
}

int main() {
    measureTimeWithSteadyClock();
    return 0;
}
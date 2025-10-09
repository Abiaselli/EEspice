#include <iostream>
#include <sys/prctl.h>
#include <cstdlib>
#include <chrono>
#include <iomanip>
#include <vector>

#ifndef PR_TASK_PERF_EVENTS_DISABLE
#define PR_TASK_PERF_EVENTS_DISABLE 31
#define PR_TASK_PERF_EVENTS_ENABLE  32
#endif

int main(int argc, char* argv[]){
    // Parse command-line arguments
    size_t size = 10000;        // Default array size
    size_t iterations = 10000;  // Default number of iterations

    if (argc > 1) {
        size = std::atol(argv[1]);
    }
    if (argc > 2) {
        iterations = std::atol(argv[2]);
    }

    // Validate inputs
    if (size == 0 || iterations == 0) {
        std::cerr << "Usage: " << argv[0] << " [array_size] [iterations]" << std::endl;
        std::cerr << "  array_size: Size of array to test (default: 10000)" << std::endl;
        std::cerr << "  iterations: Number of times to repeat summation (default: 10000)" << std::endl;
        return 1;
    }

    std::cout << "Array size: " << size << ", Iterations: " << iterations << std::endl;

    double c=0;
    std::vector<double> a(size);

    // Initialize array
    for (size_t i=0; i<size; i++){
        a[i]=i*1.0;
    }
    auto start = std::chrono::high_resolution_clock::now();
    // Run the summation loop many times for measurable cache statistics
    prctl(PR_TASK_PERF_EVENTS_ENABLE, 0, 0, 0, 0);
    for (size_t iter = 0; iter < iterations; iter++) {
        for (size_t i = 0; i < size; i++) {
            c = a[i] + c;
        }
    }
    prctl(PR_TASK_PERF_EVENTS_DISABLE, 0, 0, 0, 0);
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "c = " << c << std::endl;

    std::chrono::duration<double> total_time = end - start;
    std::cout << "Total time: " << std::fixed << std::setprecision(6)
              << total_time.count() << " seconds\n";

    return 0;
}
#include <iostream>
#include <cstdlib>
#include <chrono>
#include <iomanip>
#include <vector>

void benchmark(const size_t size, const size_t iterations, double &c, std::vector<double> &a) {
    for (size_t iter = 0; iter < iterations; iter++) {
        for (size_t i = 0; i < size; i++) {
            c = a[i] + c;
        }
    }
}

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
    benchmark(size, iterations, c, a);
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "c = " << c << std::endl;

    std::chrono::duration<double, std::milli> elapsed = end - start;
    double time_ms = elapsed.count();

    std::cout << "benchmark Time taken: " << time_ms << " ms" << std::endl;

    return 0;
}
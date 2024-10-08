#define ARMA_DONT_USE_WRAPPER
#define ARMA_USE_MKL_ALLOC
// #include "Transient_code_parser.hpp"

#include <armadillo>
#include <iostream>
#include <chrono>
#include <vector>
// #include "BS_thread_pool.hpp"
// #include <mutex>


int main(){

    double a{};

    for(int i = 0; i < 5; i++){
        a += 1;
        double b = std::move(a);
        std::cout << b << std::endl;
    }

    return 0;
}


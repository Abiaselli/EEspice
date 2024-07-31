#define ARMA_DONT_USE_WRAPPER
#define ARMA_USE_MKL_ALLOC
#include <armadillo>
#include <iostream>
// #include "Transient_code_parser.hpp"

int main() {
    arma::mat A = arma::randu<arma::mat>(100,100);
    arma::mat B = arma::randu<arma::mat>(10,1);

    arma::mat X;
    arma::vec xx;
    arma::vec bb = arma::vec(B);
    arma::vec cc = bb;

    // if (arma::approx_equal(B, bb, "absdiff", 1e-15)) {
    //     std::cout << "B and bb are exactly equal." << std::endl;
    // } else {
    //     std::cout << "B and bb are not equal." << std::endl;
    // }

    // auto start = std::chrono::high_resolution_clock::now();
    // // X = arma::solve(A,B, arma::solve_opts::fast);
    // bb.row(2).col(0) = 1.0;
    // auto end = std::chrono::high_resolution_clock::now();

    // auto t1 = std::chrono::high_resolution_clock::now();
    // // xx = arma::solve(A,bb, arma::solve_opts::fast);
    // cc(2) = 1.0;
    // auto t2 = std::chrono::high_resolution_clock::now();

    // auto duration = std::chrono::duration<double, std::milli>(end - start);
    // auto duration2 = std::chrono::duration<double, std::milli>(t2 - t1);

    // std::cout << "Time bb.row(2): " << duration.count() << " microseconds" << std::endl;
    // std::cout << "Time bb(2): " << duration2.count() << " microseconds" << std::endl;

    arma::vec pre_voltages = bb.submat(0, 0, 6, 0);
    bb.print("bb:");
    pre_voltages.print("pre_voltages:");

    return 0;
}



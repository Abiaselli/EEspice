#include "Transient_code_parser.hpp"
// #include <armadillo>
// #include <iostream>



int main() {
    // arma::mat solution = arma::zeros<arma::mat>(3, 1);
    double v = V_pulse_value(0, 5, 1e-10, 0, 5e-6, 5e-6, 2e-6, 20e-6);
    std::cout << v << std::endl;


    return 0;
}




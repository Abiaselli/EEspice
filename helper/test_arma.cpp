// Build smoke test for Armadillo: compile and run this file to verify
// that the Armadillo headers and library are available.
#include <armadillo>
#include <iostream>


int main() {
    std::cout << "Armadillo version: " << arma::arma_version::as_string() << '\n';
    return 0;
}


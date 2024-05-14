#include <iostream>
#include <iomanip>
#include <armadillo>
// #include "Transient_code_parser.hpp"

int main() {
    // arma::mat A = {
    //     {1000, -1000, 1},
    //     {-1000, 1000, 0},
    //     {1, 0, 0}
    // };
    // A.print("A:");

    // arma::vec b = {{2}, 
    //                {12},
    //                {0.0001}};
    // b.print("b:");

    // arma::mat x = arma::solve(A, b);
    // x.print("x:");

    // arma::vec c = {{5}, 
    //                {3},
    //                {0.000112343243424}};

    // // Printing using a loop to maintain precision
    // std::cout << std::scientific << std::setprecision(15);
    // std::cout << "c =\n";
    // for (size_t i = 0; i < c.n_elem; ++i) {
    //     std::cout << c(i) << std::endl;
    // }

    // std::cout << std::fixed << std::setprecision(15);
    // std::cout << "c (fixed format) =\n";
    // for (size_t i = 0; i < c.n_elem; ++i) {
    //     std::cout << c(i) << std::endl;
    // }

    // double a = 10 * 1e-15;
    // std::cout << std::scientific << std::setprecision(20) << a << std::endl;
    
    arma::vec next_voltages = {
        {1},
        {1},
        {10}};
    arma::vec current_voltages = {
        {5},
        {5},
        {5}};

        // arma::mat Vmax_next_current = arma::max(arma::abs(next_voltages), arma::abs(current_voltages));
        // Vmax_next_current.print("Vmax_next_current:");
        // arma::mat Vdelta_next_current = arma::abs(next_voltages - current_voltages);
        // Vdelta_next_current.print("Vdelta_next_current:");
        // std::cout << " tolerance: " << RELTOL * Vmax_next_current + VNTOL << std::endl;


        // If |v(k+1)-v(k)| is bigger arma::any will return true and the function will return false.
        // if (arma::any(next_voltages > current_voltages))
        // {   
        //     std::cout << "False" << std::endl; 
        //     return false;
        // }
        // std::cout << "True" << std::endl;

        arma::mat lte;
        lte = arma::max(next_voltages, current_voltages);
        lte.print("lte:");

        double asd = next_voltages.min();
        std::cout << "asd: " << asd << std::endl;

    
    return 0;
}



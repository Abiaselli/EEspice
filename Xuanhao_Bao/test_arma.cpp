#include <iostream>
#include <iomanip>
#include <armadillo>
// #include "Transient_code_parser.hpp"

int main() {
    arma::mat A = arma::zeros<arma::mat>(68,68);
    arma::mat B = arma::zeros<arma::mat>(68,68);

    
    for(int i=0; i<8243479; i++){
        // A = A + B;
        A.row(3).col(6) = 3.3;
        A.row(2).col(5) = -3.3;
        A.row(4).col(11) = 3.3;
        A.row(26).col(6) = -3.3;
    }
    

    
    return 0;
}



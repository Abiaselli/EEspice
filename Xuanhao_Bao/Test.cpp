#include "Transient_code_parser.hpp"
// #include <iostream>
// #include <armadillo>


int main() {
    arma::mat A;
    A.randn(68,68);

    arma::sp_mat B(68,68);
    B(1,1) = 3.3;
    B(2,2) = -3.3;
    B(4,11) = 3.3;
    B(26,6) = -3.3;
    B.print("B: ");

    arma::mat C(68,68);
    C(1,1) = 3.3;
    C(2,2) = -3.3;
    C(4,11) = 3.3;
    C(26,6) = -3.3;

    auto t1 = std::chrono::high_resolution_clock::now(); // Start time
    for(int i=0; i<8243479; i++){
       
        A = A + B; //Sparse matrix

    }
    auto t2 = std::chrono::high_resolution_clock::now(); // End time

    auto t3 = std::chrono::high_resolution_clock::now(); // Start time
    for(int i=0; i<8243479; i++){

        A.row(1).col(1) += 3.3;
        A.row(2).col(2) += -3.3;
        A.row(4).col(11) += 3.3;
        A.row(26).col(6) += -3.3;
    }
    auto t4 = std::chrono::high_resolution_clock::now(); // End time

    auto t5 = std::chrono::high_resolution_clock::now(); // Start time
    for(int i=0; i<8243479; i++){
       
        A = A + C; //Dense matrix

    }
    auto t6 = std::chrono::high_resolution_clock::now(); // Start time



    std::chrono::duration<double, std::milli> time_span = (t2 - t1);
    std::chrono::duration<double, std::milli> time = (t4 - t3);
    std::chrono::duration<double, std::milli> time_span2 = (t6 - t5);
    std::cout << "Total SP adding time:" <<  time_span.count() << "ms\n";
    std::cout << "Total modification time:" <<  time.count() << "ms\n";
    std::cout << "Total DENSE adding time:" <<  time_span2.count() << "ms\n";
    return 0;
}

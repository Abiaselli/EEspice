#include "Transient_code_parser.hpp"
// #include <armadillo>
// #include <iostream>

// Custom test function to compare matrices and print result

int main() {
    CKTcircuit ckt;
    DenseMatrix dematrix;
    Transient trans_op;

    CircuitParser parser("Inverter.cir");
    parser.parser();

    CKTsetup(ckt, parser, dematrix);            // Pass the parser to the ckt and the initialise LHS and RHS matrices
    ckt.setcktmatrix(dematrix);

    CKTload(ckt);
  
    return 0;
}




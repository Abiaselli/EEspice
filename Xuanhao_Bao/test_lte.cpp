#include "Transient_code_parser.hpp"

int main() {

arma::mat A = arma::randu(2,1);
arma::mat B = arma::randu(2,1);
A.print("A:");
B.print("B:");

arma::mat C = A % B;
C.print("C:");



return 0;
}
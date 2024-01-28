#include <iostream>
#include <armadillo>

// Voltage pulse assigner
double V_pulse_new(double V1, double V2, double t1, double td, double tr, double tf, double tpw, double tper)
{   
    double tnorm = fmod((t1-td),tper);
    double v = 0;
    const double EPSILON = 1e-18;


    if(t1 < 0){
        std::cout << "Simulation time in pulse voltage it wrong" << std::endl;
        exit(1);
    }
    else if(0 < t1 < td) {
        v = V1;
    }
    else if(0 < tnorm && tnorm < tr){
        v = (V2 - V1)/tr * tnorm;
    }
    else if(tr < tnorm && tnorm < tr + tpw){
        v = V2;
    }
    else if(tr + tpw < tnorm && tnorm < tper){
        v = V2 + (V1-V2) / tf * (tnorm - (tpw + tr));
    }
    else {
        std::cout << "Pulse voltage Error" << std::endl;
        exit(1);
    }
    
    return v;
}



int main(){
    double v = V_pulse_new(0, 5, 2.4e-7, 0, 1e-5, 1e-5, 1e-3, 1.02e-3);
    std::cout << v << std::endl;

    return 0;
}






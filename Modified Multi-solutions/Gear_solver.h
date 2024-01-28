#include <iostream>
#include <vector>
#include <cmath>


// double RELTOL = 10e-3;  // Relative voltage tolerance
// double ABSTOL = 10e-12; // ABSTOL is the absolute current tolerance
// double CHGTOL = 10e-14; // Introduced for the absolute charge or flux error
// double TRTOL = 7;   // Scale down the divided LTE(Local truncation error)


class GearSolver {
private:
    double h;                   // Current Timestep
    // const double h_min, h_max;  // Minimum and maximum allowable timesteps
    double RELTOL = 10e-3;      // Relative voltage tolerance
    double ABSTOL = 10e-12;     // ABSTOL is the absolute current tolerance
    double CHGTOL = 10e-14;     // Introduced for the absolute charge or flux error
    double TRTOL = 7;           // Scale down the divided LTE(Local truncation error)
    double eps_val;             // Absolute Error
    double error_val;           // Relative tolerance error
    double x_n[4];              // voltage of capacitor; x(t_n), x(t_(n+1)), x(t_(n+2), x(t_(n+3)) 
    double i_n[4];              // current of capacitor
    double h_history[2];        // h_(n-1), h_(n-2)
    // std::deque<arma::vec> voltages_history;
    // std::deque<double> h_history;
    int n_op;                   // Indicate the number of operation point

    double truncation_error(){
        return RELTOL * std::max(std::abs(x_n[1]), std::abs(x_n[0])) + ABSTOL;
    }

    double flux_error(){
        return RELTOL * std::max({std::abs(x_n[1]), std::abs(x_n[0]), CHGTOL}) / h;
    }

    // LTE
    double E(){
        return std::max(truncation_error(), flux_error());
    }

    // First divided difference
    double DD_1(int i){
        if (i + 1 >= 3) return 0;
        return (x_n[1 + i] - x_n[i]) / h;
    }

    // Second divided difference
    double DD_2(int i){
        if (i + 1 >= 3) return 0;
        return (DD_1(1 + i) - DD_1(0 + i)) / (h + h_history[0]);
    }

    // Third divided difference
    double DD_3(){
        return (DD_2(1) - DD_2(0)) / (h + h_history[0] + h_history[1]);
    }

    double h_next(int n_op){
        // Fot the first operation point
        if (n_op == 0){
            return sqrt(TRTOL * E() / std::max(DD_1(0) / 12 , eps_val));
        }
        // For the second operation point
        else if (n_op == 1){
            return sqrt(TRTOL * E() / std::max(DD_2(0) / 12 , eps_val));
        }
        else{
            return sqrt(TRTOL * E() / std::max(DD_3() / 12 , eps_val));
        }
    }

public:
    GearSolver(double timestep, double eps_val, double error_val, double input_x_n[4], double input_h_history[2], int n_op)
        : h(timestep), eps_val(eps_val), error_val(error_val) {
            for (int i = 0; i < 4; ++i) {
                x_n[i] = input_x_n[i];
            }

            for (int i = 0; i < 2; ++i) {
                h_history[i] = input_h_history[i];
            }
    }

    double getNextTimeStep() {
        return h_next(n_op);
    }
};

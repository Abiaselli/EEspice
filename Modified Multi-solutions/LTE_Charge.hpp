#pragma once
#include <iostream>
#include <fstream>
#include <armadillo>
#include <tuple>
#include <deque>
#include <vector>
#include <cmath>
#include <chrono>

// Convergence settings
double MINR = 1e-3; // 1 milli ohm

// Transient analysis settings
/* This setting is to turn on or off the adaptive steps*/
int ENADAPT = 1;

/* Mode for timestep, 1: LTE; 2:Iteration Count; */
double STEPMODE = 1;

/* Scale factor for predicting*/
double TRTOL = 7;

/* Error Charge Factor*/
double CHGTOL = 1e-14;

/* Relative error*/
double RELTOL = 1e-3;

/* Absolute error for voltage*/
double VNTOL = 1e-6;
/* Absolute error for current*/
double ABSTOL = 1e-12;

// Relative charge tolerance for truncation error
double relq = 0.01;
// Relative current tolerance for truncation error
double lteretol = 0.01;
// Absolute tolerance for truncation error
double lteabstol = 1e-6;

/* Min and Max timestep */
double TMIN = 1e-9 * 1e-5;  // rmin * tstep
double TMAX = 5 * 1e-4;     // rmax * step
// double TMAX = t_end/50;
/* Threshold for iteration count-based algorithm */
int ITL3 = 4;
int ITL4 = 10;
bool ITL4_reach = false;

/*If timestep reach the limit*/
bool TMAX_reach = false;
bool TMIN_reach = false;

// History timestep list
std::vector<double> h_history;      // Practical timestep
std::vector<double> h_1_history;    // Expected timestep

struct x_mid_up_down;
class capacitor;
class Truncation_error;
std::vector<capacitor> C_list_up;
std::vector<capacitor> C_list_mid;
std::vector<capacitor> C_list_down;



// Matrix stamps assigner using Modified Nodal Analysis
std::pair<arma::mat, arma::mat> DynamicNonLinear(arma::mat &LHS, arma::mat &RHS,
                                                                       arma::mat solution, std::deque<arma::vec> &history_voltages,
                                                                       double h, int mode);
void R_assigner(double node_x, double node_y, double R, arma::mat &LHS, arma::mat &RHS);
void Is_assigner(double node_x, double node_y, double I, arma::mat &LHS, arma::mat &RHS);
double Vs_assigner(int node_x, int node_y, double V_value, arma::mat &LHS, arma::mat &RHS);
double V_pulse(double V1, double V2, double t1, double td, double tr, double tf, double tpw, double tper);
void C_assigner(int node_x, int node_y, double C, double h, arma::mat &LHS, arma::mat &RHS, arma::mat solution, int mode);
void Diode_assigner(int node_x, int node_y, double Is, double VT, double cd, double h, arma::mat &LHS, arma::mat &RHS, arma::mat solution, int mode);
void VCCS_assigner(int node_x, int node_y, int node_cx, int node_cy, double R, arma::mat &LHS);
void C_assigner_3(int node_x, int node_y, double C, double h, arma::mat &LHS, arma::mat &RHS,
                  std::deque<arma::mat> &NR_solutions,
                  int mode);
double Diode_assigner_2(int node_x, int node_y, double Is, double VT, double h, arma::mat &LHS, arma::mat &RHS, 
                      arma::mat solution, int mode);
bool isConverge(const std::deque<arma::mat> &NR_solutions);

void UpdateStates(arma::mat &LHS, arma::mat &RHS, 
                    std::deque<arma::mat> &NR_solutions, 
                    double h, int mode);

void history_voltages_update(arma::mat &solution, std::deque<arma::vec> &history_voltages);
arma::mat NewtonRaphson_system(arma::mat &LHS, arma::mat &RHS, arma::mat solution, std::deque<arma::vec> &history_voltages,
                               double h, int mode);
x_mid_up_down  multi_solution_solver(arma::mat &temp_LHS_mid, arma::mat &temp_LHS_up, arma::mat &temp_LHS_down, arma::mat const init_RHS, arma::mat LHS, arma::mat solution, std::deque<arma::vec> &history_voltages, 
                            arma::mat &solution_mid, arma::mat &solution_up, arma::mat &solution_down, std::vector<double> RHS_locate, arma::mat &temp_RHS_mid, 
                            arma::mat &temp_RHS_up, arma::mat &temp_RHS_down, double & temp_h, double V1, double V2, double & t1_pulse, double td, double tr, 
                            double tf, double tpw, double tper, int mode, double time_trans, x_mid_up_down mid_up_down);
int h_reach_new(Truncation_error LTE_h, int &last_decision);


// Sum the matrices inside the vector
arma::mat mat_sum(std::vector<arma::mat> vector_of_matrices)
{
    arma::mat a = vector_of_matrices[0];
    for (long unsigned int i = 1; i < vector_of_matrices.size(); ++i)
    {
        a = a + vector_of_matrices[i];
    }
    return a;
}

// Testing LU solve using the triangular pivoting system
arma::mat LU_solve(arma::mat A, arma::mat b)
{
    arma::mat L, U, P;
    arma::lu(L, U, P, A);
    arma::vec x1 = arma::solve(trimatu(U), solve(trimatl(L), P * b));

    return x1;
}

// Changing resistance to conductance
double cond(double R)
{
    return 1 / R;
}

// extra -  adding arange() function here to return the time vector
arma::vec arange(double tstart, std::deque<double> history_steps)
{
    arma::vec time = arma::zeros(history_steps.size(), 1);
    time[0] = tstart;
    for (size_t i = 1; i < history_steps.size(); i++)
    {
        tstart = tstart + history_steps[history_steps.size() - i];
        time[i] = tstart;
    }
    return time;
}

// For testing purposes
arma::vec arange_V(double tstart, double h, double vec_size){
    arma::vec time = arma::zeros(vec_size,1);
    for(int i=0; i<vec_size; i++)
    {
        time[i] = tstart;
        tstart = tstart + h;
    }
    return time;
}

// assigning the matrix stamps for the VCCS
void VCCS_assigner(int node_x, int node_y, int node_cx, int node_cy, double R, arma::mat &LHS)
{
    int maxi = LHS.n_cols;
    int maxj = LHS.n_rows;
    arma::mat a = arma::zeros(maxi, maxj);

    if (node_x == 0)
    {
        if (node_cx == 0)
        {
            if (node_cy > 0)
            {
                a.row(node_y - 1).col(node_cy - 1) = R;
            }
        }
        else if (node_cy == 0)
        {
            if (node_cx > 0)
            {
                a.row(node_y - 1).col(node_cx - 1) = -R;
            }
        }
        else
        {
            a.row(node_y - 1).col(node_cx - 1) = -R;
            a.row(node_y - 1).col(node_cy - 1) = R;
        }
    }
    else if (node_y == 0)
    {
        if (node_cx == 0)
        {
            if (node_cy > 0)
            {
                a.row(node_x - 1).col(node_cy - 1) = -R;
            }
        }
        else if (node_cy == 0)
        {
            if (node_cx > 0)
            {
                a.row(node_x - 1).col(node_cx - 1) = R;
            }
        }
        else
        {
            a.row(node_x - 1).col(node_cx - 1) = R;
            a.row(node_x - 1).col(node_cy - 1) = -R;
        }
    }
    else
    {
        if (node_cx == 0)
        {
            if (node_cy > 0)
            {
                a.row(node_x - 1).col(node_cy - 1) = -R;
                a.row(node_y - 1).col(node_cy - 1) = R;
            }
        }
        else if (node_cy == 0)
        {
            if (node_cx > 0)
            {
                a.row(node_x - 1).col(node_cx - 1) = R;
                a.row(node_y - 1).col(node_cx - 1) = -R;
            }
        }
        else
        {
            a.row(node_x - 1).col(node_cx - 1) = R;
            a.row(node_x - 1).col(node_cy - 1) = -R;
            a.row(node_y - 1).col(node_cx - 1) = -R;
            a.row(node_y - 1).col(node_cy - 1) = R;
        }
    }
    LHS = LHS + a;
}

// create resistor matrix stamp
void R_assigner(double node_x, double node_y, double R, arma::mat &LHS, arma::mat &RHS)
{
    int maxi = LHS.n_cols;
    int maxj = LHS.n_rows;
    double x = 0;
    arma::mat a = arma::zeros(maxi, maxj);

    // // Compare with MINR
    // if (R <= MINR)
    //     R = MINR;

    if (R == 0)
        x = 0;
    else
        x = cond(R);

    if ((node_x == 0) && (node_y == 0))
    {
        x = 0;
    }
    else
    {
        if (node_x == 0)
        {
            a.row(node_y - 1).col(node_y - 1) = x;
        }
        else if (node_y == 0)
        {
            a.row(node_x - 1).col(node_x - 1) = x;
        }
        else
        {
            a.row(node_x - 1).col(node_x - 1) = x;
            a.row(node_x - 1).col(node_y - 1) = -x;
            a.row(node_y - 1).col(node_x - 1) = -x;
            a.row(node_y - 1).col(node_y - 1) = x;
        }
    }
    LHS = LHS + a;
}

// create a current matrix stamp
void Is_assigner(double node_x, double node_y, double I, arma::mat &LHS, arma::mat &RHS)
{
    int maxi = LHS.n_cols;
    int maxj = 1;
    arma::mat a = arma::zeros(maxi, maxj);
    if ((node_x == 0) && (node_y == 0))
    {
        I = 0;
    }
    else
    {
        if (node_x == 0)
        {
            a.row(node_y - 1).col(0) = I;
        }
        else if (node_y == 0)
        {
            a.row(node_x - 1).col(0) = -I;
        }
        else
        {
            a.row(node_x - 1).col(0) = -I;
            a.row(node_y - 1).col(0) = I;
        }
    }
    RHS = RHS + a;
}
// branch extender function
arma::mat branch_ext(arma::mat M, int node_x, int node_y)
{
    int M_cols = M.n_cols;
    int M_rows = M.n_rows;
    arma::mat va = arma::zeros(M_rows, 1);
    if (node_x == 0)
    {
        va.row(node_y - 1).col(0) = 1;
    }
    else if (node_y == 0)
    {
        va.row(node_x - 1).col(0) = 1;
    }
    else
    {
        va.row(node_x - 1).col(0) = -1;
        va.row(node_y - 1).col(0) = 1;
    }

    arma::mat zero_ext = arma::zeros(1, 1);
    arma::mat ha = va.as_row();
    arma::mat haz = arma::join_rows(ha, zero_ext);

    arma::mat M1 = arma::join_rows(M, va);
    arma::mat M2 = arma::join_cols(M1, haz);

    return M2;
}



// Voltage source stamp assigner
double Vs_assigner(int node_x, int node_y, double V_value, arma::mat &LHS, arma::mat &RHS)
{

    arma::vec value(1);
    value = V_value;
    // Extending the branch at the LHS matrix
    LHS = branch_ext(LHS, node_x, node_y);

    // Assigning the value at RHS
    RHS = arma::join_cols(RHS, value);

    // location of voltage value for transient simulation
    double V_locate1 = RHS.n_rows;

    return V_locate1;
}

// Capacitor stamp assigner using FE(Forward Eulur) method
void C_assigner(int node_x, int node_y, double C, double h, arma::mat &LHS, arma::mat &RHS, arma::mat solution, int mode)
{

    double x = 0;
    double x1 = 0;

    // this if else statment uses trapezoidal formula
    if (mode == 1)
    {
        x = C / h;
        if (node_x == 0)
        {
            x1 = C * (solution(node_y - 1, 0)) / h;
        }
        else if (node_y == 0)
        {
            x1 = C * (solution(node_x - 1, 0)) / h;
        }
        else
        {
            x1 = C * (solution(node_x - 1, 0) - solution(node_y - 1, 0)) / h;
        }
    }
    else
    {
        x = 0;
        x1 = 0;
    }
    // Matrix stamp for a capacitor on LHS
    R_assigner(node_x, node_y, cond(x), LHS, RHS);
    // Matrix stamp for a capacitor on RHS
    Is_assigner(node_x, node_y, -x1, LHS, RHS);
}




// Capacitor stamp assigner using BE(Backward Eulur) method
void C_assigner_3(int node_x, int node_y, double C, double h, arma::mat &LHS, arma::mat &RHS,
                  std::deque<arma::mat> &NR_solutions,
                  int mode)
{
    // Ensure the vector is large enough to contain the capacitor at the given index
    // if (index >= capacitorStates.size()) {
    //     capacitorStates.resize(index + 1);
    // }

    double x = 0;
    double x1 = 0;
    double vol = 0;

    // CapacitorState &state = capacitorStates[index];
    // this if else statment uses trapezoidal formula
    if (mode == 1)
    {
        // v_{n-1} 
        arma::mat pre_solution = NR_solutions.back();
        pre_solution.print("pre_solution in c is ");
        
        x = C / h; // x is the equivalent resistance, x1 is the equivalent current source
        x1 = C / h;

        if (node_x == 0)
        {
            vol = pre_solution(node_y - 1, 0);
        }
        else if (node_y == 0)
        {
            vol = pre_solution(node_x - 1, 0);
        }
        else
        {
            vol = pre_solution(node_x - 1, 0) - pre_solution(node_y - 1, 0);
        }

        // print current and pre_solution
        // std::cout << "previous_current: " << previous_current << std::endl;
        // std::cout << "pre_solution: " << vol << std::endl;
    }
    else
    {
        x = 0;
        x1 = 0;
    }

    // Matrix stamp for a capacitor on LHS
    // R_assigner(node_x, node_y, cond(x), LHS, RHS);
    R_assigner(node_x, node_y, cond(x), LHS, RHS);
    // Matrix stamp for a capacitor on RHS
    Is_assigner(node_x, node_y, (-(x1 * vol)), LHS, RHS);
    
}

class capacitor{
    public:
        int id{};
        int node_x{};
        int node_y{};
        double C{};
        double current{};
        double pre_current{};
        double charge{};
        double pre_charge{};

        std::vector<double> temp_current_up;
        std::vector<double> temp_charge_up;
};
// Voltage pulse assigner
double V_pulse(double V1, double V2, double t1, double td, double tr, double tf, double tpw, double tper){
    double v{};
    t1 = fmod(t1, tper);
    if(tper < tr + tpw + tf){
        std::cout << "Period is incorrect" << std::endl;
        exit(1);
    }

    if(t1 < 0){
        std::cout << "Simulation time in pulse voltage it wrong" << std::endl;
        exit(1);
    }
    else if(t1 >= 0 && t1 < td) {
        v = V1;
    }
    else if(t1 >= td && t1 < td+tr){
        v = V1 + (V2-V1)*(t1-td)/tr;
    }
    else if(t1 >= td+tr && t1 < td+tr+tpw){
        v = V2;
    }
    else if(t1 >= td+tr+tpw && t1 < td+tr+tpw+tf){
        v = V2 + (V1-V2)*(t1-td-tr-tpw)/tf;
    }
    else if(t1 >= td+tr+tpw+tf && t1 < tper){
        v = V1;
    }
    else {
        std::cout << "Pulse voltage Error" << std::endl;
        exit(1);
    }
    return v;
}

// Voltage pulse assigner
double V_pulse_new(double V1, double V2, double t1, double td, double tr, double tf, double tpw, double tper)
{   
    double tnorm = fmod((t1-td),tper);
    double v = 0;
    const double EPSILON = 1e-18;
    if(tper < tr + tpw + tf){
        std::cout << "Period is incorrect" << std::endl;
        exit(1);
    }

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

double Vstep_Source(double &V_start, double V_step, double h)
{
    V_start = V_start + h;
    return V_start;
}

double Vsin_Source(double V, double omega, double time)
{
    // V = v * sin(wt)
    V = V * sin(omega * time);
    return V;
}

// Diode stamp assigner
void Diode_assigner(int node_x, int node_y, double Is, double VT, double cd, double h, arma::mat &LHS, arma::mat &RHS, arma::mat solution, int mode)
{
    int col_size = LHS.n_cols;
    int row_size = LHS.n_rows;

    double x = 0;
    double x1 = 0;
    double x2 = 0;
    double val_nodex = 0;
    double val_nodey = 0;

    if (mode == 0)
    {
        cd = 0;
    }

    if ((node_x == 0 && node_y == 0))
    {
        x = 0;
        x1 = 0;
    }
    else if (node_x == 0)
    {
        x = (Is / VT) * (exp(-val_nodey / VT));
        x1 = (x * (-val_nodey) - Is * (exp((-val_nodey) / VT) - 1));
    }
    else if (node_y == 0)
    {
        x = (Is / VT) * (exp((val_nodex) / VT));
        x1 = (x * (val_nodex)-Is * (exp((val_nodex) / VT) - 1));
    }
    else
    {
        x = (Is / VT) * (exp((val_nodex - val_nodey) / VT));
        x1 = (x * (val_nodex - val_nodey) - Is * (exp((val_nodex - val_nodey) / VT) - 1));
    }

    // Matrix stamp for a diode on RHS
    Is_assigner(node_x, node_y, x1, LHS, RHS);
    // Matrix stamp for a diode on LHS
    R_assigner(node_x, node_y, cond(x), LHS, RHS);
    C_assigner(node_x, node_y, cd, h, LHS, RHS, solution, mode);
}

// Correct Diode stamp assigner (Original Diode_assigner is wrong with the node voltage assignment)
double Diode_assigner_2(int node_x, int node_y, double Is, double VT, double h, arma::mat &LHS, arma::mat &RHS, 
                      arma::mat solution, int mode)
{
    int maxi = LHS.n_cols;
    int maxj = 1;

    double x = 0;
    double x1 = 0;
    double val_nodex = 0;
    double val_nodey = 0;

    double Id = 0;
    double vol = 0;

    if ((node_x == 0 && node_y == 0))
    {
        x = 0;
        x1 = 0;
        vol = 0;
    }
    else if (node_x == 0)
    {
        val_nodey = solution(node_y - 1, 0);
        vol = -val_nodey;
        x = (Is / VT) * (exp(-val_nodey / VT));
        x1 = (x * (-val_nodey) - Is * (exp((-val_nodey) / VT) - 1));
    }
    else if (node_y == 0)
    {
        val_nodex = solution(node_x - 1, 0);
        vol = val_nodex;
        x = (Is / VT) * (exp((val_nodex) / VT));
        x1 = (x * (val_nodex)-Is * (exp((val_nodex) / VT) - 1));
    }
    else
    {
        val_nodex = solution(node_x - 1, 0);
        val_nodey = solution(node_y - 1, 0);
        vol = val_nodex - val_nodey;
        x = (Is / VT) * (exp((val_nodex - val_nodey) / VT));
        x1 = (x * (val_nodex - val_nodey) - Is * (exp((val_nodex - val_nodey) / VT) - 1));
    }

    // Matrix stamp for a diode on RHS
    Is_assigner(node_x, node_y, x1, LHS, RHS);
    // Matrix stamp for a diode on LHS
    R_assigner(node_x, node_y, cond(x), LHS, RHS);

    Id = Is * (exp(vol / VT) - 1);
    return Id;
}

// Newton Raphson system solver for non-linear and dynamic elements
arma::mat NewtonRaphson_system(arma::mat &LHS, arma::mat &RHS, arma::mat solution, std::deque<arma::vec> &history_voltages,
                               double h, int mode)
{
    std::cout << "Enter NewtonRaphson system" << std::endl;

    int NR_iteration_counter = 0;
    bool isconverge = false;

    std::deque<arma::mat> NR_solutions;
    NR_solutions.push_back(history_voltages.at(0));  //input the solution from the OP analysis

    // make sure the solution is last solution.
    solution = NR_solutions.back();

    arma::mat NR_RHS = RHS;
    arma::mat NR_LHS = LHS;
    arma::mat delta; // x(k+1) - x(k)

    // Update the model of time-related elements (Capacitor)
    // Note that pre_solution is just the temporary solution or voltage during the NR iteration
    UpdateStates(NR_LHS, NR_RHS, NR_solutions, h, mode);
    NR_LHS.print("NR_LHS is ");
    NR_RHS.print("NR_RHS is ");

    
    do
    {  
        // Update the non-linear elements (Diode, MOSFET) with updated solution
        auto matrices = DynamicNonLinear(NR_LHS, NR_RHS, solution, history_voltages, h, mode);

        // f(x) = Ax - b = 0
        arma::mat f = (matrices.first * solution) - matrices.second;
        // A*[v(k+1)-v(k)] = -f(x)
        delta = arma::solve(matrices.first, -f); 
        // delta.print("delta is ");   
        solution = solution + delta;
        
        NR_iteration_counter += 1;

        //update the NR_solutions(current_solution, pre_solution, next_solution) 
        NR_solutions.push_back(solution);
        if (NR_solutions.size() > 3)
        {
            NR_solutions.pop_front();
        }

        // Check for convergence only if you have 3 solutions
        if (NR_solutions.size() == 3) {
            isconverge = isConverge(NR_solutions);
            if(isconverge){
                LHS = matrices.first;
                RHS = matrices.second;
            }
        }

        if (NR_iteration_counter > 10 && mode == 1)
        {   
            std::cout << "Not Converge at 10 iterations" << std::endl;
            ITL4_reach = true;
            break;
        }


        if (NR_iteration_counter > 100)
        {
            std::cout << "Not Converge at 100 iterations" << std::endl;
            break;
        }

    } while (!isconverge);

    // NR_LHS.print("NR_LHS is ");
    // NR_RHS.print("NR_RHS is ");

    // print iteration count
    // std::cout << "The NR Iteration count is: " << NR_iteration_counter << std::endl;

    return solution;
}

void history_voltages_update(arma::mat &solution, std::deque<arma::vec> &history_voltages)
{
    // Update voltage history
    history_voltages.push_front(solution);
    
    if (history_voltages.size() >= 6)
    {
        history_voltages.pop_back();
        history_voltages.pop_back();
        history_voltages.pop_back();
    }
}

// Updates the value of RHS
arma::mat RHS_update(std::vector<double> RHS_locate, arma::mat RHS, std::vector<double> &val)
{
    for (int i = 0; i < RHS_locate.size(); i++)
    {
        RHS.row(RHS_locate[i] - 1).col(0) = val[i];
    }

    return RHS;
}

// Checking the convergence of the NR method
bool isConverge(const std::deque<arma::mat> &NR_solutions)
{   
    // Validate inputs
    if (NR_solutions.size() < 3) {
        throw std::runtime_error("Not enough solutions in isConverge function for convergence check.");
    }

    arma::mat pre_solution = NR_solutions.front();
    arma::mat current_solution = NR_solutions.at(1);
    arma::mat next_solution = NR_solutions.back();
    // pre_solution.print("pre_solution is ");
    // current_solution.print("current_solution is ");
    // next_solution.print("next_solution is ");
    if(pre_solution.n_rows != current_solution.n_rows || pre_solution.n_rows != next_solution.n_rows || current_solution.n_rows != next_solution.n_rows){
        std::cout << "The size of pre_solution, current_solution, next_solution is not the same in isConverge function." << std::endl;
        exit(1);
    }

    arma::mat pre_voltages = pre_solution.submat(0, 0, T_nodes-1, 0);
    arma::mat current_voltages = current_solution.submat(0, 0, T_nodes-1, 0);
    arma::mat next_voltages = next_solution.submat(0, 0, T_nodes-1, 0);

    arma::mat pre_current = pre_solution.submat(T_nodes, 0, pre_solution.n_rows - 1, 0);
    arma::mat current_current = current_solution.submat(T_nodes, 0, current_solution.n_rows - 1, 0);
    arma::mat next_current = next_solution.submat(T_nodes, 0, next_solution.n_rows - 1, 0);


    // |v(k+1)-v(k)| <= RELTOL*max(|v(k+1)|,|v(k)|) + VNTOL

        arma::mat Vmax_next_current = arma::max(arma::abs(next_voltages), arma::abs(current_voltages));
        arma::mat Vdelta_next_current = arma::abs(next_voltages - current_voltages);


        // If |v(k+1)-v(k)| is bigger arma::any will return true and the function will return false.
        if (arma::any(arma::vectorise(Vdelta_next_current > RELTOL * Vmax_next_current + VNTOL)))
        {
            return false;
        }
    
 
    // |i(k+1)-i(k)| <= RELTOL*max(|i(k+1)|,|i(k)|) + ABSTOL

        arma::mat Imax_next_current = arma::max(arma::abs(next_current), arma::abs(current_current));
        arma::mat Idelta_next_current = arma::abs(next_current - current_current);

        if (arma::any(arma::vectorise(Idelta_next_current > RELTOL * Imax_next_current + ABSTOL)))
        {
            return false;
        }
    

    // |v(k) - v(k-1)| <= RELTOL*max(|v(k)|,|v(k-1)|) + VNTOL

        arma::mat Vmax_current_pre = arma::max(arma::abs(current_voltages), arma::abs(pre_voltages));
        arma::mat Vdelta_current_pre = arma::abs(current_voltages - pre_voltages);

        if (arma::any(arma::vectorise(Vdelta_current_pre > RELTOL * Vmax_current_pre + VNTOL)))
        {
            return false;
        }
    

    // |i(k) - i(k-1)| <= RELTOL*max(|i(k)|,|i(k-1)|) + ABSTOL

        arma::mat Imax_current_pre = arma::max(arma::abs(current_current), arma::abs(pre_current));
        arma::mat Idelta_current_pre = arma::abs(current_current - pre_current);

        if (arma::any(arma::vectorise(Idelta_current_pre > RELTOL * Imax_current_pre + ABSTOL)))
        {
            return false;
        }
    

    // |v(k+1) - v(k-1)| <= √|v(k) - v(k-1)|^2 + |v(k+1) - v(k)|^2
    
        arma::mat Vdelta_next_pre = arma::abs(next_voltages - pre_voltages);
        arma::mat V_Perpendicular = arma::sqrt(arma::pow(arma::abs(Vdelta_current_pre),2) + arma::pow(arma::abs(Vdelta_next_current),2));
        
        if (arma::any(arma::vectorise(V_Perpendicular)) )
        {
            return false;
        }
    

    // |i(k+1) - i(k-1)| <= √|i(k) - i(k-1)|^2 + |i(k+1) - i(k)|^2

        arma::mat Idelta_next_pre = arma::abs(next_current - pre_current);
        arma::mat I_Perpendicular = arma::sqrt(arma::pow(arma::abs(Idelta_current_pre),2) + arma::pow(arma::abs(Idelta_next_current),2));

        if (arma::any(arma::vectorise(I_Perpendicular)))
        {
            return false;
        }
    

    return true;
}

class Truncation_error{
    public:
        double LTE_current_mid{};
        double LTE_current_up{};
        double LTE_current_down{};

        double LTE_charge_mid{};
        double LTE_charge_up{};
        double LTE_charge_down{};

        double LTE_bound_mid{};
        double LTE_bound_up{};
        double LTE_bound_down{};

        double LTE_BE_mid{};
        double LTE_BE_up{};
        double LTE_BE_down{};
};

// Adatpive timestep adjustment
/* STEPMODE = 1 is the count-based adjustment: still some buggy performance in this kind of adjustment */
/* STEPMODE = 2 is the LTE-based adjustment for BE method */

struct x_mid_up_down{
    double v_mid{};  // 
    double v_up{};
    double v_down{};
    double pre_v_up{};
    double v_mid_max{};
    double v_up_max{};
    double v_down_max{};
     
    double i_mid{};
    double i_up{};
    double i_down{};
    double pre_i_up{};
    double i_mid_max{};
    double i_up_max{};
    double i_down_max{};

    double t_mid{};
    double t_up{};
    double t_down{};  
    double pre_t_up{}; 

    arma::mat RHS_mid{};
    arma::mat RHS_up{};
    arma::mat RHS_down{}; 
    arma::mat pre_RHS_up{};

    arma::mat LHS_mid{};
    arma::mat LHS_up{};
    arma::mat LHS_down{};
    arma::mat pre_LHS_up{};

    arma::mat solution_mid{};
    arma::mat solution_up{};
    arma::mat solution_down{};
    arma::mat pre_solution_up{};


};


/* 
    0 = temp_h / 2
    1 = down
    2 = mid
    3 = up
    4 = temp_h * 2
    5 = previous up
*/ 
// last_decision = 1 means all true for previous time step
// last_decision = 2 means all false for previous time step
int h_reach_new(Truncation_error LTE_h, int &last_decision){
    int reach = 0;
    bool mid = LTE_h.LTE_BE_mid <= TRTOL * LTE_h.LTE_bound_mid;
    bool up = LTE_h.LTE_BE_up <= TRTOL * LTE_h.LTE_bound_up;
    bool down = LTE_h.LTE_BE_down <= TRTOL * LTE_h.LTE_bound_down;

    if(mid == true && up == true && down == true){
        if(last_decision == 2){ // if last time step is all false
            last_decision = 1; // All true
            reach = 3;  // up
            return reach;
        }
        else{
            last_decision = 1; // All true
            if(TMAX_reach){
                reach = 3;  // up
                return reach;
            }
            reach = 4;  // temp_h * 2
            return reach;
        }
    }

    if(mid == false && up == false && down == false){
        if(last_decision == 1){ // if last time step is all true
            last_decision = 2;
            reach = 5;  // previous up
            return reach;
        }
        else{
            last_decision = 2;
            reach = 0;  // temp_h / 2
            return reach;
        }
    }

    if(mid == true && up == false && down == true){
        last_decision = 0;
        reach = 2;  // mid
        return reach;
    }

    if(mid == false && up == false && down == true){
        last_decision = 0;
        reach = 1;  // down
        return reach;
    }

    std::cout << "h_reach_new function error" << std::endl;
    exit(1);

}


// New timestep options
void timestep_options(double & temp_h, double & next_h_up, double & next_h_down) {
    if(temp_h > TMAX) {
        temp_h = TMAX;
        TMAX_reach = true;
        
    }
    if(temp_h < TMIN) {
        temp_h = TMIN;
        
    }

    next_h_up = temp_h * 1.2;
    next_h_down = temp_h / 1.2;

    if(next_h_up > TMAX) {
        next_h_up = TMAX;
        TMAX_reach = true;
        
    }
    if(next_h_up < TMIN) {
        next_h_down = TMIN;
        
    }
    if(next_h_down > TMAX) {
        next_h_up = TMAX;
        TMAX_reach = true;
        
    }
    if(next_h_down < TMIN) {
        next_h_down = TMIN;
        
    }
}

x_mid_up_down  multi_solution_solver(arma::mat &temp_LHS_mid, arma::mat &temp_LHS_up, arma::mat &temp_LHS_down, arma::mat const init_RHS, arma::mat LHS, arma::mat solution, std::deque<arma::vec> &history_voltages, 
                            arma::mat &solution_mid, arma::mat &solution_up, arma::mat &solution_down, std::vector<double> RHS_locate, arma::mat &temp_RHS_mid, 
                            arma::mat &temp_RHS_up, arma::mat &temp_RHS_down, double & temp_h, double V1, double V2, double & t1_pulse, double td, double tr, 
                            double tf, double tpw, double tper, int mode, double time_trans, x_mid_up_down mid_up_down) {
    std::cout << "multi_solution_solver" << std::endl;
    double next_h_up;
    double next_h_down;

    
    timestep_options(temp_h, next_h_up, next_h_down);
    x_mid_up_down x = mid_up_down;

    // std::cout << "temp_h:" << temp_h ;
    // std::cout << " next_h_up:" << next_h_up;
    // std::cout << " next_h_down:" << next_h_down << std::endl;
    //  exit(0);
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    double t1_mid = time_trans + temp_h;
    double t1_up = time_trans + next_h_up;
    double t1_down = time_trans + next_h_down;
    std::cout << "t1_mid" << t1_mid;
    std::cout << " t1_up" << t1_up;
    std::cout << " t1_down" << t1_down << std::endl;
    

    std::cout << "The time steps are:" << next_h_up << " " << temp_h << " " << next_h_down << std::endl;

    std::vector<double> RHS_value = {
        V_pulse(V1, V2, t1_mid, td, tr, tf, tpw, tper)
    };
    // std::cout << "RHS_value_mid" << RHS_value[0] << std::endl;

    std::vector<double> RHS_value_up = {
        V_pulse(V1, V2, t1_up, td, tr, tf, tpw, tper)
    };
    // std::cout << "RHS_value_up" << RHS_value_up[0] << std::endl;
    std::vector<double> RHS_value_down = {
        V_pulse(V1, V2, t1_down, td, tr, tf, tpw, tper)
    };
    // std::cout << "RHS_value_down" << RHS_value_down[0] << std::endl;

    temp_RHS_mid = RHS_update(RHS_locate, init_RHS, RHS_value);

    temp_RHS_up = RHS_update(RHS_locate, init_RHS, RHS_value_up);

    temp_RHS_down = RHS_update(RHS_locate, init_RHS, RHS_value_down);

    // history_voltages_update(solution, history_voltages);


    // temp_RHS_mid.print("RHS_mid matrix =");
    solution_mid = NewtonRaphson_system(temp_LHS_mid, temp_RHS_mid, solution, history_voltages, temp_h, mode); // Soluition with current time step
    // I = C/h*v(n+1) - C/h*v(n)
    C_list_mid.at(0).current = ((C_list_mid.at(0).C / temp_h) * solution_mid(1,0)) - ((C_list_mid.at(0).C / temp_h) * solution(1,0));
    // q = C*v
    C_list_mid.at(0).charge = C_list_mid.at(0).C * solution_mid(1,0);
    solution_mid.print("Solution_mid:");
   
    // temp_RHS_up.print("RHS_up matrix =");
    solution_up = NewtonRaphson_system(temp_LHS_up, temp_RHS_up, solution, history_voltages, next_h_up, mode);  // Soluition with larger time step
    // I = C/h*v(n+1) - C/h*v(n)
    C_list_up.at(0).current = ((C_list_up.at(0).C / next_h_up) * solution_up(1,0)) - ((C_list_up.at(0).C / next_h_up) * solution(1,0));
    C_list_up.at(0).temp_current_up.push_back(C_list_up.at(0).current);
    // q = C*v
    C_list_up.at(0).charge = C_list_up.at(0).C * solution_up(1,0);
    C_list_up.at(0).temp_charge_up.push_back(C_list_up.at(0).charge);
    
    solution_up.print("Solution_up:");

    // temp_RHS_down.print("RHS_down matrix =");
    solution_down = NewtonRaphson_system(temp_LHS_down, temp_RHS_down, solution, history_voltages, next_h_down, mode);  // Soluition with smaller time step
    // I = C/h*v(n+1) - C/h*v(n)
    C_list_down.at(0).current = ((C_list_down.at(0).C / next_h_down) * solution_down(1,0)) - ((C_list_down.at(0).C / next_h_down) * solution(1,0));
    // q = C*v
    C_list_down.at(0).charge = C_list_down.at(0).C * solution_down(1,0);
    solution_down.print("Solution_down:");


    // Calculate the difference between the solutions with different time steps
    arma::mat voltage_mid = solution_mid.submat(0, 0, T_nodes - 1, 0);
    arma::mat voltage_up = solution_up.submat(0, 0, T_nodes - 1, 0);
    arma::mat voltage_down = solution_down.submat(0, 0, T_nodes - 1, 0);
    arma::mat solution_voltages = solution.submat(0, 0, T_nodes - 1, 0);
    // x.v_mid = arma::abs(voltage_mid - solution_voltages).max();
    // x.v_up = arma::abs(voltage_up - solution_voltages).max();
    // x.v_down = arma::abs(voltage_down - solution_voltages).max();
    // x.v_mid_max = arma::max(voltage_mid,solution_voltages).max();
    // x.v_up_max = arma::max(voltage_up,solution_voltages).max();
    // x.v_down_max = arma::max(voltage_down,solution_voltages).max();
    

    // Calculate the difference between the solutions with different time steps
    arma::mat current_mid = solution_mid.submat(T_nodes, 0, solution_mid.n_rows - 1, 0);
    arma::mat current_up = solution_up.submat(T_nodes, 0, solution_up.n_rows - 1, 0);
    arma::mat current_down = solution_down.submat(T_nodes, 0, solution_down.n_rows - 1, 0);
    arma::mat solution_current = solution.submat(T_nodes, 0, solution.n_rows - 1, 0);
    // x.i_mid = arma::abs(current_mid - solution_current).max();
    // x.i_up = arma::abs(current_up - solution_current).max();
    // x.i_down = arma::abs(current_down - solution_current).max();
    // x.i_mid_max = arma::max(current_mid,solution_current).max();
    // x.i_up_max = arma::max(current_up,solution_current).max();
    // x.i_down_max = arma::max(current_down,solution_current).max();

    x.t_mid = temp_h;
    x.t_up = next_h_up;
    x.t_down = next_h_down;

    x.RHS_mid = temp_RHS_mid;
    x.RHS_up = temp_RHS_up;
    x.RHS_down = temp_RHS_down;

    x.LHS_mid = temp_LHS_mid;
    x.LHS_up = temp_LHS_up;
    x.LHS_down = temp_LHS_down;

    x.solution_mid = solution_mid;
    x.solution_up = solution_up;
    x.solution_down = solution_down;

   
    return x;
}

/* 
    0 = temp_h / 8
    1 = down
    2 = mid
    3 = up
    4 = temp_h * 8
    5 = previous up
*/ 
// The wrapper function of multi-solution method
arma::mat multi_next_h(arma::mat const init_LHS, arma::mat const init_RHS, arma::mat &LHS, arma::mat & RHS, arma::mat solution, 
                        std::deque<arma::vec> &history_voltages, double &h, std::deque<double> &history_steps, int mode, 
                        std::vector<double> RHS_locate, double V1, double V2, double & t1_pulse, double td, double tr, double tf, double tpw, double tper, int Maxi, double time_trans){

    std::cout << "Enter multi_next_h" << std::endl;
 
    double temp_h = h;
    int index_h = 0;  // 0 is temp_h / 8, 1 is down, 2 is mid, 3 is up, 4 is temp_h * 8, 5 is previous up
    int pre_decision = 0;
    Truncation_error LTE;
    //TMAX_reach = false;

    

    x_mid_up_down mid_up_down{};
    do{
        TMAX_reach = false;
        arma::mat temp_solution = solution;
        // solution with difference time step options
        arma::mat solution_mid{};
        arma::mat solution_up{};
        arma::mat solution_down{};
        // RHS for difference time step options
        arma::mat temp_RHS = init_RHS;
        arma::mat temp_RHS_mid = init_RHS;
        arma::mat temp_RHS_up = init_RHS;
        arma::mat temp_RHS_down = init_RHS;

        //LHS for difference time step options
        //arma::mat temp_LHS = init_LHS;
        arma::mat temp_LHS_mid = init_LHS;
        arma::mat temp_LHS_up = init_LHS;
        arma::mat temp_LHS_down = init_LHS;

    
        mid_up_down = multi_solution_solver(temp_LHS_mid, temp_LHS_up, temp_LHS_down, init_RHS, LHS, temp_solution, history_voltages, solution_mid, solution_up, solution_down, RHS_locate, temp_RHS_mid, temp_RHS_up, 
            temp_RHS_down, temp_h, V1, V2, t1_pulse, td, tr, tf, tpw, tper, mode, time_trans, mid_up_down);
        // if(ITL4_reach){
        //     temp_h = temp_h / 8;
        //     continue;
        // }
        
        // LTE current = h(k)*(lteabstol+lteretol*max(|ik+1|,|ik|))
        double last_step = history_steps.front();
        LTE.LTE_current_mid = last_step * (ABSTOL + RELTOL * std::max(std::abs(C_list_mid.at(0).pre_current), std::abs(C_list_mid.at(0).current)));
        LTE.LTE_current_up = last_step * (ABSTOL + RELTOL * std::max(std::abs(C_list_up.at(0).pre_current), std::abs(C_list_up.at(0).current)));
        LTE.LTE_current_down = last_step * (ABSTOL + RELTOL * std::max(std::abs(C_list_down.at(0).pre_current), std::abs(C_list_down.at(0).current)));
        // LTE current in spice book: e = RELTOL*max(|ik+1|,|ik|) + ABSTOL
        // LTE.LTE_current_mid = 1 * (ABSTOL + RELTOL * std::max(std::abs(C_list_mid.at(0).pre_current), std::abs(C_list_mid.at(0).current)));
        // LTE.LTE_current_up = 1 * (ABSTOL + RELTOL * std::max(std::abs(C_list_up.at(0).pre_current), std::abs(C_list_up.at(0).current)));
        // LTE.LTE_current_down = 1 * (ABSTOL + RELTOL * std::max(std::abs(C_list_down.at(0).pre_current), std::abs(C_list_down.at(0).current)));

        // LTE charge = relq * max(|qk+1|,|qk|, chgtol)
        LTE.LTE_charge_mid = relq * std::max({std::abs(C_list_mid.at(0).charge), std::abs(C_list_mid.at(0).pre_charge), CHGTOL});
        LTE.LTE_charge_up = relq * std::max({std::abs(C_list_up.at(0).charge), std::abs(C_list_up.at(0).pre_charge), CHGTOL});
        LTE.LTE_charge_down = relq * std::max({std::abs(C_list_down.at(0).charge), std::abs(C_list_down.at(0).pre_charge), CHGTOL});

        // LTE charge in spice book: e = RELTOL*max(|qk+1|,|qk|, chgtol) / hn
        // LTE.LTE_charge_mid = RELTOL * std::max({std::abs(C_list_mid.at(0).charge), std::abs(C_list_mid.at(0).pre_charge),CHGTOL}) / last_step;
        // LTE.LTE_charge_up = RELTOL * std::max({std::abs(C_list_up.at(0).charge), std::abs(C_list_up.at(0).pre_charge), CHGTOL}) / last_step;
        // LTE.LTE_charge_down = RELTOL * std::max({std::abs(C_list_down.at(0).charge), std::abs(C_list_down.at(0).pre_charge), CHGTOL}) / last_step;

        // LTE bound = max(LTE_current, LTE_charge)
        LTE.LTE_bound_mid = std::min(LTE.LTE_current_mid, LTE.LTE_charge_mid);
        LTE.LTE_bound_up = std::min(LTE.LTE_current_up, LTE.LTE_charge_up);
        LTE.LTE_bound_down = std::min(LTE.LTE_current_down, LTE.LTE_charge_down);
        // LTE BE = |h^2/2 * x''|     x'' = 2!DD2
        double DD0 = (solution_mid(1,0) * C_list_mid.at(0).C - solution(1,0) * C_list_mid.at(0).C) / mid_up_down.t_mid;
        double DD1 = (solution(1,0) * C_list_mid.at(0).C - history_voltages.at(1)(1,0) * C_list_mid.at(0).C) / history_steps.at(0);
        double DD2 = (DD0 - DD1) / (mid_up_down.t_mid + history_steps.at(0));
        LTE.LTE_BE_mid = std::abs(std::pow(mid_up_down.t_mid,2)/2 * std::max(lteabstol*last_step, std::abs(2*DD2))); // std::max(lteabstol*last_step, 2*DD2)

        DD0 = (solution_up(1,0) * C_list_up.at(0).C - solution(1,0) * C_list_up.at(0).C) / mid_up_down.t_up;
        DD1 = (solution(1,0) * C_list_up.at(0).C - history_voltages.at(1)(1,0) * C_list_up.at(0).C) / history_steps.at(0);
        DD2 = (DD0 - DD1) / (mid_up_down.t_up + history_steps.at(0));
        LTE.LTE_BE_up = std::abs(std::pow(mid_up_down.t_up,2)/2 * std::max(lteabstol*last_step, std::abs(2*DD2)));

        DD0 = (solution_down(1,0) * C_list_down.at(0).C - solution(1,0) * C_list_down.at(0).C) / mid_up_down.t_down;
        DD1 = (solution(1,0) * C_list_down.at(0).C - history_voltages.at(1)(1,0) * C_list_down.at(0).C) / history_steps.at(0);
        DD2 = (DD0 - DD1) / (mid_up_down.t_down + history_steps.at(0));
        LTE.LTE_BE_down = std::abs(std::pow(mid_up_down.t_down,2)/2 * std::max(lteabstol*last_step, std::abs(2*DD2)));

        index_h = h_reach_new(LTE, pre_decision);
        /*  index_h:
            0 = temp_h / 2
            1 = down
            2 = mid
            3 = up
            4 = temp_h * 2
            5 = previous up
        */ 
        switch(index_h){
        case 0:
            temp_h = temp_h / 2;
            mid_up_down.pre_t_up = mid_up_down.t_up;
            mid_up_down.pre_RHS_up = mid_up_down.RHS_up;
            mid_up_down.pre_LHS_up = mid_up_down.LHS_up;
            mid_up_down.pre_solution_up = mid_up_down.solution_up;
            break;
        case 1:                                                         
            h = mid_up_down.t_down;
            RHS = mid_up_down.RHS_down;
            LHS = mid_up_down.LHS_down;
            solution = mid_up_down.solution_down;
            C_list_down.at(0).pre_current = C_list_down.at(0).current;
            C_list_down.at(0).pre_charge = C_list_down.at(0).charge;
            C_list_mid.at(0).pre_current =  C_list_down.at(0).current;
            C_list_mid.at(0).pre_charge =  C_list_down.at(0).charge;
            C_list_up.at(0).pre_current =  C_list_down.at(0).current;
            C_list_up.at(0).pre_charge =  C_list_down.at(0).charge;
            break;
        case 2:
            h = mid_up_down.t_mid;
            RHS = mid_up_down.RHS_mid;
            LHS = mid_up_down.LHS_mid;
            solution = mid_up_down.solution_mid;
            C_list_mid.at(0).pre_current = C_list_mid.at(0).current;
            C_list_mid.at(0).pre_charge = C_list_mid.at(0).charge;
            C_list_up.at(0).pre_current = C_list_mid.at(0).current;
            C_list_up.at(0).pre_charge = C_list_mid.at(0).charge;
            C_list_down.at(0).pre_current = C_list_mid.at(0).current;
            C_list_down.at(0).pre_charge = C_list_mid.at(0).charge;
            break;
        case 3:
            h = mid_up_down.t_up;
            RHS = mid_up_down.RHS_up;
            LHS = mid_up_down.LHS_up;
            solution = mid_up_down.solution_up;
            C_list_up.at(0).pre_current = C_list_up.at(0).current;
            C_list_up.at(0).pre_charge = C_list_up.at(0).charge;
            C_list_mid.at(0).pre_current =  C_list_up.at(0).current;
            C_list_mid.at(0).pre_charge =  C_list_up.at(0).charge;
            C_list_down.at(0).pre_current =  C_list_up.at(0).current;
            C_list_down.at(0).pre_charge =  C_list_up.at(0).charge;
            break;
        case 4:
            temp_h = temp_h * 2;           
            mid_up_down.pre_t_up = mid_up_down.t_up;
            mid_up_down.pre_RHS_up = mid_up_down.RHS_up;
            mid_up_down.pre_LHS_up = mid_up_down.LHS_up;
            mid_up_down.pre_solution_up = mid_up_down.solution_up;
            break;
        case 5:
            h = mid_up_down.pre_t_up;
            RHS = mid_up_down.pre_RHS_up;
            LHS = mid_up_down.pre_LHS_up;
            solution = mid_up_down.pre_solution_up;
            C_list_up.at(0).pre_current = C_list_up.at(0).temp_current_up.at(C_list_up.at(0).temp_current_up.size() - 1);
            C_list_up.at(0).pre_charge = C_list_up.at(0).temp_charge_up.at(C_list_up.at(0).temp_charge_up.size() - 1);
            C_list_mid.at(0).pre_current =  C_list_up.at(0).pre_current;
            C_list_mid.at(0).pre_charge =  C_list_up.at(0).pre_charge;
            C_list_down.at(0).pre_current =  C_list_up.at(0).pre_current;
            C_list_down.at(0).pre_charge =  C_list_up.at(0).pre_charge;
            break;
        default:
            std::cout << "index_h is wrong in multi_next_h function." << "The value of index_h is: " << index_h <<std::endl;
            exit(1);
        }




    }while(index_h == 0 || index_h == 4);

    pre_decision = 0;
   


    // std::cout << "Fianl solution is:" << solution << std::endl;
    // std::cout << "Final time step is:" << h << std::endl;

    return solution;
}

void update_preClist(arma::mat solution, int tran_count, double h, std::deque<arma::vec> history_voltages){
    if(tran_count == 0){
        C_list_mid.at(0).pre_current = ((C_list_mid.at(0).C / h) * solution(1,0)) - ((C_list_mid.at(0).C / h) * history_voltages.at(0)(1,0));
        C_list_up.at(0).pre_current = C_list_mid[0].pre_current;
        C_list_down.at(0).pre_current = C_list_mid[0].pre_current;

        C_list_mid.at(0).pre_charge = C_list_mid.at(0).C * solution(1,0);
        C_list_up.at(0).pre_charge = C_list_mid[0].pre_charge;
        C_list_down.at(0).pre_charge = C_list_mid[0].pre_charge;
    }

    else{
        std::cout << "tran_count is wrong in update_preClist function." << "The value of tran_count is: " << tran_count <<std::endl;
        exit(1);
    }
}



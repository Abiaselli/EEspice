#ifndef Transient_code
#define Transient_code
#include <iostream>
#include "armadillo"
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
// double RELTOL = 1e-3;
double RELTOL = 1e-4;

/* Absolute error for voltage*/
double VNTOL = 1e-6;
/* Absolute error for current*/
double ABSTOL = 1e-12;
// double ABSTOL = 1e-13;

/* Min and Max timestep */
double TMIN = 1e-18;
double TMAX = t_end / 10; 

/* Threshold for iteration count-based algorithm */
int ITL3 = 4;
int ITL4 = 10;

/*If timestep reach the limit*/
bool TMAX_reach = false;
bool TMIN_reach = false;

// History timestep list
std::vector<double> h_history;      // Practical timestep
std::vector<double> h_1_history;    // Expected timestep

struct x_mid_up_down;

//History Current List (nodex, nodey, Currents) for time-related elements
struct CapacitorState;
std::vector<CapacitorState> capacitorStates;

// Matrix stamps assigner using Modified Nodal Analysis
std::pair<arma::mat, arma::mat> DynamicNonLinear(arma::mat &LHS, arma::mat &RHS,
                                                                       arma::mat solution, std::deque<arma::vec> &history_voltages,
                                                                       double h, int mode);
void R_assigner(double node_x, double node_y, double R, arma::mat &LHS, arma::mat &RHS);
void Is_assigner(double node_x, double node_y, double I, arma::mat &LHS, arma::mat &RHS);
double Vs_assigner(int node_x, int node_y, double V_value, arma::mat &LHS, arma::mat &RHS);
void C_assigner(int node_x, int node_y, double C, double h, arma::mat &LHS, arma::mat &RHS, arma::mat solution, int mode);
void Diode_assigner(int node_x, int node_y, double Is, double VT, double cd, double h, arma::mat &LHS, arma::mat &RHS, arma::mat solution, int mode);
void VCCS_assigner(int node_x, int node_y, int node_cx, int node_cy, double R, arma::mat &LHS);
void C_assigner_3(int node_x, int node_y, double C, double h, arma::mat &LHS, arma::mat &RHS,
                  std::deque<arma::mat> &NR_solutions,
                  int mode);
double Diode_assigner_2(int node_x, int node_y, double Is, double VT, double h, arma::mat &LHS, arma::mat &RHS, 
                      arma::mat solution, int mode);
bool isConverge(std::deque<arma::mat> &NR_solutions);

void UpdateStates(arma::mat &LHS, arma::mat &RHS, 
                    std::deque<arma::mat> &NR_solutions, 
                    double h, int mode);

void history_voltages_update(arma::mat &solution, std::deque<arma::vec> history_voltages);
arma::mat NewtonRaphson_system(arma::mat &LHS, arma::mat &RHS, arma::mat solution, std::deque<arma::vec> &history_voltages,
                               double &h, int mode, std::vector<CapacitorState> &capacitorStates);
void clear_capacitorStates(std::vector<CapacitorState> &capacitorStates);
x_mid_up_down  multi_solution_solver(arma::mat &temp_LHS_mid, arma::mat &temp_LHS_up, arma::mat &temp_LHS_down, arma::mat const init_RHS, arma::mat LHS, arma::mat solution, std::deque<arma::vec> &history_voltages, 
                            arma::mat &solution_mid, arma::mat &solution_up, arma::mat &solution_down, std::vector<double> RHS_locate, arma::mat &temp_RHS_mid, 
                            arma::mat &temp_RHS_up, arma::mat &temp_RHS_down, double & temp_h, double V1, double V2, double & t1_pulse, double td, double tr, 
                            double tf, double tpw, double tper, int mode, double time_trans, std::vector<CapacitorState> &capacitorStates_mid, 
                            std::vector<CapacitorState> &capacitorStates_up, std::vector<CapacitorState> &capacitorStates_down);


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




struct CapacitorState {
    double pre_x = 0;
    double pre_x1 = 0;
    double pre_vol = 0;
};


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
        //pre_solution.print("pre_solution in c is ");

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
    // state.pre_x = x;
    // state.pre_x1 = x1;
    // state.pre_vol = vol;
    
}
// Voltage pulse assigner
double V_pulse(double V1, double V2, double &t1, double td, double tr, double tf, double tpw, double tper, double h)
{
    double V_t1 = V1;
    if ((t1 >= 0) && (t1 < td))
    {
        V_t1 = V1;
    }
    else if ((t1 >= td) && (t1 < (td + tr)))
    {
        V_t1 = V1 + (V2 - V1) * (t1 - td) / tr;
    }
    else if ((t1 >= (td + tr)) && (t1 < (td + tr + tpw)))
    {
        V_t1 = V2;
    }
    else if ((t1 >= (td + tr + tpw)) && (t1 < (td + tr + tpw + tf)))
    {
        V_t1 = V2 + (V1 - V2) * (t1 - (td + tr + tpw)) / tf;
    }
    else if ((t1 >= (td + tr + tpw + tf)) && (t1 <= (td + tper)))
    {
        V_t1 = V1;
    }
    else
    {
        t1 = td;
    }
    t1 = t1 + h;

    return V_t1;
}

// Voltage pulse assigner
double V_pulse_new(double V1, double V2, double t1, double td, double tr, double tf, double tpw, double tper)
{   
    double tnorm = fmod((t1-td),tper);
    double v = 0;
    const double EPSILON = 1e-18;

    // if(std::fabs(tper - (tr + tpw + tf)) > EPSILON) {
    //     std::cout << tper << std::endl;
    //     std::cout << tr + tpw + tf << std::endl;
    //     std::cout << "Pulse voltage period is wrong" << std::endl;
    //     exit(1);
    // }

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
                               double &h, int mode, std::vector<CapacitorState> &capacitorStates)
{
    std::cout << "Enter NewtonRaphson system" << std::endl;

    int col_size = LHS.n_cols;
    int row_size = LHS.n_rows;
    int iteration_counter = 0;
    bool isconverge = false;
    std::deque<arma::mat> NR_solutions;
    NR_solutions.push_back(history_voltages.at(0));  //input the solution from the OP analysis
    arma::mat NR_RHS;
    arma::mat NR_LHS;

    // using relative error to check for convergence
    // Note that pre_solution is just the temporary solution or voltage during the NR iteration
    do
    {   //std::cout<<"iteration_counter is "<< iteration_counter <<std::endl;
        NR_LHS = LHS;
        NR_RHS = RHS;
        // Update the model of time-related elements (Capacitor)
        UpdateStates(NR_LHS, NR_RHS, NR_solutions, h, mode);

        // Update the non-linear elements (Diode, MOSFET) with updated solution
        auto matrices = DynamicNonLinear(NR_LHS, NR_RHS, solution, history_voltages, h, mode);
        // std::cout << "h_111: " << h << std::endl;
        solution = arma::solve(matrices.first, matrices.second);    // solution is a list of all the current node voltages
        iteration_counter += 1;

        //update the NR_solutions(current_solution, pre_solution, next_solution) 
        NR_solutions.push_back(solution);
        if (NR_solutions.size() > 3)
        {
            NR_solutions.pop_front();
        }

        // Check for convergence only if you have 3 solutions
        if (NR_solutions.size() == 3) {
            isconverge = isConverge(NR_solutions);
        }


        if (iteration_counter > 100)
        {
            std::cout << "Not Converge at 100 iterations" << std::endl;
            break;
        }

    } while (!isconverge);

    // NR_LHS.print("NR_LHS is ");
    // NR_RHS.print("NR_RHS is ");

    LHS = NR_LHS;
    RHS = NR_RHS;

    // print iteration count
    std::cout << "The NR is converge and Iteration count is: " << iteration_counter << std::endl;

    return solution;
}

void history_voltages_update(arma::mat &solution, std::deque<arma::vec> history_voltages)
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
    if(pre_solution.n_rows != current_solution.n_rows || pre_solution.n_rows != next_solution.n_rows || current_solution.n_rows != next_solution.n_rows){
        std::cout << "The size of pre_solution, current_solution, next_solution is not the same in isConverge function." << std::endl;
        exit(1);
    }
    if(pre_solution.n_rows != T_nodes + num_v){
        std::cout << "The size of pre_solution is not the same as T_nodes + num_v in isConverge function." << std::endl;
        exit(1);
    }

    arma::mat pre_voltages = pre_solution.submat(0, 0, T_nodes-1, 0);
    arma::mat current_voltages = current_solution.submat(0, 0, T_nodes-1, 0);
    arma::mat next_voltages = next_solution.submat(0, 0, T_nodes-1, 0);

    arma::mat pre_current = pre_solution.submat(T_nodes, 0, T_nodes + num_v - 1, 0);
    arma::mat current_current = current_solution.submat(T_nodes, 0, T_nodes + num_v - 1, 0);
    arma::mat next_current = next_solution.submat(T_nodes, 0, T_nodes + num_v - 1, 0);


    // |v(k+1)-v(k)| <= RELTOL*max(|v(k+1)|,|v(k)|) + VNTOL
    for (int i = 0; i < next_voltages.n_rows; i++)
    {
        double maxVoltage = std::max(std::abs(next_voltages(i, 0)), std::abs(current_voltages(i, 0)));
        double Veps = std::abs(next_voltages(i, 0) - current_voltages(i, 0));

        if (Veps > RELTOL * maxVoltage + VNTOL)
        {
            return false;
        }
    }
 
    // |i(k+1)-i(k)| <= RELTOL*max(|i(k+1)|,|i(k)|) + ABSTOL
    for (int i = 0; i < next_current.n_rows; i++)
    {
        double maxCurrent = std::max(std::abs(next_current(i, 0)), std::abs(current_current(i, 0)));
        double Ieps = std::abs(next_current(i, 0) - current_current(i, 0));

        if (Ieps > RELTOL * maxCurrent + ABSTOL)
        {
            return false;
        }
    }

    // |v(k) - v(k-1)| <= RELTOL*max(|v(k)|,|v(k-1)|) + VNTOL
    for (int i = 0; i < current_voltages.n_rows; i++)
    {
       double maxVoltage = std::max(std::abs(current_voltages(i, 0)), std::abs(pre_voltages(i, 0)));
        double Veps = std::abs(current_voltages(i, 0) - pre_voltages(i, 0));

        if (Veps > RELTOL * maxVoltage + VNTOL)
        {
            return false;
        }
    }

    // |i(k) - i(k-1)| <= RELTOL*max(|i(k)|,|i(k-1)|) + ABSTOL
    for (int i = 0; i < current_current.n_rows; i++)
    {
        double maxCurrent = std::max(std::abs(current_current(i, 0)), std::abs(pre_current(i, 0)));
        double Ieps = std::abs(current_current(i, 0) - pre_current(i, 0));

        if (Ieps > RELTOL * maxCurrent + ABSTOL)
        {
            return false;
        }
    }

    // |v(k+1) - v(k-1)| <= √|v(k) - v(k-1)|^2 + |v(k+1) - v(k)|^2
    for (int i = 0; i < next_voltages.n_rows; i++)
    {
        double Veps = std::abs(next_voltages(i, 0) - pre_voltages(i, 0));

        if (Veps > std::sqrt(std::pow(std::abs(next_voltages(i, 0) - current_voltages(i, 0)), 2) + std::pow(std::abs(current_voltages(i, 0) - pre_voltages(i, 0)), 2)))
        {
            return false;
        }
    }

    // |i(k+1) - i(k-1)| <= √|i(k) - i(k-1)|^2 + |i(k+1) - i(k)|^2
    for (int i = 0; i < next_current.n_rows; i++)
    {
        double Ieps = std::abs(next_current(i, 0) - pre_current(i, 0));

        if (Ieps > std::sqrt(std::pow(std::abs(next_current(i, 0) - current_current(i, 0)), 2) + std::pow(std::abs(current_current(i, 0) - pre_current(i, 0)), 2)))
        {
            return false;
        }
    }

    return true;
}

// Adatpive timestep adjustment
/* STEPMODE = 1 is the count-based adjustment: still some buggy performance in this kind of adjustment */
/* STEPMODE = 2 is the LTE-based adjustment for BE method */

struct x_mid_up_down{
    double v_mid{};  // 
    double v_up{};
    double v_down{};
    
    double i_mid{};
    double i_up{};
    double i_down{};

    double t_mid{};
    double t_up{};
    double t_down{};   

    arma::mat RHS_mid{};
    arma::mat RHS_up{};
    arma::mat RHS_down{}; 

    arma::mat LHS_mid{};
    arma::mat LHS_up{};
    arma::mat LHS_down{};

    arma::mat solution_mid{};
    arma::mat solution_up{};
    arma::mat solution_down{};
};

// 0 is temp_h / 8, 1 is down, 2 is mid, 3 is up, 4 is temp_h * 8
int h_reach(x_mid_up_down mid_up_down , double &temp_h, int counter){
    
    int reach = 0;
    bool a = mid_up_down.v_up > VNTOL || mid_up_down.i_up > ABSTOL;
    bool b =mid_up_down.v_mid > VNTOL || mid_up_down.i_mid > ABSTOL;
    bool c = mid_up_down.v_down > VNTOL || mid_up_down.i_down > ABSTOL;
    
    if(mid_up_down.v_up > VNTOL || mid_up_down.i_up > ABSTOL){
        
        reach = 0;
        return reach;
    }
    else if(mid_up_down.v_mid > VNTOL || mid_up_down.i_mid > ABSTOL){
        reach = 3;
        return reach;
    }
    else if(mid_up_down.v_down > VNTOL || mid_up_down.i_down > ABSTOL){
        reach = 2;
        return reach;
    }
    else if(counter > 2){
        reach = 3;
        return reach;
    }else if(mid_up_down.v_up > VNTOL || mid_up_down.i_up > ABSTOL){
        reach = 4;
        return reach;
    }

    // if (c){
    // }A
    // if(!a){}
}



// New timestep options
void timestep_options(double & temp_h, double & next_h_up, double & next_h_down) {
    if(temp_h > TMAX) {
        temp_h = TMAX;
        TMAX_reach = true;
    }
    if(temp_h < TMIN) {
        temp_h = TMIN;
        TMIN_reach = true;
    }

    next_h_up = temp_h * 2;
    next_h_down = temp_h / 2;

    if(next_h_up > TMAX) {
        next_h_up = TMAX;
        TMAX_reach = true;
    }
    if(next_h_up < TMIN) {
        next_h_down = TMIN;
        TMIN_reach = true;
    }
    if(next_h_down > TMAX) {
        next_h_up = TMAX;
        TMAX_reach = true;
    }
    if(next_h_down < TMIN) {
        next_h_down = TMIN;
        TMIN_reach = true;
    }
}

x_mid_up_down  multi_solution_solver(arma::mat &temp_LHS_mid, arma::mat &temp_LHS_up, arma::mat &temp_LHS_down, arma::mat const init_RHS, arma::mat LHS, arma::mat solution, std::deque<arma::vec> &history_voltages, 
                            arma::mat &solution_mid, arma::mat &solution_up, arma::mat &solution_down, std::vector<double> RHS_locate, arma::mat &temp_RHS_mid, 
                            arma::mat &temp_RHS_up, arma::mat &temp_RHS_down, double & temp_h, double V1, double V2, double & t1_pulse, double td, double tr, 
                            double tf, double tpw, double tper, int mode, double time_trans, std::vector<CapacitorState> &capacitorStates_mid, 
                            std::vector<CapacitorState> &capacitorStates_up, std::vector<CapacitorState> &capacitorStates_down) {
    std::cout << "multi_solution_solver" << std::endl;
    double next_h_up;
    double next_h_down;

    
    timestep_options(temp_h, next_h_up, next_h_down);

    std::cout << "temp_h:" << temp_h ;
    std::cout << " next_h_up:" << next_h_up;
    std::cout << " next_h_down:" << next_h_down << std::endl;
    //  exit(0);
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    double t1_mid = time_trans + temp_h;
    double t1_up = time_trans + next_h_up;
    double t1_down = time_trans + next_h_down;
    std::cout << "t1_mid" << t1_mid;
    std::cout << " t1_up" << t1_up;
    std::cout << " t1_down" << t1_down << std::endl;
    

    // std::cout << "The time steps are:" << next_h_up << " " << temp_h << " " << next_h_down << std::endl;

    std::vector<double> RHS_value = {
        V_pulse_new(V1, V2, t1_mid, td, tr, tf, tpw, tper)
    };
    // std::cout << "RHS_value_mid" << RHS_value[0] << std::endl;

    std::vector<double> RHS_value_up = {
        V_pulse_new(V1, V2, t1_up, td, tr, tf, tpw, tper)
    };
    // std::cout << "RHS_value_up" << RHS_value_up[0] << std::endl;
    std::vector<double> RHS_value_down = {
        V_pulse_new(V1, V2, t1_down, td, tr, tf, tpw, tper)
    };
    // std::cout << "RHS_value_down" << RHS_value_down[0] << std::endl;

    temp_RHS_mid = RHS_update(RHS_locate, init_RHS, RHS_value);

    temp_RHS_up = RHS_update(RHS_locate, init_RHS, RHS_value_up);

    temp_RHS_down = RHS_update(RHS_locate, init_RHS, RHS_value_down);

    history_voltages_update(solution, history_voltages);


    temp_RHS_mid.print("RHS_mid matrix =");
    solution_mid = NewtonRaphson_system(temp_LHS_mid, temp_RHS_mid, solution, history_voltages, temp_h, mode, capacitorStates_mid); // Soluition with current time step
    std::cout << "Solution_mid:"<< std::endl << solution_mid << std::endl;
   
    temp_RHS_up.print("RHS_up matrix =");
    solution_up = NewtonRaphson_system(temp_LHS_up, temp_RHS_up, solution, history_voltages, next_h_up, mode, capacitorStates_up);  // Soluition with larger time step
    std::cout << "Solution_up:" << std::endl << solution_up << std::endl;

    temp_RHS_down.print("RHS_down matrix =");
    solution_down = NewtonRaphson_system(temp_LHS_down, temp_RHS_down, solution, history_voltages, next_h_down, mode, capacitorStates_down);  // Soluition with smaller time step
    std::cout << "Solution_down:" << std::endl << solution_down << std::endl;

    x_mid_up_down x;

    // Calculate the difference between the solutions with different time steps
    arma::mat voltage_mid = solution_mid.submat(0, 0, T_nodes - 1, 0);
    arma::mat voltage_up = solution_up.submat(0, 0, T_nodes - 1, 0);
    arma::mat voltage_down = solution_down.submat(0, 0, T_nodes - 1, 0);
    arma::mat solution_voltages = solution.submat(0, 0, T_nodes - 1, 0);
    x.v_mid = arma::abs(voltage_mid - solution_voltages).max();
    x.v_up = arma::abs(voltage_up - solution_voltages).max();
    x.v_down = arma::abs(voltage_down - solution_voltages).max();

    // Calculate the difference between the solutions with different time steps
    arma::mat current_mid = solution_mid.submat(T_nodes, 0, T_nodes + num_v - 1, 0);
    arma::mat current_up = solution_up.submat(T_nodes, 0, T_nodes + num_v - 1, 0);
    arma::mat current_down = solution_down.submat(T_nodes, 0, T_nodes + num_v - 1, 0);
    arma::mat solution_current = solution.submat(T_nodes, 0, T_nodes + num_v - 1, 0);
    x.i_mid = arma::abs(current_mid - solution_current).max();
    x.i_up = arma::abs(current_up - solution_current).max();
    x.i_down = arma::abs(current_down - solution_current).max();

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
    //exit(0);
}



// The wrapper function of multi-solution method
arma::mat multi_next_h(arma::mat const init_LHS, arma::mat const init_RHS, arma::mat &LHS, arma::mat & RHS, arma::mat solution, 
                        std::deque<arma::vec> &history_voltages, double &h, std::deque<double> &history_steps, int mode, 
                        std::vector<double> RHS_locate, double V1, double V2, double & t1_pulse, double td, double tr, double tf, double tpw, double tper, int Maxi, double time_trans){

    std::cout << "Enter multi_next_h" << std::endl;
 
    double temp_h = h;
    int index_h = 0;  // 0 is false, 1 is down, 2 is mid, 3 is up
    int counter = 0;
    x_mid_up_down mid_up_down{};
    do{
        clear_capacitorStates(capacitorStates);
        // solution with difference time step options
        arma::mat solution_mid{};
        arma::mat solution_up{};
        arma::mat solution_down{};
        // RHS for difference time step options
        arma::mat temp_RHS = init_RHS;
        arma::mat temp_RHS_mid = init_RHS;
        arma::mat temp_RHS_up = init_RHS;
        arma::mat temp_RHS_down = init_RHS;
        // Capacitor state for difference time step options
        std::vector<CapacitorState> capacitorStates_up = capacitorStates;
        std::vector<CapacitorState> capacitorStates_down = capacitorStates;
        std::vector<CapacitorState> capacitorStates_mid = capacitorStates;
        //LHS for difference time step options
        //arma::mat temp_LHS = init_LHS;
        arma::mat temp_LHS_mid = init_LHS;
        arma::mat temp_LHS_up = init_LHS;
        arma::mat temp_LHS_down = init_LHS;

        history_voltages_update(solution, history_voltages);
    
        mid_up_down = multi_solution_solver(temp_LHS_mid, temp_LHS_up, temp_LHS_down, init_RHS, LHS, solution, history_voltages, solution_mid, solution_up, solution_down, RHS_locate, temp_RHS_mid, temp_RHS_up, 
            temp_RHS_down, temp_h, V1, V2, t1_pulse, td, tr, tf, tpw, tper, mode, time_trans, capacitorStates_mid, capacitorStates_up, capacitorStates_down);

        index_h = h_reach(mid_up_down,temp_h, counter);
        counter += 1;

        // 0 is temp_h / 8, 1 is down, 2 is mid, 3 is up, 4 is temp_h * 8
        switch(index_h){
        case 0:
            temp_h = temp_h / 8;
            break;
        case 1:
            h = mid_up_down.t_down;
            RHS = mid_up_down.RHS_down;
            LHS = mid_up_down.LHS_down;
            solution = mid_up_down.solution_down;
            break;
        case 2:
            h = mid_up_down.t_mid;
            RHS = mid_up_down.RHS_mid;
            LHS = mid_up_down.LHS_mid;
            solution = mid_up_down.solution_mid;
            break;
        case 3:
            h = mid_up_down.t_up;
            RHS = mid_up_down.RHS_up;
            LHS = mid_up_down.LHS_up;
            solution = mid_up_down.solution_up;
            break;
        case 4:
            temp_h = temp_h * 8;
            break;
        default:
            std::cout << "index_h is wrong in multi_next_h function." << "The value of index_h is: " << index_h <<std::endl;
            exit(1);
    }
    }while(index_h == 0 || index_h == 4);

    
    clear_capacitorStates(capacitorStates);


    std::cout << "Fianl solution is:" << solution << std::endl;
    std::cout << "Final time step is:" << h << std::endl;

    return solution;
}

void clear_capacitorStates(std::vector<CapacitorState> &capacitorStates) {
    for (auto &state : capacitorStates) {
        state.pre_x = 0;
        state.pre_x1 = 0;
        state.pre_vol = 0;
    }
}


#endif
/*  This code could run a Direct Current Operating Point (DC OP) Analysis and Transient Analysis of a circuit. The circuit could be consisting of
    resistors, capacitors, diodes, current source, voltage source, voltage-controlled current source (VCCS), pulsed voltage source, and N-MOS transistor.

    The assignments are given as following (with solution, LHS, and RHS being constant):
    1) Resistor - R_assigner(node_x, node_y, cond(R), LHS, RHS);
    2) Capacitors - C_assigner(node_x, node_y, C, h, LHS, RHS, solution, mode);
    3) Diodes - Diode_assigner( node_x, node_y, Is, VT, cd, h, LHS, RHS, solution, mode);
    4) Current Source - Is_assigner( node_x, node_y, I, LHS, RHS);
    5) Voltage Source - Vs_assigner(node_x, node_y, V_value, LHS, RHS);
    6) VCCS - VCCS_assigner(node_x, node_y, node_cx, node_cy, R, LHS);
    7) Pulsed Voltage Source - V_pulse(V1, V2, t1, td, tr, tf, tpw, tper, h);
    8) N-MOS Transistor - NMOS_assigner(number, node_vd, node_vg, node_vs, node_vb, h, solution, LHS, RHS, mode);
    9) P-MOS Transistor - PMOS_assigner(number, node_vs, node_vg, node_vd, node_vb, h, solution, LHS, RHS, mode);
    10) Ring Oscillator - RingOscillatorStages(W, L, R, C, LHS, RHS, solution, h, mode);

    The linear components could be assigned inside the main function while the non-linear and dynamic components can be assigned inside the DynamicNonLinear function.

    V_pulse also has a special assignment where the normal voltage source needs to be assigned with the initial voltage (V1) inside the RHS_locate list.
    The V_pulse function will then be called inside the transient simulation loop.
*/
int const code = 0; // choosing code for transient simulation

// Setting for the Mosfet
double W = 500e-6;
double L = 50e-6;

// Settings for the Ring Oscillator
double W_oscillator = 500e-6;
double L_oscillator = 50e-6;
double C_oscillator = 1e-15;
double R_oscillator = 1e3;
//int const cascaded_level = 3;      // Number of cascaded ring oscillators
int const cascaded_level = 0;      // Number of cascaded ring oscillators
int const supply_voltage_node = 1; // supply voltage node for the ring oscillator


// // TRANSIENT SIMULATION SETTINGS
// double t_start = 0;
// double t_end = 1e-4;
// double h = t_end / 2000; // t_end/2000 is the default value

// TRANSIENT SIMULATION SETTINGS
double t_start = 0;

double t_end = 1.2e-3;

double h = t_end/5000; // t_end/5000 is the default value
int num_v = 1; // number of voltage sources

int if_step_down = 0;
double last_step = 0;

int ref_node_x;
int ref_node_y;

/*  TOTAL NUMBER OF NODES EXCLUDING GROUND
    Two port components such as resistors, initially adds 2 nodes. If more than 1 component is added, it then adds 1 node per component.
    Number of MOSFETs and cascaded levels are assigned above too. External nodes are the nodes which  */

// int const external_nodes = 1;                                    // Number of external nodes (excluding ground and ring oscillator loop nodes)
// int const external_mosfets = 0;                                  // Number of standalone mosfets (excluding mosfets from ring oscillator)
// int const no_of_mosfets = external_mosfets + 2 * cascaded_level; // Total number of MOSFETs
// int const T_nodes = external_nodes + 4 * no_of_mosfets + 2 * cascaded_level;
int const external_nodes = 2;                                    // Number of external nodes (excluding ground and ring oscillator loop nodes)
int const external_mosfets = 0;                                  // Number of standalone mosfets (excluding mosfets from ring oscillator)
int const no_of_mosfets = 0; // Total number of MOSFETs
int const T_nodes = external_nodes + 4 * no_of_mosfets + 2 * cascaded_level;


#include "Transient_code.h"

/*  The line above the code means that the section can be changed by the user which analyses circuit simulation
P
    --fixed--           : can't be edited for circuit simulation purposes, but can be edited for code debugging.
    --can be changed--  : can be edited for circuit simulation purposes and edited for code debugging.

*/

// Assigning the stamp matrices for dynamic and non-linear components
std::pair<arma::mat, arma::mat> DynamicNonLinear(arma::mat &LHS, arma::mat &RHS,
                                                                    arma::mat solution, std::deque<arma::vec> &history_voltages, 
                                                                    double h, int mode)
{
    // (All the circuit assigners can be used except for voltage and current sources)
    /*--------------------------------------------can be changed-------------------------------------------------*/
    // RingOscillatorStages(W_oscillator, L_oscillator, R_oscillator, C_oscillator, LHS, RHS, solution, history_voltages, h, mode);
    // NMOS_assigner(1,1,2,3,3, W, L, h, solution, history_voltages, LHS, RHS, mode);

    arma::mat J_x = LHS;
    arma::mat Z_x = RHS;

    return {J_x, Z_x};
}

// Updates the stamp matrices for time-related components
void UpdateStates(arma::mat &LHS, arma::mat &RHS, 
                    std::deque<arma::mat> &NR_solutions, 
                    double h, int mode)
{
    // (Time-realted assigners (RingOscillatorStages, Capacitor, Inductor, etc) can be used except for voltage and current sources)
    /*--------------------------------------------can be changed-------------------------------------------------*/
    // RingOscillatorStages(W_oscillator, L_oscillator, R_oscillator, C_oscillator, LHS, RHS, solution, history_voltages, h, mode);
    // Voltage_divider(W_oscillator, L_oscillator, R_oscillator, C_oscillator, LHS, RHS, solution, history_voltages, h, mode);

    C_assigner_3(2, 0, 1e-9, h, LHS, RHS, NR_solutions, mode);
    // if(mode == 1){
    //     RHS.print("RHS matrix after UpdateStates =");
    //     LHS.print("LHS matrix after UpdateStates =");
        
    // }
    // NMOS_assigner(1,1,2,3,3, W, L, h, solution, history_voltages, LHS, RHS, mode);
}

// Main function for the circuit simulation
int main(int argc, const char **argv)
{
    auto t1 = std::chrono::high_resolution_clock::now(); // Start time
    /*----------------------------------------------fixed--------------------------------------------------------*/
    // Size of matrix
    int Maxi{T_nodes}; // defined in header file (Transient_code.h)
    int Maxj{Maxi};

    /*----------------------------------------------fixed--------------------------------------------------------*/
    // ASSIGNING THE STAMPS TO THE LHS AND RHS MATRICES

    // default state
    arma::mat LHS = arma::zeros(Maxi, Maxj); // LHS matrix
    arma::mat RHS = arma::zeros(Maxi, 1);    // RHS matrix
    arma::mat Last_LHS = arma::zeros(Maxi, Maxj); // LHS matrix
    arma::mat Last_RHS = arma::zeros(Maxi, 1);    // RHS matrix
    arma::mat Last_LHS_2 = arma::zeros(Maxi, Maxj); // LHS matrix
    arma::mat Last_RHS_2 = arma::zeros(Maxi, 1);    // RHS matrix

    /*--------------------------------------------can be changed-------------------------------------------------*/
    // ASSIGNING THE RESISTOR STAMP (R_assigner)
    // R_assigner(1,2,3e3,LHS,RHS);
    // R_assigner(1,2,3e3,LHS,RHS);
    // R_assigner(2,0,3e3,LHS,RHS);

    R_assigner(1,2,100000,LHS,RHS);
    // R_assigner(2,0,3,LHS,RHS);

    //R_assigner(3,0,3,LHS,RHS);

    /*--------------------------------------------can be changed-------------------------------------------------*/
    // ASSIGNING THE CURRENT STAMP (Is_assigner, VCCS_assigner)

    /*--------------------------------------------can be changed-------------------------------------------------*/
    // ASSIGNING THE VOLTAGE SOURCES (Vs_assigner)
    // Pulse voltage settings
    double t1_pulse = 0; // time used for the loop
    double V1 = 0;
    double V2 = 5;
    double td = 0;
    double tr = 1e-5;
    double tf = 1e-5;
    double tpw = 1e-3;
    double tper = 1.02e-3;

    std::vector<double> RHS_locate = {
        // Assigning the voltage matrix on LHS and RHS for the pulse voltage
        Vs_assigner(1, 0, V1, LHS, RHS)
        // Vs_assigner(2,0, V1, LHS, RHS)
    };
    // Assigning DC voltage sources
    // Vs_assigner(supply_voltage_node, 0, 5, LHS, RHS); // supply voltage for the ring oscillator, vdd

    // Assigning the stamps that would affect the RHS in transient simulation
    // (only for  time-dependent voltage, e.g. pulse voltages)

    /*----------------------------------------------fixed--------------------------------------------------------*/
    // Checking the LHS and RHS matrices
    LHS.print("init_LHS matrix =");
    RHS.print("init_RHS matrix =");
    /*----------------------------------------------fixed--------------------------------------------------------*/
    // OPERATING POINT ANALYSIS SYSTEM
    int mode = 0; // 0 to do OP analysis, 1 to do transient simulation
    // zero as initial condition
    Maxi = RHS.n_rows;
    Maxj = RHS.n_cols;

    // The initial LHS and RHS values to be used in the NR-algorithm
    arma::mat const init_LHS = LHS;
    arma::mat const init_RHS = RHS;

    // OP analysis used as initial condition for next evaluation
    arma::vec solution = arma::zeros(Maxi, Maxj);

    // Store the history voltages(solutions)
    std::deque<arma::vec> history_voltages;
    history_voltages.push_front(solution);
    
    // Store the history steps(solutions)
    std::deque<double> history_steps;
    history_steps.push_front(h);

    // Benchmarking for OP analysis
    auto tstart_op = std::chrono::high_resolution_clock::now();
    solution = NewtonRaphson_system(LHS, RHS, solution, history_voltages, h, mode);
    // std::cout << "DC_Solution:"<< std::endl << solution << std::endl;
    // exit(true);

    // add the solution(voltage) to the history
    history_voltages.push_front(solution);

    auto tstop_op = std::chrono::high_resolution_clock::now();

    solution.print("The OP analysis of the circuit is: ");
    
    /*----------------------------------------------fixed--------------------------------------------------------*/
    // The solution csv that is going to be plotted which contains the values of nodal voltages
    // and voltage source currents
    arma::vec solution_csv = solution;
    arma::vec Max_I = arma::zeros(1, 1);
    Max_I.row(0).col(0) = Maxi;
    arma::mat zero_ext = arma::zeros(1, 1);
    /*--------------------------------------------can be changed-------------------------------------------------*/
    // ADDING TRANSIENT SIMULATION LOOP (includes V_pulse or any time dependent sources)
    double time_trans = t_start;
    mode = 1;
    int count = 0;
    int tran_count = 0;
    auto tstart_trans = std::chrono::high_resolution_clock::now();

    history_voltages_update(solution, history_voltages);
    std::cout << "transient simulation start" << std::endl;
    while (time_trans < t_end)
    {

        LHS = init_LHS;
        RHS = init_RHS;

        std::cout << "-------------------------------------------------------" << std::endl;
        std::cout << "Step one" << std::endl;

        solution = multi_next_h(init_LHS, init_RHS, LHS, RHS, solution, history_voltages, h, history_steps, mode, RHS_locate, V1, V2, t1_pulse, td, tr, tf,  tpw, tper, Maxi, time_trans);
        solution.print("The solution is:");
        history_voltages_update(solution, history_voltages); // save the solution to the history

        double step = h;

        // if(step == 0) exit(true);

        // Assigning the variables that will be plotted and analysed as seen in a circuit simulator
        solution_csv = arma::join_cols(solution_csv, solution);
        time_trans += step;
    
        // add the current step to the history steps
        history_steps.push_front(step);
        
        count += 1;

        // if(count == 1) exit(true);
        // std::cout<< "This is the :" << count  << step << std::endl;
        std::cout<< "current step is:" << step << std::endl;
        std::cout << "h_step size is:" << history_steps.size() << std::endl;
        std::cout<< "time_trans is: " << time_trans << std::endl;

    

        // if (history_steps.size() == 30)
        //     exit(true);
        tran_count = tran_count + 1;
    }
    
    auto tstop_trans = std::chrono::high_resolution_clock::now();
    // time vector to be inputted in plot for python analysis
    arma::mat time = arange(t_start, history_steps);

    /* Getting number of milliseconds as a double. */
    std::chrono::duration<double, std::milli> OP_time = (tstop_op - tstart_op);
    std::chrono::duration<double, std::milli> trans_time = (tstop_trans - tstart_trans);

    std::cout << "DC OP time:" << OP_time.count() << "ms\n";
    std::cout << "Transient time:" << trans_time.count() << "ms\n";

    /*-----------------------------------------------------------------------------------------------------------*/
    // SAVING THE SOLUTION AND TIME MATRICES INTO CSV FILES
    std::ofstream file("solution.csv");
    file << "X_matrix" << std::endl;
    solution_csv.save(file, arma::csv_ascii);
    file.close();
    std::ofstream file2("MaxI.csv");
    file2 << "Max_I" << std::endl;
    Max_I.save(file2, arma::csv_ascii);
    file2.close();
    std::ofstream file3("time.csv");
    file3 << "time" << std::endl;
    time.save(file3, arma::csv_ascii);
    file3.close();

    arma::mat h_h = h_history;
    std::ofstream file4("h_history.csv");
    file4 << "timestep_h" << std::endl;
    h_h.save(file4, arma::csv_ascii);
    file4.close();

    arma::mat h_1_h = h_1_history;
    std::ofstream file5("h_1_history.csv");
    file5 << "timestep_h_1" << std::endl;
    h_1_h.save(file5, arma::csv_ascii);
    file5.close();

    std::ofstream file6("h_step.csv");
    file6 << "timestep_h_1" << std::endl;
    h_1_h.save(file6, arma::csv_ascii);
    file6.close();

    auto t2 = std::chrono::high_resolution_clock::now(); // End time
    std::chrono::duration<double, std::milli> time_span = (t2 - t1);
    std::cout << "Total time:" << time_span.count() << "ms\n";

    return 0;
    /*-----------------------------------------------------------------------------------------------------------*/
}
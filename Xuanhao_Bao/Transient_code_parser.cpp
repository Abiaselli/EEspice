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



// TRANSIENT SIMULATION SETTINGS
// double t_start = 0;
// double t_end = 2e-8;
// double h = t_end/5000; // t_end/5000 is the default value

/*  TOTAL NUMBER OF NODES EXCLUDING GROUND
    Two port components such as resistors, initially adds 2 nodes. If more than 1 component is added, it then adds 1 node per component.
    Number of MOSFETs and cascaded levels are assigned above too. External nodes are the nodes which  */



#include "Transient_code_parser.hpp"

/*  The line above the code means that the section can be changed by the user which analyses circuit simulation

    --fixed--           : can't be edited for circuit simulation purposes, but can be edited for code debugging.
    --can be changed--  : can be edited for circuit simulation purposes and edited for code debugging.
    
*/


// Main function for the circuit simulation
int main(int argc, const char ** argv)
{
    auto t1 = std::chrono::high_resolution_clock::now(); // Start time
    CKTcircuit ckt;
    DenseMatrix dematrix;

    CircuitParser parser("single_RC.cir");
    parser.parser();

    CKTsetup(ckt, parser, dematrix);
    ckt.setcktmatrix(dematrix);

    CKTload(ckt);
    std::cout << "h = " << ckt.h << std::endl;

    /*--------------------------------------------can be changed-------------------------------------------------*/
    // std::vector<double> RHS_locate = {
    //     // Assigning the voltage matrix on LHS and RHS for the pulse voltage
    //     Vs_assigner(1,0,0,LHS,RHS)
    // };
    // Assigning DC voltage sources
    // Vs_assigner(supply_voltage_node,0,5,LHS,RHS); // supply voltage for the ring oscillator, vdd

    // Assigning the stamps that would affect the RHS in transient simulation 
    // (only for  time-dependent voltage, e.g. pulse voltages)
    
    /*----------------------------------------------fixed--------------------------------------------------------*/
    // Checking the LHS and RHS matrices
    dematrix.LHS.print("LHS matrix =");
    dematrix.RHS.print("RHS matrix =");
    // ckt.cktdematrix->LHS.print("cktLHS matrix =");
    // ckt.cktdematrix->RHS.print("cktRHS matrix =");
    /*----------------------------------------------fixed--------------------------------------------------------*/
    // OPERATING POINT ANALYSIS SYSTEM
    ckt.mode = 0; // 0 to do OP analysis, 1 to do transient simulation

    // The initial LHS and RHS values to be used in the NR-algorithm
    arma::mat const init_LHS = dematrix.LHS;
    arma::mat const init_RHS = dematrix.RHS;
    
    // OP analysis used as initial condition for next evaluation
    arma::vec solution= arma::zeros(dematrix.RHS.n_rows,dematrix.RHS.n_cols);
    // Benchmarking for OP analysis
    auto tstart_op = std::chrono::high_resolution_clock::now();
    solution = NewtonRaphson_system(init_LHS,init_RHS,solution,ckt);
    auto tstop_op = std::chrono::high_resolution_clock::now();

    solution.print("The OP analysis of the circuit is: ");
    /*----------------------------------------------fixed--------------------------------------------------------*/
    // The solution csv that is going to be plotted which contains the values of nodal voltages
    // and voltage source currents
    arma::vec solution_csv = solution;
    arma::vec Max_I = arma::zeros(1,1);
    Max_I.row(0).col(0) = dematrix.RHS.n_rows;
    arma::mat zero_ext = arma::zeros(1,1);
    /*--------------------------------------------can be changed-------------------------------------------------*/
    // ADDING TRANSIENT SIMULATION LOOP (includes V_pulse or any time dependent sources)
    int i = 0;
    double time_trans = ckt.t_start;
    ckt.mode = 1;
    auto tstart_trans = std::chrono::high_resolution_clock::now();
    while(time_trans < ckt.t_end){
        
        dematrix.LHS = init_LHS;
        dematrix.RHS = init_RHS;

        if(ckt.pulse_num > 0){
            std::vector<double> RHS_value;
            dematrix.RHS = RHS_update(ckt,init_RHS, RHS_value);
        }

        // std::vector<double> RHS_value = {
        //     V_pulse(V1,V2,t1,td,tr,tf,tpw,tper,h)
        // };
        // RHS = RHS_update(RHS_locate, init_RHS, RHS_value);

        // Calling the Newton-Raphson system here
        solution = NewtonRaphson_system(init_LHS,init_RHS,solution,ckt);
        // Assigning the variables that will be plotted and analysed as seen in a circuit simulator
        solution_csv = arma::join_cols(solution_csv,solution);
        i++;
        time_trans += ckt.h;
    }

    auto tstop_trans = std::chrono::high_resolution_clock::now();
    // time vector to be inputted in plot for python analysis
    arma::mat time = arange(ckt.t_start,ckt.h,i);
    
    /* Getting number of milliseconds as a double. */
    std::chrono::duration<double, std::milli> OP_time = (tstop_op - tstart_op) ;
    std::chrono::duration<double, std::milli> trans_time = (tstop_trans - tstart_trans);

    std::cout << "DC OP time:" <<  OP_time.count() << "ms\n";
    std::cout << "Transient time:" <<  trans_time.count() << "ms\n";

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

    auto t2 = std::chrono::high_resolution_clock::now(); // End time
    std::chrono::duration<double, std::milli> time_span = (t2 - t1) ;
    std::cout << "Total time:" <<  time_span.count() << "ms\n";

    return 0;
    /*-----------------------------------------------------------------------------------------------------------*/
}


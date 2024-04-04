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

// Global variables
// std::vector<Transient> history_trans;
// std::deque<arma::vec> history_voltages;
// std::deque<double> history_steps;


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
    Transient trans_op;

    CircuitParser parser("Inverter.cir");
    parser.parser();

    CKTsetup(ckt, parser, dematrix);            // Pass the parser to the ckt and the initialise LHS and RHS matrices
    ckt.setcktmatrix(dematrix);

    CKTload(ckt);
    ckt.cktdematrix -> set_initmatrix();        // Set the initial LHS and RHS matrices
    Transsetup(trans_op, parser, ckt);
    trans_op.LHS = dematrix.get_init_LHS();
    trans_op.RHS = dematrix.get_init_RHS();
    trans_op.C_list = ckt.C_list;               // Pass the capacitance list to the transient analysis
    trans_op.h = 0;
    // history_steps.push_front(trans_op.h);
    
    

    /*--------------------------------------------can be changed-------------------------------------------------*/


    /*----------------------------------------------fixed--------------------------------------------------------*/
    // Checking the LHS and RHS matrices
    dematrix.LHS.print("LHS matrix =");
    dematrix.RHS.print("RHS matrix =");
    /*----------------------------------------------fixed--------------------------------------------------------*/
    // OPERATING POINT ANALYSIS SYSTEM
    trans_op.mode = 0; // 0 to do OP analysis, 1 to do transient simulation

    // OP analysis used as initial condition for next evaluation
    arma::vec solution= arma::zeros(dematrix.RHS.n_rows,dematrix.RHS.n_cols);
    
    // Benchmarking for OP analysis
    auto tstart_op = std::chrono::high_resolution_clock::now();

    auto matrixs_op = DynamicNonLinear(ckt, 0, solution, 0, 0, trans_op.C_list);  // mode = 0 for DC OP analysis
    solution = NewtonRaphson_system(matrixs_op.first, matrixs_op.second, ckt, solution);
    if(trans_op.C_list.size() > 0){
        auto cur_vol_op = get_currents_voltages(trans_op.C_list, 1e-25, solution, arma::zeros(dematrix.RHS.n_rows,dematrix.RHS.n_cols));
        trans_op.Capacitance = get_capacitance(trans_op.C_list);                    // Get the capacitance matrix from c_list
        trans_op.C_current = cur_vol_op.first;
        trans_op.C_voltage = cur_vol_op.second;
        trans_op.C_charge = trans_op.C_voltage % trans_op.Capacitance;

    }
   
    auto tstop_op = std::chrono::high_resolution_clock::now();

    trans_op.solution = solution;
    history_trans_update(trans_op);
    
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


    auto tstart_trans = std::chrono::high_resolution_clock::now();
    std::cout << "transient simulation start" << std::endl;

    std::deque<double> breakpoints;

    if(ckt.pulse_num != 0){
        breakpoints = get_breakpoints(ckt, trans_op);  // Get the time of breakpoints for the transient simulation
        std::cout << "The breakpoints are: ";
        for(double num : breakpoints) {
            std::cout << num << " ";
        }
        std::cout << std::endl;

    }
    
    
    do{
        std::cout << "-------------------------------------------------------" << std::endl;
        std::cout << "Next step" << std::endl;

        Transient trans;
        trans.t_end = trans_op.t_end;
        trans.h_MAX = trans_op.h_MAX;
        trans.h_MIN = trans_op.h_MIN;
        trans.mode = 1; // 0 to do OP analysis, 1 to do transient simulation
        trans.C_list = trans_op.C_list;

        // Fixed time step
        if(ckt.pulse_num == 0 || trans.C_list.size() == 0){
            trans.h = trans_op.init_h;  // Initial time step
            trans.time_trans = vec_trans.back().time_trans;
            trans.C_list = vec_trans.back().C_list;

            // If circuit includes V_pulse.
            if(trans.time_trans + trans.h > breakpoints.front() && ckt.pulse_num > 0){
                std::cout<< breakpoints.front()<< std::endl;
                Transient breakpoints_trans;
                breakpoints_trans.t_end = trans_op.t_end;
                breakpoints_trans.h_MAX = trans_op.h_MAX;
                breakpoints_trans.h_MIN = trans_op.h_MIN;
                breakpoints_trans.mode = 1; // 0 to do OP analysis, 1 to do transient simulation
                breakpoints_trans.C_list = vec_trans.back().C_list;
            

                std::cout << "The time step is too large, back to breakpoint" << std::endl;
                breakpoints_trans.time_trans = breakpoints.front();
                breakpoints_trans.h = breakpoints_trans.time_trans - vec_trans.back().time_trans;
                auto matrixs = DynamicNonLinear(ckt, breakpoints_trans.h, vec_trans.back().solution, 1, breakpoints_trans.time_trans, breakpoints_trans.C_list);
                solution = NewtonRaphson_system(matrixs.first, matrixs.second, ckt, vec_trans.back().solution);

                breakpoints_trans.solution = solution;
                breakpoints_trans.LHS = matrixs.first;
                breakpoints_trans.RHS = matrixs.second;
                breakpoints_trans.trans_count = vec_trans.back().trans_count + 1;

                history_steps.push_front(breakpoints_trans.h);
                history_trans_update(breakpoints_trans);

                solution.print("The transient analysis of the circuit is: ");
                std::cout << "time step: " << breakpoints_trans.h << std::endl;
                std::cout << "time_trans: " << breakpoints_trans.time_trans << std::endl;

                solution_csv = arma::join_cols(solution_csv, solution);
                breakpoints.pop_front();
                continue;

            }

            auto matrixs = DynamicNonLinear(ckt, trans.h, vec_trans.back().solution, trans.mode, trans.time_trans, trans.C_list);  // Update LHS and RHS matrices
            solution = NewtonRaphson_system(matrixs.first, matrixs.second, ckt, vec_trans.back().solution);
            solution.print("The transient analysis of the circuit is: ");

            if(trans.C_list.size() > 0){

                auto cur_vol = get_currents_voltages(trans.C_list, trans.h, solution, vec_trans.back().solution);
                trans.Capacitance = get_capacitance(trans.C_list);
                
                
                trans.C_current = cur_vol.first;
                trans.C_voltage = cur_vol.second;
                trans.C_charge = trans.C_voltage % trans.Capacitance;
                trans.C_list_update();
            }
            
            trans.time_trans = vec_trans.back().time_trans + trans.h;
            trans.solution = solution;
            trans.LHS = matrixs.first;
            trans.RHS = matrixs.second;
            trans.trans_count = vec_trans.back().trans_count + 1;
            // trans.C_current.print("The current matrix is: ");
            std::cout << "time trans: " << trans.time_trans << std::endl;
            std::cout << "time step: " << trans.h << std::endl;

            history_steps.push_front(trans.h);
            history_trans_update(trans);
            // solution_csv.print("The solution matrix is: ");
            solution_csv = arma::join_cols(solution_csv, solution);

            continue;

        }



        if(vec_trans.size() == 1){
            // The first step of the transient simulation
            trans.h = trans_op.init_h;  // Initial time step
            trans.time_trans = trans.t_start + trans.h;

            auto matrixs = DynamicNonLinear(ckt, trans.h, trans_op.solution, trans.mode, trans.time_trans, trans.C_list);  // Update LHS and RHS matrices
            solution = NewtonRaphson_system(matrixs.first, matrixs.second, ckt, trans_op.solution);
            solution.print("The transient analysis of the circuit is: ");

            if(trans.C_list.size() > 0){

                auto cur_vol = get_currents_voltages(trans.C_list, trans.h, solution, trans_op.solution);
                trans.Capacitance = get_capacitance(trans.C_list);
                
                trans.C_current = cur_vol.first;
                trans.C_voltage = cur_vol.second;
                trans.C_charge = trans.C_voltage % trans.Capacitance;
                trans.C_list_update();
            }
            
            trans.solution = solution;
            trans.LHS = matrixs.first;
            trans.RHS = matrixs.second;
            trans.trans_count += 1;
            // trans.C_current.print("The current matrix is: ");

            history_steps.push_front(trans.h);
            history_trans_update(trans);
            // solution_csv.print("The solution matrix is: ");
            solution_csv = arma::join_cols(solution_csv, solution);

            // solution_csv.print("The solution matrix is: ");
        }

       

        else{

            trans.C_list = vec_trans.back().C_list;

            trans.time_trans = vec_trans.back().time_trans;

            solution = multi_next_h(trans, ckt);
            
            trans.time_trans = trans.time_trans + trans.h;


            if(trans.time_trans > breakpoints.front()){
                
                Transient breakpoints_trans;
                breakpoints_trans.t_end = trans_op.t_end;
                breakpoints_trans.h_MAX = trans_op.h_MAX;
                breakpoints_trans.h_MIN = trans_op.h_MIN;
                breakpoints_trans.mode = 1; // 0 to do OP analysis, 1 to do transient simulation
                breakpoints_trans.C_list = vec_trans.back().C_list;
            

                std::cout << "The time step is too large, back to breakpoint" << std::endl;
                breakpoints_trans.time_trans = breakpoints.front();
                breakpoints_trans.h = breakpoints_trans.time_trans - vec_trans.back().time_trans;
                auto matrixs = DynamicNonLinear(ckt, breakpoints_trans.h, vec_trans.back().solution, 1, breakpoints_trans.time_trans, breakpoints_trans.C_list);
                solution = NewtonRaphson_system(matrixs.first, matrixs.second, ckt, vec_trans.back().solution);
                auto cur_vol = get_currents_voltages(breakpoints_trans.C_list, breakpoints_trans.h, solution, vec_trans.back().solution);
                breakpoints_trans.Capacitance = get_capacitance(breakpoints_trans.C_list);

                breakpoints_trans.solution = solution;
                breakpoints_trans.C_current = cur_vol.first;
                breakpoints_trans.C_voltage = cur_vol.second;
                breakpoints_trans.C_charge = breakpoints_trans.C_voltage % breakpoints_trans.Capacitance;
                breakpoints_trans.C_list_update();

                breakpoints_trans.LHS = matrixs.first;
                breakpoints_trans.RHS = matrixs.second;
                breakpoints_trans.trans_count = vec_trans.back().trans_count + 1;

                history_steps.push_front(breakpoints_trans.h);
                history_trans_update(breakpoints_trans);

                solution.print("The transient analysis of the circuit is: ");
                std::cout << "time step: " << breakpoints_trans.h << std::endl;
                std::cout << "time_trans: " << breakpoints_trans.time_trans << std::endl;

                solution_csv = arma::join_cols(solution_csv, solution);
                breakpoints.pop_front();
                continue;

            }

            else{
                solution.print("The transient analysis of the circuit is: ");
                std::cout <<"the time step is: "<< trans.h << std::endl;
                std::cout << "time_trans: " << trans.time_trans << std::endl;


                trans.C_list_update();
                trans.trans_count = vec_trans.back().trans_count + 1;
                history_trans_update(trans);
                history_steps.push_front(trans.h);
                // solution_csv.print("The solution matrix is: ");
                solution_csv = arma::join_cols(solution_csv, solution);
                // solution_csv.print("The solution matrix is: ");

            }
            

        }

        // Update the data
        // solution_csv = arma::join_cols(solution_csv, solution);

    }while (vec_trans.back().time_trans < vec_trans.back().t_end);

    auto tstop_trans = std::chrono::high_resolution_clock::now();
    // time vector to be inputted in plot for python analysis
    arma::mat time = arange(trans_op.t_start, history_steps);

    
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

    save_csv();

    return 0;
    /*-----------------------------------------------------------------------------------------------------------*/
}


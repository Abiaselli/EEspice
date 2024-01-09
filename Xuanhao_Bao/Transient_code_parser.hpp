#pragma once

#include <iostream>
#include "armadillo"
#include <fstream>
#include <tuple>
#include <cmath>
#include <chrono>
#include <sstream>
#include <string>
#include <variant>
#include <vector>
#include <algorithm>

// Forward declarations
class CircuitElement;
void R_assigner(double node_x, double node_y, double R, arma::mat &LHS, arma::mat &RHS);
void Vs_assigner(int node_x, int node_y, double V_value, arma::mat &LHS, arma::mat &RHS);
double convertToValue(const std::string& valueStr);
arma::mat branch_ext(arma::mat M, int node_x, int node_y);
void Is_assigner(double node_x, double node_y, double I, arma::mat &LHS, arma::mat &RHS);
void C_assigner(int node_x,int node_y,double C, double h, arma::mat &LHS, arma::mat &RHS, arma::mat solution, int mode);
double V_pulse_assigner(int node_x, int node_y, double V_value, arma::mat &LHS, arma::mat &RHS);
double V_pulse_value(double V1, double V2, double &t1,double td, double tr, double tf, double tpw, double tper, double h);
arma::vec arange(double tstart, double h, double vec_size);
/*-----------------------------------------------------------------------------------------------------------------------------------------------------*/


double cond(double R){
    return 1/R;
}

// extra -  adding arange() function here to return the time vector
arma::vec arange(double tstart, double h, double vec_size){
    arma::vec time = arma::zeros(vec_size,1);
    for(int i=0; i<vec_size; i++)
    {
        time[i] = tstart;
        tstart = tstart + h;
    }
    return time;
}

class VoltageSource {
public:
    int id;
    int nodePos, nodeNeg;
    double value;
};

class Pulsevoltage{
public:
    int id;
    int nodePos, nodeNeg;
    double t1_pulse;
    double V1;
    double V2;
    double td;
    double tr;
    double tf;
    double pw;
    double per;
};

class CurrentSource {
public:
    int id;
    int nodePos, nodeNeg;
    double value;
};
 
class Resistor {
public:
    int id;
    int nodePos, nodeNeg;
    double value;
};
 
class Capacitor {
public:
    int id;
    int nodePos, nodeNeg;
    double value;
};

class CircuitElement {
public:
  std::variant<VoltageSource, CurrentSource, Resistor, Capacitor, Pulsevoltage> element;
};

// Parser class

class CircuitParser {
public:
    std::string filename;
    std::vector<CircuitElement> elements;

    CircuitParser(const std::string& filename) : filename(filename) {}
 
    void parser() {

        std::ifstream file(filename);

        if (!file.is_open()) {

            std::cerr << "Error opening file: " << filename << std::endl;

            return;

        }

        std::string line;

        std::getline(file, line);  // Read and discard the first line

        while (std::getline(file, line)) {

            size_t commentPos = line.find('*');

            if (commentPos != std::string::npos) {

                line = line.substr(0, commentPos);  // Remove comment

            }
 
            if (!line.empty()) {

                parseLine(line);

            }

        }

    }
 
    const std::vector<CircuitElement>& getCircuitElements() const {

        return elements;

    }

    int getMaxNode() const {
        int maxNode = 0;
        for (const auto& element : elements) {
            std::visit([&maxNode](auto&& arg) {
                maxNode = std::max(maxNode, arg.nodePos);
                maxNode = std::max(maxNode, arg.nodeNeg);
            }, element.element);
        }
        return maxNode;
    }


    double double_t_end;  // This double_t_end can be passed to the CKTcircuit class
 
    void parseLine(const std::string& line) {

        std::istringstream iss(line);
        std::string type, valueStr;
        int v_nodePos, v_nodeNeg;
        std::string v_type, t1_pulse, v1, v2, td, tr, tf, pw, per; 
        
        

        iss >> type;  // Automatically skips leading whitespace before reading type
        std::cout<<"type is "<<type<<std::endl;

        if (type.empty()) {
            return;  // Skip empty lines or lines with only whitespaces
        }

 
        if (type[0] == 'V') {
            int v_id = std::stoi(type.substr(1));

            iss >> v_nodePos >> v_nodeNeg >> v_type;
            if (v_type == "pulse"){   // Pulse voltage settings
                Pulsevoltage pv;

                pv.id = v_id;
                pv.nodePos = v_nodePos;
                pv.nodeNeg = v_nodeNeg;

               iss >> t1_pulse >> v1 >> v2 >> td >> tr >> tf >> pw >> per;

                pv.t1_pulse = convertToValue(t1_pulse);
                pv.V1 = convertToValue(v1);
                pv.V2 = convertToValue(v2);
                pv.td = convertToValue(td);
                pv.tr = convertToValue(tr);
                pv.tf = convertToValue(tf);
                pv.pw = convertToValue(pw);
                pv.per = convertToValue(per);

                elements.push_back(CircuitElement{pv});
            }
            else{

                VoltageSource vs;

                vs.id = v_id;
                vs.nodePos = v_nodePos;
                vs.nodeNeg = v_nodeNeg;
               
                vs.value = convertToValue(v_type);

                elements.push_back(CircuitElement{vs});
            }

        } else if (type[0] == 'R') {

            Resistor r;

            r.id = std::stoi(type.substr(1));

            iss >> r.nodePos >> r.nodeNeg >> valueStr;

            r.value = convertToValue(valueStr);

            elements.push_back(CircuitElement{r});

        } else if (type[0] == 'C') {

            Capacitor c;

            c.id = std::stoi(type.substr(1));

            iss >> c.nodePos >> c.nodeNeg >> valueStr;

            c.value = convertToValue(valueStr);

            elements.push_back(CircuitElement{c});

        }else if (type[0] == 'I') {

            CurrentSource cs;    

            cs.id = std::stoi(type.substr(1));

            iss >> cs.nodePos >> cs.nodeNeg >> valueStr;

            cs.value = convertToValue(valueStr);

            elements.push_back(CircuitElement{cs});

        }else if (type == ".tran") {
            std::string ignord, string_t_end;

            iss >> ignord >> string_t_end;

            double_t_end = convertToValue(string_t_end);
        
        }else {
                
            std::cerr << "Error: Unknown element type: " << type << std::endl;
        }

    }

};

double convertToValue(const std::string& valueStr) {

    size_t unitPos = valueStr.find_first_not_of("0123456789.-eE");  // Find the position of the first non-numeric character
    double value = 0;

    try{
        value = std::stod(valueStr.substr(0, unitPos));  // Convert the numeric part of the string to double
    }
    catch(const std::invalid_argument& ia){
        std::cerr << "Error: Invalid argument: " << ia.what() << std::endl;
        return 0;
    }
    catch(const std::out_of_range& oor){
        std::cerr << "Error: Out of Range error: " << oor.what() << std::endl;
        return 0;
    }

    if(unitPos != std::string::npos){
        char unit = valueStr[unitPos];          // Get the unit character
        switch (unit) {

            case 'k': return value * 1000;      // Kilo 
            case 'u': return value * 1e-6;      // Micro          
            case 'n': return value * 1e-9;      // Nano
            case 'p': return value * 1e-12;     // Pico
            case 'f': return value * 1e-15;     // Femto
            default: std::cerr << "Error: Unknown unit: " << unit << std::endl;
        }  
    }

    return value; // No unit or unrecognized unit, assume the value is in base units

}

class DenseMatrix{
public:
    int Maxi{};                             // Size of matrix
    int Maxj{};

    arma::mat LHS;                          // LHS matrix
    arma::mat RHS;                          // RHS matrix
    // arma::mat init_LHS;                     // The initial LHS and RHS values to be used in the NR-algorithm
    // arma::mat init_RHS;

    // void set_initmatrix(){
    //     init_LHS = LHS;
    //     init_RHS = RHS;
    // }
};

class CKTcircuit{
public:

    int const code = 0;                         // Choosing code for transient simulation
    int const supply_voltage_node = 1;          // supply voltage node for the ring oscillator
    double t_start = 0;                         // TRANSIENT SIMULATION SETTINGS
    double t_end;
    double h{};                                 // time step
    int mode{};                                   // 0 to do OP analysis, 1 to do transient simulation
    int pulse_num{};                            // Number of pulse voltages
    std::vector<double> RHS_locate;                      // Locations of RHS need to be changed in transient simulation. 
    // uword is a typedef for an unsigned integer type; it is used for matrix indices as well as all internal counters and loops                                           
    
    std::vector<CircuitElement> CKTelements;    // Vector of circuit elements
    int external_nodes{};                       // Number of external nodes (excluding ground and ring oscillator loop nodes)
    //int external_mosfets{};                   // Number of standalone mosfets (excluding mosfets from ring oscillator)
    int no_of_mosfets{};                        // Total number of MOSFETs
    int T_nodes{};                              // Total number of nodes excluding ground

    DenseMatrix* cktdematrix;                   // Dense matrix class object
    void setcktmatrix(DenseMatrix& DenseMatrix){
        cktdematrix = &DenseMatrix;
    }
};


void CKTsetup(CKTcircuit &ckt, const CircuitParser &parser, DenseMatrix &DenseMatrix){
    ckt.t_end = parser.double_t_end;
    ckt.h = ckt.t_end/5000;                   // t_end/5000 is the default value
    ckt.mode = 0;                             // 0 to do OP analysis, 1 to do transient simulation   
    
    // Careful! getCircuitElements function is const, so it can't be used to modify the elements vector
    //ckt.elements = parser.getCircuitElements(); 
    ckt.CKTelements = parser.elements;  
    ckt.external_nodes = parser.getMaxNode();
    ckt.no_of_mosfets = 0;
    ckt.T_nodes = ckt.external_nodes + 4*ckt.no_of_mosfets;

    // Size of matrix
    DenseMatrix.Maxi = ckt.T_nodes;
    DenseMatrix.Maxj = DenseMatrix.Maxi;
    DenseMatrix.LHS = arma::zeros(DenseMatrix.Maxi,DenseMatrix.Maxj);   // LHS matrix
    DenseMatrix.RHS = arma::zeros(DenseMatrix.Maxi,1);                  // RHS matrix
}

void CKTload(CKTcircuit &ckt){
    // ASSIGNING THE STAMPS TO THE LHS AND RHS MATRICES

    for (auto& element : ckt.CKTelements) {
        std::visit([&](auto&& arg) {
            if constexpr(std::is_same_v<std::decay_t<decltype(arg)>, Resistor>){
                std::cout << "R Element ID: " << arg.id << ", Node Pos: " << arg.nodePos << ", Node Neg: " << arg.nodeNeg << "value: "<< arg.value << std::endl;

                R_assigner(arg.nodePos, arg.nodeNeg, arg.value, ckt.cktdematrix->LHS, ckt.cktdematrix->RHS);

            }else if constexpr(std::is_same_v<std::decay_t<decltype(arg)>, VoltageSource>){
                std::cout << "VS Element ID: " << arg.id << ", Node Pos: " << arg.nodePos << ", Node Neg: " << arg.nodeNeg << "value: "<< arg.value << std::endl;

                Vs_assigner(arg.nodePos, arg.nodeNeg, arg.value, ckt.cktdematrix->LHS, ckt.cktdematrix->RHS);

            }else if constexpr(std::is_same_v<std::decay_t<decltype(arg)>, Pulsevoltage>){
                std::cout << "VPulse Element ID: " << arg.id << ", Node Pos: " << arg.nodePos << ", Node Neg: " << arg.nodeNeg << std::endl;

                ckt.pulse_num++;

                ckt.RHS_locate.push_back(V_pulse_assigner(arg.nodePos, arg.nodeNeg, arg.V1, ckt.cktdematrix->LHS, ckt.cktdematrix->RHS));

            }else if constexpr(std::is_same_v<std::decay_t<decltype(arg)>, CurrentSource>){
                std::cout << "I Element ID: " << arg.id << ", Node Pos: " << arg.nodePos << ", Node Neg: " << arg.nodeNeg << "value: "<< arg.value << std::endl;

                Is_assigner(arg.nodePos, arg.nodeNeg, arg.value, ckt.cktdematrix->LHS, ckt.cktdematrix->RHS);

            }
            
        }, element.element);
    }
}

// create resistor matrix stamp
void R_assigner(double node_x, double node_y, double R, arma::mat &LHS, arma::mat &RHS){
    int maxi = LHS.n_cols;
    int maxj = LHS.n_rows;
    double x = 0;
    arma::mat a = arma::zeros(maxi,maxj);

    if(R == 0)
        x = 0;
    else
        x = cond(R);

    if((node_x == 0) && (node_y == 0)){
        x = 0;
    }
    else{
        if(node_x == 0){
            a.row(node_y-1).col(node_y-1) = x;
        }
        else if(node_y == 0){
            a.row(node_x-1).col(node_x-1) = x;
        }
        else{
            a.row(node_x-1).col(node_x-1) = x;
            a.row(node_x-1).col(node_y-1) = -x;
            a.row(node_y-1).col(node_x-1) = -x;
            a.row(node_y-1).col(node_y-1) = x;
        }
    }
    LHS = LHS + a;
}

// Voltage source stamp assigner
void Vs_assigner(int node_x, int node_y, double V_value, arma::mat &LHS, arma::mat &RHS){

    arma::vec value(1);
    value = V_value;
    // Extending the branch at the LHS matrix
    LHS = branch_ext(LHS, node_x, node_y);

    // Assigning the value at RHS
    RHS = arma::join_cols(RHS, value);

    // location of voltage value for transient simulation
    // double V_locate1 = RHS.n_rows;

    // return V_locate1;
}

// branch extender function
arma::mat branch_ext(arma::mat M, int node_x, int node_y){
    int M_cols = M.n_cols;
    int M_rows = M.n_rows;
    arma::mat va = arma::zeros(M_rows,1);
    if(node_x == 0){
        va.row(node_y-1).col(0) = 1;
    }
    else if(node_y == 0){
        va.row(node_x-1).col(0) = 1;
    }
    else{
        va.row(node_x-1).col(0) = -1;
        va.row(node_y-1).col(0) = 1;
    }

    arma::mat zero_ext = arma::zeros(1,1);
    arma::mat ha = va.as_row();
    arma::mat haz = arma::join_rows(ha,zero_ext);

    arma::mat M1 = arma::join_rows(M,va);
    arma::mat M2 = arma::join_cols(M1,haz);

    return M2;
}

// create a current matrix stamp
void Is_assigner(double node_x, double node_y, double I, arma::mat &LHS, arma::mat &RHS){
    int maxi = LHS.n_cols;
    int maxj = 1;
    arma::mat a = arma::zeros(maxi,maxj);
    if((node_x == 0) && (node_y == 0)){
        I = 0;
    }
    else{
        if(node_x == 0){
            a.row(node_y-1).col(0) = I;
        }
        else if(node_y == 0){
            a.row(node_x-1).col(0) = -I;
        }
        else{
            a.row(node_x-1).col(0) = -I;
            a.row(node_y-1).col(0) = I;
        }
    }
    RHS = RHS + a;
}

// Capacitor stamp assigner
void C_assigner(int node_x,int node_y,double C, double h, arma::mat &LHS, arma::mat &RHS, arma::mat solution, int mode){
    
    double x = 0;
    double x1 = 0;
    // this if else statment uses trapezoidal formula
    if(mode>0){
        x = C/h;
        if(node_x == 0){
            x1 = C*(solution(node_y-1,0))/h;
        }else if(node_y == 0){
            x1 = C*(solution(node_x-1,0))/h;
        }else{
            x1 = C*(solution(node_x-1,0)-solution(node_y-1,0))/h;
        }
    }
    else{
        x = 0;
        x1 = 0;
    }
    // Matrix stamp for a capacitor on LHS
    R_assigner(node_x,node_y,cond(x),LHS,RHS);
    //std::cout<<"cond x value is "<<cond(x)<<std::endl;
    // Matrix stamp for a capacitor on RHS
    //std::cout<<"x1 value is "<<x1<<std::endl;
    Is_assigner(node_x,node_y,-x1,LHS,RHS);
}

double V_pulse_assigner(int node_x, int node_y, double V_value, arma::mat &LHS, arma::mat &RHS){

    arma::vec value(1);
    value = V_value;
    // Extending the branch at the LHS matrix
    LHS = branch_ext(LHS, node_x, node_y);

    // Assigning the value at RHS
    RHS = arma::join_cols(RHS, value);

    // location of voltage value for transient simulation
    double Pulse_RHS_locate1 = RHS.n_rows;

    return Pulse_RHS_locate1;
}

double V_pulse_value(double V1, double V2, double &t1,double td, double tr, double tf, double tpw, double tper, double h){
    double V_t1 = V1;
    if((t1 >= 0) && (t1 < td)){
        V_t1 = V1;
    }else if((t1 >= td) && (t1 < (td + tr))){
        V_t1 = V1 + (V2-V1)*(t1-td)/tr;
    }else if((t1 >= (td + tr)) && (t1 < (td + tr + tpw))){
        V_t1 = V2;
    }else if((t1 >= (td + tr + tpw)) && (t1 < (td + tr + tpw + tf))){
        V_t1 = V2 + (V1 - V2)*(t1-(td+tr+tpw))/tf;
    }else if((t1 >= (td + tr + tpw + tf)) && (t1 <= (td + tper))){
        V_t1 = V1;
    }else{
        t1 = td;
    }
    t1 = t1 + h;
    return V_t1;
}

// Updates the value of RHS
arma::mat RHS_update(CKTcircuit &ckt, arma::mat RHS, std::vector<double> &val){
    for (auto& element : ckt.CKTelements) {
        std::visit([&](auto&& arg) {
            if constexpr(std::is_same_v<std::decay_t<decltype(arg)>, Pulsevoltage>){
                val.push_back(V_pulse_value(arg.V1, arg.V2, arg.t1_pulse, arg.td, arg.tr, arg.tf, arg.pw, arg.per, ckt.h));
                //std::cout<<"V_pulse_value is called"<<std::endl;
            }
        }, element.element);
    }
   
    for(int i = 0; i < ckt.RHS_locate.size(); i++){
        RHS.row(ckt.RHS_locate[i]-1).col(0) = val[i];
    }

    return RHS;
}

// Assigning the stamp matrices for dynamic and non-linear components
std::pair<arma::mat,arma::mat> DynamicNonLinear( arma::mat solution, CKTcircuit &ckt){
    // (All the circuit assigners can be used except for voltage and current sources)
    /*--------------------------------------------can be changed-------------------------------------------------*/
    // RingOscillatorStages(W_oscillator, L_oscillator, R_oscillator, C_oscillator, LHS, RHS, solution, h, mode);
    for (auto& element : ckt.CKTelements) {
        std::visit([&](auto&& arg) {
            if constexpr(std::is_same_v<std::decay_t<decltype(arg)>, Capacitor>){
                //std::cout << "C Element ID: " << arg.id << ", Node Pos: " << arg.nodePos << ", Node Neg: " << arg.nodeNeg << "value: "<< arg.value << std::endl;
                C_assigner(arg.nodePos, arg.nodeNeg, arg.value, ckt.h, ckt.cktdematrix->LHS, ckt.cktdematrix->RHS, solution, ckt.mode);
            }
            
        }, element.element);
    }

    arma::mat J_x = ckt.cktdematrix->LHS;
    arma::mat F_x = ckt.cktdematrix->RHS;
    
    return {J_x,F_x};
}

// Newton Raphson system solver for non-linear and dynamic elements
arma::mat NewtonRaphson_system(arma::mat const init_LHS, arma::mat const init_RHS, arma::mat &solution, CKTcircuit &ckt){
    int col_size = ckt.cktdematrix->LHS.n_cols;
    int row_size = ckt.cktdematrix->LHS.n_rows;
    double eps_val = 1e-8;
    arma::mat error = arma::zeros(1,1);
    double error_val = 9e9;
    error.row(0) = error_val;
    int iteration_counter = 0;
    arma::mat delta = arma::zeros(row_size,1);
    while((error(0,0) > eps_val) && (iteration_counter < 50)){ // iteration counter can be changed depending on the non-linearity of the circuit
        auto matrices = DynamicNonLinear(solution, ckt);
        delta = arma::solve(matrices.first,(matrices.first*solution) - matrices.second);
        error.row(0) = arma::max(arma::abs(delta));
        solution -= delta;
        iteration_counter += 1;
    }
    return solution;
}
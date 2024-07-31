#pragma once
#define ARMA_DONT_USE_WRAPPER
#define ARMA_USE_MKL_ALLOC
#define ARMA_USE_SUPERLU

#include <iostream>
#include <armadillo>
#include <fstream>
#include <tuple>
#include <cmath>
#include <chrono>
#include <sstream>
#include <string>
#include <variant>
#include <vector>
#include <algorithm>
#include <deque>
#include <iomanip>
#include <typeinfo>
#include "Global_variables.hpp"
#include "BS_thread_pool.hpp"
#include "BS_thread_pool_utils.hpp"
#include "XB_timer.hpp"

// Macro to enable debug code
#define DEBUG_PRINT(x) if (isDebugMode()) { std::cout << x << std::endl; }
#define ARMA_PRINT(var, label) if (isDebugMode()) { (var).print(label); }


// Forward declarations
class XB_Timer;
class CircuitElement;
class Transient;
class CKTcircuit;
class multi_timestep;
class Truncation_error;
class Capacitor;

void setDebugMode(bool mode);
bool isDebugMode();
void R_assigner(int node_x, int node_y, double G, arma::sp_mat &LHS);
void Vs_assigner(int node_x, int node_y, double V_value, arma::sp_mat &LHS, arma::mat &RHS);
double convertToValue(const std::string& valueStr);
// arma::mat branch_ext(const arma::mat &M, int node_x, int node_y);
arma::sp_mat branch_ext(const arma::sp_mat &M, int node_x, int node_y);
void Is_assigner(double node_x, double node_y, double I, arma::mat &RHS);
void C_assigner_BE(int node_x,int node_y,double C, double h, arma::mat &LHS, arma::mat &RHS,const arma::mat &pre_solution, int mode);
int V_pulse_assigner(int node_x, int node_y, double V_value, arma::sp_mat &LHS, arma::mat &RHS);
double V_pulse_value(double V1, double V2, double t1, double td, double tr, double tf, double tpw, double tper);
void VCCS_assigner(int node_x, int node_y, int node_cx, int node_cy, double R, arma::sp_mat &LHS);
double Diode_assigner(int node_x, int node_y, double Is, double VT, arma::mat &LHS, arma::mat &RHS, arma::mat solution, int mode);
std::pair<arma::sp_mat,arma::mat> DynamicNonLinear(const CKTcircuit &ckt, double h, arma::mat pre_solution, int mode,const double time_trans, std::vector<Capacitor> &C_list, int NR_iteration_counter);
arma::mat NewtonRaphson_system(const CKTcircuit &ckt, const arma::mat &pre_solution, const double &h, const int &mode, const double time_trans, std::vector<Capacitor> &C_list);
void dummy_task();
/*  Global variables:

*/  
std::vector<Transient> vec_trans;
BS::thread_pool pool(3);
BS::timer tmr;
BS::synced_stream sync_out;
XB_Timer timer;
XB_Timer multi_timer;

/*--------------------------Timeing----------------------------------------------------------------*/
std::vector<std::chrono::time_point<std::chrono::high_resolution_clock>> begin_mid;
std::vector<std::chrono::time_point<std::chrono::high_resolution_clock>> begin_up;
std::vector<std::chrono::time_point<std::chrono::high_resolution_clock>> begin_down;
std::vector<std::chrono::time_point<std::chrono::high_resolution_clock>> end_mid;
std::vector<std::chrono::time_point<std::chrono::high_resolution_clock>> end_up;
std::vector<std::chrono::time_point<std::chrono::high_resolution_clock>> end_down;


/*-----------------------------------------------------------------------------------------------------------------------------------------------------*/

// Function to control debug mode
void setDebugMode(bool mode) {
    debugMode = mode;
}

bool isDebugMode() {
    return debugMode;
}

double cond(double R){
    return 1/R;
}

void history_trans_update(Transient &trans){
    // Update transient history
    vec_trans.push_back(trans);
}


class Truncation_error{
public:
    arma::mat LTE_bound_mid{};
    arma::mat LTE_bound_up{};
    arma::mat LTE_bound_down{};

    arma::mat LTE_BE_mid_h{};
    arma::mat LTE_BE_up_h{};
    arma::mat LTE_BE_down_h{};
};


class VoltageSource {
public:
    int id;
    int nodePos, nodeNeg;
    double value;
};

class Pulsevoltage{
public:
    int id{};
    int nodePos{}, nodeNeg{};
    double t1_pulse{};
    double V1{};
    double V2{};
    double td{};
    double tr{};
    double tf{};
    double pw{};
    double per{};

    int RHS_locate{};
};

class Diode {
public:
    int id{};
    int nodePos{}, nodeNeg{};
    double Is{};
    double VT{};

};

class VCCS{
public:
    int id{};
    int node_x{}, node_y{}, node_cx{}, node_cy{};
    double value{};
};


class NMOS{
public:
    int id{};
    int node_vd{}, node_vg{}, node_vs{}, node_vb{};
    double W{}, L{};
    
};

class PMOS{
public:
    int id{};
    int node_vd{}, node_vg{}, node_vs{}, node_vb{};
    double W{}, L{};
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
    int id{};
    std::string name{};      // It's used in MOSFETs Eg: M1.1, M1.2, M1.3, M1.4
    int nodePos{}, nodeNeg{};
    double value{};

    double current{};
    double voltage{};
    double charge{};
};

class multi_timestep{
public:
    // This class is used to store the data during multi-timestep algorithm

    arma::mat C_current_mid;
    arma::mat C_voltage_mid;
    arma::mat C_charge_mid;
    arma::mat Capacitance_mid;

    arma::mat C_current_up;
    arma::mat C_voltage_up;
    arma::mat C_charge_up;
    arma::mat Capacitance_up;

    arma::mat C_current_down;
    arma::mat C_voltage_down;
    arma::mat C_charge_down;
    arma::mat Capacitance_down;

    double h_mid{};
    double h_up{};
    double h_down{};   

    double t_mid{};
    double t_up{};
    double t_down{};

    arma::mat RHS_mid;
    arma::mat RHS_up;
    arma::mat RHS_down; 

    arma::mat LHS_mid;
    arma::mat LHS_up;
    arma::mat LHS_down;

    arma::mat solution_mid;
    arma::mat solution_up;
    arma::mat solution_down;

    bool TMAX_reach;
    bool TMIN_reach;

    bool ITE_reach_mid;
    bool ITE_reach_up;
    bool ITE_reach_down;

    std::vector<Capacitor> C_list_mid;
    std::vector<Capacitor> C_list_up;
    std::vector<Capacitor> C_list_down;

};

class CircuitElement {
public:
  std::variant<VoltageSource, CurrentSource, Resistor, Capacitor, Pulsevoltage, Diode, NMOS, PMOS, VCCS> element;
};

// Parser class

class CircuitParser {
public:
    std::string filename;
    std::vector<CircuitElement> elements;
    double double_t_end;  // This double_t_end can be passed to the CKTcircuit class
    double double_init_h;
    int num_mosfets{};

    CircuitParser(const std::string& filename) : filename(filename) {}
 
    void parser() {

        std::ifstream file(filename);

        if (!file.is_open()) {

            std::cerr << "Error opening file: " << filename << std::endl;

            exit(1);

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
                if constexpr(std::is_same_v<std::decay_t<decltype(arg)>, NMOS>){
                    maxNode = std::max({maxNode, arg.node_vd, arg.node_vg, arg.node_vs, arg.node_vb});
                
                }else if constexpr(std::is_same_v<std::decay_t<decltype(arg)>, PMOS>){
                    maxNode = std::max({maxNode, arg.node_vd, arg.node_vg, arg.node_vs, arg.node_vb});

                }else if constexpr(std::is_same_v<std::decay_t<decltype(arg)>, VCCS>){
                    maxNode = std::max({maxNode, arg.node_x, arg.node_y, arg.node_cx, arg.node_cy});
                }
                else{
                    maxNode = std::max(maxNode, arg.nodePos);
                    maxNode = std::max(maxNode, arg.nodeNeg);
                }

            }, element.element);
        }
        return maxNode;
    }
 
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
                std::string pulseParamsString;
                std::getline(iss, pulseParamsString);
                // Remove the parentheses
                pulseParamsString.erase(std::remove(pulseParamsString.begin(), pulseParamsString.end(), '('), pulseParamsString.end());
                pulseParamsString.erase(std::remove(pulseParamsString.begin(), pulseParamsString.end(), ')'), pulseParamsString.end());

                // Split the pulseParamsString into individual parameters
                std::istringstream pulseParamsStream(pulseParamsString);
                std::string v1, v2, td, tr, tf, pw, per;


                pv.id = v_id;
                pv.nodePos = v_nodePos;
                pv.nodeNeg = v_nodeNeg;

               pulseParamsStream >> v1 >> v2 >> td >> tr >> tf >> pw >> per;

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

        }else if(type[0] == 'D'){
            
            Diode d;

            d.id = std::stoi(type.substr(1));

            iss >> d.nodePos >> d.nodeNeg >> valueStr;

            d.Is = convertToValue(valueStr);

            iss >> valueStr;

            d.VT = convertToValue(valueStr);

            elements.push_back(CircuitElement{d});


        }else if(type[0] == 'G'){
            VCCS g;
            g.id = std::stoi(type.substr(1));

            iss >> g.node_x >> g.node_y >> g.node_cx >> g.node_cy >> valueStr;

            g.value = convertToValue(valueStr);

            elements.push_back(CircuitElement{g});
        }
        
        else if(type[0]=='M'){

            int M_id = std::stoi(type.substr(1));

            int M_node_vd, M_node_vg, M_node_vs, M_node_vb;

            std::string M_model, parameter;

            num_mosfets = num_mosfets + 1;

            iss >> M_node_vd >> M_node_vg >> M_node_vs >> M_node_vb >> M_model;

            if(M_model == "NMOS"){

                NMOS mn;

                mn.id = M_id;

                mn.node_vd = M_node_vd;

                mn.node_vg = M_node_vg;

                mn.node_vs = M_node_vs;

                mn.node_vb = M_node_vb;

                // Read and parse the W and L parameters with their prefixes
                while(iss >> parameter){
                    size_t pos = parameter.find('=');
                    if(pos != std::string::npos){
                        std::string key = parameter.substr(0, pos);
                        std::string value = parameter.substr(pos + 1);
                
                        if(key == "W"){
                            valueStr = value;
                            mn.W = convertToValue(valueStr);
                            
                        }else if(key == "L"){
                            valueStr = value;
                            mn.L = convertToValue(valueStr);

                            // std::cout<<"L is "<<mn.L<<std::endl;
                            // std::cout<<"type of L is "<<typeid(mn.L).name() << std::endl;
                        }
                    }
                }

                elements.push_back(CircuitElement{mn});
            }
            else if(M_model == "PMOS"){

                PMOS mp;

                mp.id = M_id;

                mp.node_vd = M_node_vd;

                mp.node_vg = M_node_vg;

                mp.node_vs = M_node_vs;

                mp.node_vb = M_node_vb;

                while(iss >> parameter){
                    size_t pos = parameter.find('=');
                    if(pos != std::string::npos){
                        std::string key = parameter.substr(0, pos);
                        std::string value = parameter.substr(pos + 1);
                
                        if(key == "W"){
                            valueStr = value;
                            mp.W = convertToValue(valueStr);
                            
                        }else if(key == "L"){
                            valueStr = value;
                            mp.L = convertToValue(valueStr);
                        }
                    }
                }

                elements.push_back(CircuitElement{mp});
            }
            else{
                std::cerr << "Error: Unknown MOSFET model: " << M_model << std::endl;
                exit(1);
            }
            

        }else if (type == ".tran") {
            std::string string_h, string_t_end;

            iss >> string_h >> string_t_end;

            double_init_h = convertToValue(string_h);
            double_t_end = convertToValue(string_t_end);
        
        }else {
                
            std::cerr << "Error: Unknown element type: " << type << std::endl;
            exit(1);
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

            case 'k': return value * 1000.0;      // Kilo 
            case 'm': return value * 1.0e-3;      // Milli
            case 'u': return value * 1.0e-6;      // Micro          
            case 'n': return value * 1.0e-9;      // Nano
            case 'p': return value * 1.0e-12;     // Pico
            case 'f': return value * 1.0e-15;     // Femto
            default: 
                std::cerr << "Error: Unknown unit: " << unit << std::endl;
                exit(1);
        }  
    }

    return value * 1.0; // No unit or unrecognized unit, assume the value is in base units

}

class DenseMatrix{
public:
    int Maxi{};                             // Size of matrix without branch current
    int Maxj{};

    arma::sp_mat LHS;                          // LHS matrix
    arma::mat RHS;                          // RHS matrix

    int n_rows{};                           // Number of rows and columns
    int n_cols{};
    
    void set_initmatrix(){
        init_LHS = LHS;
        init_RHS = RHS;
        n_rows = LHS.n_rows;
        n_cols = LHS.n_cols;
    }
    arma::sp_mat get_init_LHS() const{
        return init_LHS;
    }
    arma::mat get_init_RHS() const{
        return init_RHS;
    }

private:
    arma::sp_mat init_LHS;                     // The initial LHS and RHS values to be used in the NR-algorithm
    arma::mat init_RHS;

};

class Transient{
public:
    const int code = 0;                         // Choosing code for transient simulation
    const double t_start = 0;                   
    double t_end{};
    double h{};                                 // time step
    double init_h{};                            // initial time step
    double h_MAX{};
    double h_MIN{};
    double time_trans{};                        // transient time
    int mode{};                                   // 0 to do OP analysis, 1 to do transient simulation
    int trans_count{};                            // the count of the transient simulation loop
    arma::mat solution;                         // Solution matrix
    arma::sp_mat LHS;                              // LHS matrix
    arma::mat RHS;                              // RHS matrix
    arma::mat C_current;                        // Current matrix
    arma::mat C_voltage;                        // Voltage matrix
    arma::mat C_charge;                         // Charge matrix
    arma::mat Capacitance;                      // Capacitance matrix

    std::vector<Capacitor> C_list;              // Vector of capacitors in trans (including the parasitic capacitance of the MOSFETs)

    void C_list_update(){
        for(size_t i = 0; i < C_list.size(); i++){
            C_list.at(i).current = C_current(i,0);
            C_list.at(i).voltage = C_voltage(i,0);
            C_list.at(i).charge = C_list.at(i).value * C_list.at(i).voltage;
            C_list.at(i).value = Capacitance(i,0);
        }

    }

    void get_capacitance(){
        Capacitance = arma::zeros(C_list.size(),1);
        for(size_t i = 0; i < C_list.size(); i++){
            Capacitance(i,0) = C_list.at(i).value;
        }
    }
};

arma::mat get_capacitance(const std::vector<Capacitor> &C_list){
    arma::mat Capacitance = arma::zeros(C_list.size(),1);
    for(size_t i = 0; i < C_list.size(); i++){
        Capacitance(i,0) = C_list.at(i).value;
    }
    return Capacitance;
}


class CKTcircuit{
public:
    int const supply_voltage_node = 1;          // supply voltage node for the ring oscillator
    int pulse_num{};                            // Number of pulse voltages 
    // uword is a typedef for an unsigned integer type; it is used for matrix indices as well as all internal counters and loops                                           
    
    std::vector<CircuitElement> CKTelements;    // Vector of circuit elements
    int external_nodes{};                       // Number of external nodes (excluding ground and ring oscillator loop nodes)
    //int external_mosfets{};                   // Number of standalone mosfets (excluding mosfets from ring oscillator)
    int no_of_mosfets{};                        // Total number of MOSFETs
    int no_of_V_sources{};                      // Total number of voltage sources
    int T_nodes{};                              // Total number of nodes excluding ground

    DenseMatrix* cktdematrix;                   // Dense matrix class object
    void setcktmatrix(DenseMatrix& DenseMatrix){
        cktdematrix = &DenseMatrix;
    }

    std::vector<Capacitor> C_list;              // Vector of capacitors (including the parasitic capacitance of the MOSFETs)

    bool ckt_loaded{};                          // To check if the circuit is loaded or not
};


void CKTsetup(CKTcircuit &ckt, const CircuitParser &parser, DenseMatrix &DenseMatrix){  
    
    // Careful! getCircuitElements function is const, so it can't be used to modify the elements vector
    //ckt.elements = parser.getCircuitElements(); 
    ckt.CKTelements = parser.elements;  
    ckt.external_nodes = parser.getMaxNode();
    ckt.no_of_mosfets = parser.num_mosfets;
    ckt.T_nodes = ckt.external_nodes + 3 * ckt.no_of_mosfets; 
    // ckt.T_nodes = ckt.external_nodes;

    // Size of matrix
    DenseMatrix.Maxi = ckt.T_nodes;
    DenseMatrix.Maxj = DenseMatrix.Maxi;
    DenseMatrix.LHS = arma::sp_mat(DenseMatrix.Maxi,DenseMatrix.Maxj);   // LHS matrix
    DenseMatrix.RHS = arma::zeros(DenseMatrix.Maxi,1);    // RHS matrix
}

void CKTload(CKTcircuit &ckt){
    // ASSIGNING THE STAMPS TO THE LHS AND RHS MATRICES

    for (auto& element : ckt.CKTelements) {
        std::visit([&](auto&& arg) {
            if constexpr(std::is_same_v<std::decay_t<decltype(arg)>, Resistor>){
                std::cout << "R Element ID: " << arg.id << ", Node Pos: " << arg.nodePos << ", Node Neg: " << arg.nodeNeg << ", value: "<< arg.value << std::endl;

                if(arg.value == 0){

                    arg.value = 1e-3;
                }
                R_assigner(arg.nodePos, arg.nodeNeg, 1/arg.value, ckt.cktdematrix->LHS);

                

            }else if constexpr(std::is_same_v<std::decay_t<decltype(arg)>, VoltageSource>){
                std::cout << "VS Element ID: " << arg.id << ", Node Pos: " << arg.nodePos << ", Node Neg: " << arg.nodeNeg << ", value: "<< arg.value << std::endl;

                Vs_assigner(arg.nodePos, arg.nodeNeg, arg.value, ckt.cktdematrix->LHS, ckt.cktdematrix->RHS);


                ckt.no_of_V_sources++;

            }else if constexpr(std::is_same_v<std::decay_t<decltype(arg)>, Pulsevoltage>){
                std::cout << "VPulse Element ID: " << arg.id << ", Node Pos: " << arg.nodePos << ", Node Neg: " << arg.nodeNeg << std::endl;

                ckt.no_of_V_sources++;

                ckt.pulse_num++;

                arg.RHS_locate = V_pulse_assigner(arg.nodePos, arg.nodeNeg, arg.V1, ckt.cktdematrix->LHS, ckt.cktdematrix->RHS);

            }else if constexpr(std::is_same_v<std::decay_t<decltype(arg)>, CurrentSource>){
                std::cout << "I Element ID: " << arg.id << ", Node Pos: " << arg.nodePos << ", Node Neg: " << arg.nodeNeg << ", value: "<< arg.value << std::endl;

                Is_assigner(arg.nodePos, arg.nodeNeg, arg.value, ckt.cktdematrix->RHS);

            }else if constexpr(std::is_same_v<std::decay_t<decltype(arg)>, Capacitor>){
                std::cout << "C Element ID: " << arg.id << ", Node Pos: " << arg.nodePos << ", Node Neg: " << arg.nodeNeg << ", value: "<< arg.value << std::endl;

                ckt.C_list.push_back(arg);

            }else if constexpr(std::is_same_v<std::decay_t<decltype(arg)>, Diode>){
                std::cout << "D Element ID: " << arg.id << ", Node Pos: " << arg.nodePos << ", Node Neg: " << arg.nodeNeg << ", value: "<< arg.Is << std::endl;

            }else if constexpr(std::is_same_v<std::decay_t<decltype(arg)>, NMOS>){
                std::cout << "NMOS Element ID: " << arg.id << ", Node VD: " << arg.node_vd << ", Node VG: " << arg.node_vg << ", Node VS: " << arg.node_vs << ", Node VB: " << arg.node_vb << std::endl;

            }else if constexpr(std::is_same_v<std::decay_t<decltype(arg)>, PMOS>){
                std::cout << "PMOS Element ID: " << arg.id << ", Node VD: " << arg.node_vd << ", Node VG: " << arg.node_vg << ", Node VS: " << arg.node_vs << ", Node VB: " << arg.node_vb << std::endl;
            
            }else if constexpr(std::is_same_v<std::decay_t<decltype(arg)>, VCCS>){
                std::cout << "G Element ID: " << arg.id << ", Node X: " << arg.node_x << ", Node Y: " << arg.node_y << ", Node CX: " << arg.node_cx << ", Node CY: " << arg.node_cy << ", value: " << arg.value << std::endl;

                VCCS_assigner(arg.node_x, arg.node_y, arg.node_cx, arg.node_cy, arg.value, ckt.cktdematrix->LHS);
            }
               
            
        }, element.element);
    }
}

void Transsetup(Transient &trans, const CircuitParser &parser, CKTcircuit &ckt){
    trans.t_end = parser.double_t_end;
    

    if(ckt.pulse_num > 0 && (ckt.C_list.size() > 0)){
        for (const auto& element : ckt.CKTelements) {
            std::visit([&](auto&& arg) {

                if constexpr(std::is_same_v<std::decay_t<decltype(arg)>, Pulsevoltage>){
                
                    double step = std::min(arg.tr, arg.tf);
                    trans.h_MAX = parser.double_init_h;                             // 5 is rmax in spice opus
                    // trans.h_MIN = std::max(1e-9 * parser.double_init_h, 1e-14);    // 1e-9 is rmin in spice opus
                    trans.h_MIN =  1e-14; 

                    if(arg.td == 0){
                        trans.init_h = parser.double_init_h /100;   // initial time step, fs = 0.25 from spice opus
                           
                    }
                    else{
                        trans.init_h = std::min(arg.td * 0.25, parser.double_init_h * 0.25);
                    }
                    

                }
            
            }, element.element);
        }
    }
    else if(ckt.no_of_mosfets > 0){
        trans.init_h = parser.double_init_h / 100;     // initial time step
        trans.h_MAX = parser.double_init_h ;   // maximum time step
        // trans.h_MIN = std::max(1e-9 * parser.double_init_h, 1e-14);     // minimum time step
        trans.h_MIN =  1e-14; 
    }
    else{

        // trans.h_MAX = trans.t_end / 5000;
        trans.h_MAX = 1e-9;
        trans.init_h = trans.h_MAX;

    }

    
}

/*
    1. Start of rise time: t=td. This is when the voltage begins to transition from V1 to V2​.
    2. End of rise time: t=td+tr. This is when the voltage finishes transitioning to V2​.
    3. Start of fall time: t=td+tr+tpw. This is when the voltage begins to transition back to V1
    4​. End of fall time: t=td+tr+tpw+tf. This is when the voltage finishes transitioning back to V1​.
*/
std::deque<double> get_breakpoints(const CKTcircuit &ckt, const Transient &trans){

    std::deque<double> breakpoints;

    if(ckt.pulse_num == 0){
        return breakpoints;
    }

    for (const auto& element : ckt.CKTelements) {
        std::visit([&](auto&& arg) {
            if constexpr(std::is_same_v<std::decay_t<decltype(arg)>, Pulsevoltage>){
                
                for(double cycle_start = 0; cycle_start < trans.t_end; cycle_start += arg.per){

                    double cycle_times[] = {
                        arg.td+cycle_start, 
                        arg.td+arg.tr+cycle_start, 
                        arg.td+arg.tr+arg.pw+cycle_start, 
                        arg.td+arg.tr+arg.pw+arg.tf+cycle_start
                    };

                    for(double t : cycle_times){
                        if(t <= trans.t_end){
                            breakpoints.push_back(t);
                        }
                    }
                    
                }

            }
        }, element.element);
    }

    std::sort(breakpoints.begin(), breakpoints.end());

    if(breakpoints.at(0) == 0){
        breakpoints.pop_front();
    }

    if(breakpoints.back() != trans.t_end){
        breakpoints.push_back(trans.t_end);
    }

    return breakpoints;
}


// create resistor matrix stamp
void R_assigner(int node_x, int node_y, double G, arma::sp_mat &LHS){

    // if(G == 0){
    //     std::cerr << "Error: Resistor value cannot be zero" << std::endl;
    //     exit(1);
    // }

    if((node_x == 0) && (node_y == 0)){

        return;
    }
    else{
        if(node_x == 0){
            LHS(node_y-1, node_y-1) += G;
        }
        else if(node_y == 0){
            LHS(node_x-1, node_x-1) += G;
        }
        else{
            LHS(node_x-1, node_x-1) +=  G;
            LHS(node_x-1, node_y-1) += -G;
            LHS(node_y-1, node_x-1) += -G;
            LHS(node_y-1, node_y-1) +=  G;
        }
    }

}

// Voltage source stamp assigner
void Vs_assigner(int node_x, int node_y, double V_value, arma::sp_mat &LHS, arma::mat &RHS){

    arma::vec value(1);
    value = V_value;
    // Extending the branch at the LHS matrix
    LHS = branch_ext(LHS, node_x, node_y);

    // Assigning the value at RHS
    RHS = arma::join_cols(RHS, value);

}

// branch extender function
// arma::mat branch_ext(const arma::mat &M, int node_x, int node_y){
//     int M_cols = M.n_cols;
//     int M_rows = M.n_rows;
//     arma::mat va = arma::zeros(M_rows,1);
//     if(node_x == 0){
//         va.row(node_y-1).col(0) = 1;
//     }
//     else if(node_y == 0){
//         va.row(node_x-1).col(0) = 1;
//     }
//     else{
//         va.row(node_x-1).col(0) = -1;
//         va.row(node_y-1).col(0) = 1;
//     }

//     arma::mat zero_ext = arma::zeros(1,1);
//     arma::mat ha = va.as_row();
//     arma::mat haz = arma::join_rows(ha,zero_ext);

//     arma::mat M1 = arma::join_rows(M,va);
//     arma::mat M2 = arma::join_cols(M1,haz);

//     return M2;
// }
// Overload for spares matrixes (LHS)
arma::sp_mat branch_ext(const arma::sp_mat &M, int node_x, int node_y){
    // Create a sparse column vector va
    arma::sp_mat va(M.n_rows, 1);

    if (node_x == 0) {
        va(node_y - 1, 0) = 1;
    } else if (node_y == 0) {
        va(node_x - 1, 0) = 1;
    } else {
        va(node_x - 1, 0) = -1;
        va(node_y - 1, 0) = 1;
    }

    // Create a sparse row vector ha (transpose of va)
    arma::sp_mat ha = va.t();

    // Add the zero element to ha to form haz
    ha.resize(1, M.n_cols + 1);

    // Create M1 by horizontally joining M and va
    arma::sp_mat M1 = arma::join_horiz(M, va);

    // Create the final matrix M2 by vertically joining M1 and ha
    arma::sp_mat M2 = arma::join_vert(M1, ha);

    return M2;
}

// create a current matrix stamp
void Is_assigner(double node_x, double node_y, double I, arma::vec &RHS){

    if((node_x == 0) && (node_y == 0)){
        I = 0;
    }
    else{
        if(node_x == 0){
            RHS.row(node_y-1).col(0) += I;
        }
        else if(node_y == 0){
            RHS.row(node_x-1).col(0) += -I;
        }
        else{
            RHS.row(node_x-1).col(0) += -I;
            RHS.row(node_y-1).col(0) += I;
        }
    }
    
}

// Capacitor stamp assigner with backward euler method
void C_assigner_BE(int node_x,int node_y,double C, double h, arma::sp_mat &LHS, arma::mat &RHS, const arma::mat &pre_solution, int mode){
    
    if(node_x == 0 && node_y == 0){
        return;
    }
    double x{};  // x = C/h
    double vol{};
    
        if(mode > 0){
            x = C / h;
            
            if(node_x == 0){
                vol = pre_solution(node_y - 1, 0);
            }else if(node_y == 0){
                vol = pre_solution(node_x - 1, 0);
            }else{
                vol = pre_solution(node_x-1,0) - pre_solution(node_y-1,0);
            }
        }
        else{
            x = 0;
            vol = 0;
        }
    // Matrix stamp for a capacitor on LHS
    R_assigner(node_x, node_y, x, LHS);
    // Matrix stamp for a capacitor on RHS
    Is_assigner(node_x, node_y, -(x * vol), RHS);
}

// Correct Diode stamp assigner (Original Diode_assigner is wrong with the node voltage assignment)
double Diode_assigner(int node_x, int node_y, double Is, double VT, arma::sp_mat &LHS, arma::mat &RHS, 
                      arma::mat solution, int mode)
{
    int maxi = LHS.n_cols;
    int maxj = 1;

    double G_eq = 0;
    double I_eq = 0;
    double val_nodex = 0;
    double val_nodey = 0;

    double Id = 0;
    double vol = 0; // Vd(k)

    // The Companion Model of a Diode is Resistor and Current Source in parallel.
    // Id(k+1) = G_ed * Vd(k+1) + I_eq(k)
    // G_eq = Id' or g'(vd) = Isat/Vt * e^(vd(k)/Vt)
    // I_eq = Id - G_eq * Vd(k)


    if ((node_x == 0 && node_y == 0))
    {
        G_eq = 0;
        I_eq = 0;
        vol = 0;
        Id = Is * (exp(vol / VT) - 1);
        return Id;
    }
    else if (node_x == 0)
    {
        val_nodey = solution(node_y - 1, 0);
        vol = -val_nodey;
        Id = Is * (exp(vol / VT) - 1);
        G_eq = (Is / VT) * (exp(vol / VT));
        I_eq = Id - (G_eq * vol);
    }
    else if (node_y == 0)
    {
        val_nodex = solution(node_x - 1, 0);
        vol = val_nodex;
        Id = Is * (exp(vol / VT) - 1);
        G_eq = (Is / VT) * (exp(vol / VT));
        I_eq = Id - (G_eq * vol);
    }
    else
    {
        val_nodex = solution(node_x - 1, 0);
        val_nodey = solution(node_y - 1, 0);
        vol = val_nodex - val_nodey;
        Id = Is * (exp(vol / VT) - 1);
        G_eq = (Is / VT) * (exp(vol / VT));
        I_eq = Id - (G_eq * vol);
    }

    // Matrix stamp for a diode on RHS
    Is_assigner(node_x, node_y, I_eq, RHS);
    // Matrix stamp for a diode on LHS
    R_assigner(node_x, node_y, G_eq, LHS);

    
    return Id;
}



int V_pulse_assigner(int node_x, int node_y, double V_value, arma::sp_mat &LHS, arma::mat &RHS){

    arma::vec value(1);
    value = V_value;
    // Extending the branch at the LHS matrix
    LHS = branch_ext(LHS, node_x, node_y);

    // Assigning the value at RHS
    RHS = arma::join_cols(RHS, value);

    // location of voltage value for transient simulation
    int Pulse_RHS_locate1 = RHS.n_rows;

    return Pulse_RHS_locate1;
}

double V_pulse_value(double V1, double V2, double t1, double td, double tr, double tf, double tpw, double tper){
    double v{};
    t1 = fmod(t1, tper);
    if(tper < tr + tpw + tf){
        std::cerr << "Period is incorrect" << std::endl;
        exit(1);
    }

    if(t1 < 0){
        std::cerr << "Simulation time in pulse voltage it wrong" << std::endl;
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
        std::cerr << "Pulse voltage Error" << std::endl;
        exit(1);
    }
    return v;
}

// assigning the matrix stamps for the VCCS
void VCCS_assigner(int node_x, int node_y, int node_cx, int node_cy, double R, arma::sp_mat &LHS)
{

    if (node_x == 0)
    {
        if (node_cx == 0)
        {
            if (node_cy > 0)
            {
                LHS(node_y - 1, node_cy - 1) += R;
            }
        }
        else if (node_cy == 0)
        {
            if (node_cx > 0)
            {
                LHS(node_y - 1, node_cx - 1) += -R;
            }
        }
        else
        {
            LHS(node_y - 1, node_cx - 1) += -R;
            LHS(node_y - 1, node_cy - 1) += R;
        }
    }
    else if (node_y == 0)
    {
        if (node_cx == 0)
        {
            if (node_cy > 0)
            {
                LHS(node_x - 1, node_cy - 1) += -R;
            }
        }
        else if (node_cy == 0)
        {
            if (node_cx > 0)
            {
                LHS(node_x - 1, node_cx - 1) += R;
            }
        }
        else
        {
            LHS(node_x - 1, node_cx - 1) += R;
            LHS(node_x - 1, node_cy - 1) += -R;
        }
    }
    else
    {
        if (node_cx == 0)
        {
            if (node_cy > 0)
            {
                LHS(node_x - 1, node_cy - 1) += -R;
                LHS(node_y - 1, node_cy - 1) += R;
            }
        }
        else if (node_cy == 0)
        {
            if (node_cx > 0)
            {
                LHS(node_x - 1, node_cx - 1) += R;
                LHS(node_y - 1, node_cx - 1) += -R;
            }
        }
        else
        {
            LHS(node_x - 1, node_cx - 1) += R;
            LHS(node_x - 1, node_cy - 1) += -R;
            LHS(node_y - 1, node_cx - 1) += -R;
            LHS(node_y - 1, node_cy - 1) += R;
        }
    }


}

double NMOS_assigner(int number, int node_vd, int node_vg, int node_vs, int node_vb, double W, double L, double h,
                    const arma::mat &solution , int T_nodes, arma::sp_mat &LHS, arma::mat &RHS, int mode, std::vector<Capacitor> &C_list, int NR_iteration_counter)
{
    

    double vd = 0.0;
    double vg = 0.0;
    double vs = 0.0;
    double vb = 0.0;

    // The settings for the large signal analysis model
    if (node_vd == 0)
    {
        vd = 0.0;
    }
    else
        vd = solution(node_vd - 1, 0);

    if (node_vg == 0)
    {
        vg = 0.0;
    }
    else
        vg = solution(node_vg - 1, 0);

    if (node_vs == 0)
    {
        vs = 0.0;
    }
    else
        vs = solution(node_vs - 1, 0);

    if (node_vb == 0)
    {
        vb = 0.0;
    }
    else
        vb = solution(node_vb - 1, 0);
    // enhancement mode voltages
    double vgs = vg - vs;
    double vds = vd - vs;
    double vbs = vb - vs;
    // depletion mode voltages
    double vbd = vb - vd;
    double vgd = vg - vd;

    double id = 0.0;
    double gds = 0.0;
    double gm = 0.0;
    double gmb = 0.0;
    double Ld = 0.08e-6;
    double Leff = L;
    double kp = 2.0e-5; // default is 2e-5
    double mCox = kp;
    double LAMBDA = 0.1;

    double Beta = (mCox) * (W / Leff);
    double gamma = 0.37;
    double phi = 0.65;
    double vt0 = 0.7; // default is 0
    double vt = 0.0;
    double I_DSeq = 0.0;

    // Capacitances settings
    double CGSO = 4.0e-11;
    double CGDO = 4.0e-11;
    double CGBO = 2.0e-11;
    double CGS = CGSO * W;
    double CGD = CGDO * W;
    double CGB = CGBO * L;
    double CBD = 20.0e-15; // typical value for CBD
    double CBS = 20.0e-15; // typical value for CBS

    double FC = 0.5;     // Coefficient for forward-bias depletion capacitance formula
    double PB = 0.9;     // Bulk junction potential
    double CJ = 0.56e-3;    // Zero bias bulk junction capacitance per unit area
    double MJ = 0.45;     // Bulk junction bottom grading coefficient
    double AD = 200.0e-12; // Drain area
    double AS = 200.0e-12; // Source area
    double CJSW = 0.35e-11; // Zero bias bulk junction sidewall capacitance per unit periphery
    double PD = 20.0e-6;   // Drain junction potential
    double PS = 20.0e-6;   // Source junction potential
    double TT = 1.0e-9;    // Transit time
    double MJSW = 0.2;   // Bulk junction sidewall grading coefficient level 1
    double JSSW = 1.0e-9;  // Bulk junction saturation current per meter of sidewall
    double JS = 1.0e-8;    // Bulk junction saturation current per meter of junction perimeter

    if (vbd <= FC * PB)
    {
        CBD = (CJ * AD) / (pow((1.0 - vbd / PB), MJ)) + (CJSW * PD) / (pow((1.0 - vbd / PB), MJSW));
    }
    else
    {
        CBD = ((CJ * AD) * (1.0 - (1.0 + MJ) * FC + MJ * vbd / PB)) / (pow((1.0 - FC), (1.0 + MJ))) + (CJSW * PD) * (1.0 - (1.0 + MJSW) * FC + MJSW * vbd / PB) / (pow((1.0 - FC), (1.0 + MJSW)));
    }
    if (vbs <= FC * PB)
    {
        CBS = (CJ * AS) / (pow((1.0 - vbs / PB), MJ)) + (CJSW * PS) / (pow((1.0 - vbs / PB), MJSW));
    }
    else
    {
        CBS = ((CJ * AS) * (1.0 - (1.0 + MJ) * FC + MJ * vbs / PB)) / (pow((1.0 - FC), (1.0 + MJ))) + (CJSW * PS) * (1.0 - (1.0 + MJSW) * FC + MJSW * vbs / PB) / (pow((1.0 - FC), (1.0 + MJSW)));
    }

    // # the settings for fet model based on the large signal analysis
    

        R_assigner(node_vd, T_nodes - (3 * number) + 1, 1.0/RD, LHS); // # RD
        R_assigner(node_vg, T_nodes - (3 * number) + 2, 1.0/RG, LHS);   // # RG
        R_assigner(T_nodes - (3 * number) + 3, node_vs, 1.0/RS, LHS); // # RS
        // R_assigner(node_vb ,T_nodes - (4 * number) + 4, 1, LHS, RHS);   // # RB

        // // New model for Diode
        // Diode_assigner(T_nodes - (4 * number) + 3, node_vd, 1e-14, 0.05, LHS, RHS, solution, mode); // # Diode BD
        // Diode_assigner(T_nodes - (4 * number) + 4, node_vs, 1e-14, 0.05, LHS, RHS, solution, mode); // # Diode BS
        // Diode_assigner(node_vb, node_vd, 1e-14, 0.05, LHS, RHS, solution, mode); // # Diode BD
        // Diode_assigner(node_vb, node_vs, 1e-14, 0.05, LHS, RHS, solution, mode); // # Diode BS

        
    

    vt = vt0 + gamma * ((sqrt(phi - vbs) - sqrt(phi))); // # already taking into account the body effect of MOSFETs
    double n_vt = std::abs(vt); 

    if ((vds <= (vgs - vt)) && (vgs >= vt))
    { // # the transistor is in linear
        id = Beta * (vgs - vt - vds / 2.0) * vds * (1.0 + LAMBDA * vds);
        // gds = Beta/2.0 * LAMBDA * pow(vds,2) + Beta * (vgs - n_vt - vds) * (1 + 2 * LAMBDA * vds);  // From the book: Circuit simulation
        gds = Beta * (1.0 + LAMBDA * vds) * (vgs - vt - vds) + Beta * LAMBDA * vds * (vgs - vt - vds / 2.0);
        gm = Beta * vds * (1.0 + LAMBDA * vds);
        gmb = gm * gamma / (2.0 * sqrt(phi - vbs));

        CGS = 2.0/3.0 * mCox * (1.0 - pow(vgs - vds -vt,2)/pow(2.0*(vgs -vt) - vds,2)) + CGSO * W;
        CGD = 2.0/3.0 * mCox * (1.0 - pow(vgs -vt,2)/pow(2.0*(vgs -vt) - vds,2)) + CGDO * W;
        CGB = 0.0 + CGBO * L;    
    }
    else if ((vds >= (vgs - vt)) && (vgs >= vt))
    { // # the transistor is in saturation
        id = (Beta / 2.0) * pow((vgs - vt), 2) * (1.0 + LAMBDA * vds);
        gds = (Beta / 2.0) * LAMBDA * pow((vgs - vt), 2);
        gm = Beta * (1 + LAMBDA * vds) * (vgs - vt);
        gmb = gm * gamma / (2.0 * sqrt(phi - vbs));

        CGS = 2.0/3.0 * mCox + CGSO * W;
        CGD = CGDO * W;
        CGB = 0.0 + CGBO * L;
    }
    else{
        id = 0.0;
        gds = 0.0;
        gm = 0.0;
        gmb = 0.0;
        // CGS = 0 + CGSO * W;
        // CGD = 0 + CGDO * W;
        // CGB = mCox + CGBO * L;

            if(vgs - vt <= -phi){
                CGS = 0.0 + CGSO * W;
                CGD = 0.0 + CGDO * W;
                CGB = mCox + CGBO * L;
            }else if(vgs - vt > -phi && vgs - vt <= -phi/2){
                CGS = 0.0 + CGSO * W;
                CGD = 0.0 + CGDO * W;
                CGB = -mCox * ((vgs - vt)/phi) + CGBO * L;
            }else{
                CGS = 2.0/3.0 * mCox + 4.0/3.0 * mCox * (vgs - vt)/phi + CGSO * W;
                CGD = 0.0 + CGDO * W;
                CGB = -mCox * ((vgs - vt)/phi) + CGBO * L;
            }

    }

    // gds cannot be zero.
    if(gds == 0){
        gds = 1.0e-12;
    }
    // if(id == 0){
    //     id = 1.0e-12;
    // }
    if(gm == 0){
        gm = 1.0e-12;
    }
    
    I_DSeq = id - gds * vds - gm * vgs; // # 10.190 equation


    
        // C_assigner_BE(node_vg, T_nodes - (4 * number) + 2, CGD, h, LHS, RHS, solution, mode); // # Capacitor CGD
        // C_assigner_BE(node_vg, T_nodes - (4 * number) + 1, CGS, h, LHS, RHS, solution, mode); // # Capacitor CGS
        // C_assigner_BE(T_nodes - (4 * number) + 3, node_vd, CBD, h, LHS, RHS, solution, mode); // # Capacitor CBD
        // C_assigner_BE(T_nodes - (4 * number) + 4, node_vs, CBS, h, LHS, RHS, solution, mode); // # Capacitor CBS



        // C_assigner_BE(node_vg, node_vd, CGD, h, LHS, RHS, solution, mode); // # Capacitor CGD
        // C_assigner_BE(node_vg, node_vs, CGS, h, LHS, RHS, solution, mode); // # Capacitor CGS
        // if(node_vb != 0 || node_vs != 0){
        //     C_assigner_BE(node_vb, node_vs, CBS, h, LHS, RHS, solution, mode); // # Capacitor CBS
        // }
        // C_assigner_BE(node_vg, node_vb, CGB, h, LHS, RHS, solution, mode); // # Capacitor CGB
        // C_assigner_BE(node_vb, node_vd, CBD, h, LHS, RHS, solution, mode); // # Capacitor CBD
        


        Is_assigner(T_nodes - (3 * number) + 1, T_nodes - (3 * number) + 3, I_DSeq, RHS);
        R_assigner(T_nodes - (3 * number) + 1, T_nodes - (3 * number) + 3, gds, LHS);           // # assigning gds
        VCCS_assigner(T_nodes - (3 * number) + 1, T_nodes - (3 * number) + 3, T_nodes - (3 * number) + 2, T_nodes - (3 * number) + 3, gm, LHS);  // # assigning gm


        // VCCS_assigner(node_vd, node_vs, node_vb, node_vs, gmb, LHS); // assigning gmb
        // std::cout << "gds_NMOS: " << gds << std::endl;  
        // std::cout << "gm_NMOS: " << gm << std::endl;
        // std::cout << "id_eq_NMOS: " << I_DSeq << std::endl;
        // std::cout << "id_NMOS: " << id << std::endl;

        DEBUG_PRINT("NMOS gds: " << gds);
        DEBUG_PRINT("NMOS gm: " << gm);
        DEBUG_PRINT("NMOS id_eq: " << I_DSeq);
        DEBUG_PRINT("NMOS id: " << id);


    // if(mode == 0 && NR_iteration_counter == 0){

        // Capacitor c1{0, "M"+std::to_string(number)+".1", node_vg, node_vd, CGD}; // M1.1
        // Capacitor c2{0, "M"+std::to_string(number)+".2", node_vg, node_vs, CGS}; // M1.2
        // if(node_vb != 0 || node_vs != 0){
        //     Capacitor c3{0, "M"+std::to_string(number)+".3", node_vb, node_vs, CBS}; // M1.3
        //     C_list.push_back(c3);
        // }
        // Capacitor c4{0, "M"+std::to_string(number)+".4", node_vg, node_vb, CGB}; // M1.4
        // Capacitor c5{0, "M"+std::to_string(number)+".5", node_vb, node_vd, CBD}; // M1.5


        // C_list.push_back(c1);
        // C_list.push_back(c2);

        // C_list.push_back(c4);
        // C_list.push_back(c5);
    // }
    // else{

    //     for(auto& Cap: C_list){                                      // Update the Capacitance values
    //         if(Cap.name == "M"+std::to_string(number)+".1"){
    //             Cap.value = CGD;
    //         }
    //         else if(Cap.name == "M"+std::to_string(number)+".2"){
    //             Cap.value = CGS;
    //         }
    //         // else if(Cap.name == "M"+std::to_string(number)+".3"){
    //         //     Cap.value = CBS;
    //         // }
    //         else if(Cap.name == "M"+std::to_string(number)+".4"){
    //             Cap.value = CGB;
    //         }
    //         // else if(Cap.name == "M"+std::to_string(number)+".5"){
    //         //     Cap.value = CBD;
    //         // }

    //     }
    // }
    
    
    return id;
}



double PMOS_assigner(int number, int node_vd, int node_vg, int node_vs, int node_vb, double W, double L, double h,
                    const arma::mat &solution , int T_nodes, arma::sp_mat &LHS, arma::mat &RHS, int mode, std::vector<Capacitor> &C_list, int NR_iteration_counter)
{
    double vd = 0.0;
    double vg = 0.0;
    double vs = 0.0;
    double vb = 0.0;

    // The settings for the large signal analysis model
    if (node_vd == 0)
    {
        vd = 0.0;
    }
    else
        vd = solution(node_vd - 1, 0);

    if (node_vg == 0)
    {
        vg = 0.0;
    }
    else
        vg = solution(node_vg - 1, 0);

    if (node_vs == 0)
    {
        vs = 0.0;
    }
    else
        vs = solution(node_vs - 1, 0);

    if (node_vb == 0)
    {
        vb = 0.0;
    }
    else
        vb = solution(node_vb - 1, 0);

    // enhancement mode voltages
    double vgs = vg - vs;
    double vds = vd - vs;
    double vbs = vb - vs;
    // enhancement mode voltages
    double vsg = vs - vg;
    double vsd = vs - vd;
    double vsb = vs - vb;
    // depletion mode voltages
    double vdb = vd - vb;
    double vdg = vd - vg;
    double vbd = vb - vd;
    double vgd = vg - vd;

    double id = 0.0;
    double gds = 0.0;
    double gm = 0.0;
    double gmb = 0.0;
    double Ld = 0.0;
    double Leff = L;
    double kp = 2.0e-5; // default is 2e-5
    double mCox = kp;
    double LAMBDA = 0.1;

    double Beta = (mCox) * (W / Leff);
    double gamma = 0.4;
    double phi = 0.8;
    double vt0 = -0.8; // default is 0
    double vt = 0.0;
    double I_DSeq = 0.0;

    // Capacitances settings
    double CGSO = 4.0e-10;
    double CGDO = 4.0e-10;
    double CGBO = 2.0e-10;
    double CGS = CGSO * W;
    double CGD = CGDO * W;
    double CGB = CGBO * L;
    double CBD = 20.0e-15; // typical value for CBD
    double CBS = 20.0e-15; // typical value for CBS

    double FC = 0.5;     // Coefficient for forward-bias depletion capacitance formula
    double PB = 0.8;     // Bulk junction potential
    double CJ = 2.0e-4;    // Zero bias bulk junction capacitance per unit area
    double MJ = 0.5;     // Bulk junction bottom grading coefficient
    double AD = 200.0e-12; // Drain area
    double AS = 200.0e-12; // Source area
    double CJSW = 1.0e-12; // Zero bias bulk junction sidewall capacitance per unit periphery
    double PD = 20.0e-6;   // Drain junction potential
    double PS = 20.0e-6;   // Source junction potential
    double TT = 1.0e-9;    // Transit time
    double MJSW = 0.5;   // Bulk junction sidewall grading coefficient level 1
    double JSSW = 1.0e-9;  // Bulk junction saturation current per meter of sidewall
    double JS = 1.0e-6;    // Bulk junction saturation current per meter of junction perimeter

    if (vdb <= FC * PB)
    {
        CBD = (CJ * AD) / (pow((1 - vdb / PB), MJ)) + (CJSW * PD) / (pow((1 - vdb / PB), MJSW));
    }
    else
    {
        CBD = ((CJ * AD) * (1 - (1 + MJ) * FC + MJ * vdb / PB)) / (pow((1 - FC), (1 + MJ))) + (CJSW * PD) * (1 - (1 + MJSW) * FC + MJSW * vdb / PB) / (pow((1 - FC), (1 + MJSW)));
    }
    if (vsb <= FC * PB)
    {
        CBS = (CJ * AS) / (pow((1 - vsb / PB), MJ)) + (CJSW * PS) / (pow((1 - vsb / PB), MJSW));
    }
    else
    {
        CBS = ((CJ * AS) * (1 - (1 + MJ) * FC + MJ * vsb / PB)) / (pow((1 - FC), (1 + MJ))) + (CJSW * PS) * (1 - (1 + MJSW) * FC + MJSW * vsb / PB) / (pow((1 - FC), (1 + MJSW)));
    }

    // # the settings for fet model based on the large signal analysis
    
        R_assigner(T_nodes - (3 * number) + 1, node_vd, 1.0/RD, LHS); // # RD
        R_assigner(node_vg, T_nodes - (3 * number) + 2, 1.0/RG, LHS);   // # RG
        R_assigner(node_vs, T_nodes - (3 * number) + 3, 1.0/RS, LHS); // # RS
        // R_assigner(T_nodes - (4 * number) + 4, node_vb, 1, LHS, RHS);   // # RB

        
        // Diode_assigner(T_nodes - (4 * number) + 2, T_nodes - (4 * number) + 3, 1e-14, 0.05, LHS, RHS,  solution, mode); // # Diode DB
        // Diode_assigner(T_nodes - (4 * number) + 4, T_nodes - (4 * number) + 3, 1e-14, 0.05, LHS, RHS,  solution, mode); // # Diode SB
        // Diode_assigner(node_vd, node_vb, 1e-14, 0.05, LHS, RHS, solution, mode); // # Diode BD
        // Diode_assigner(node_vs, node_vb, 1e-14, 0.05, LHS, RHS, solution, mode); // # Diode BS


        
    

    vt = vt0 - gamma * ((sqrt(phi - vsb) - sqrt(phi))); // already taking into account the body effect of MOSFETs

    double n_vt = std::abs(vt);

    if ((vds >= (vgs - vt)) && (vgs <= vt))
    { // # the transistor is in linear
        id = Beta * (vsg - n_vt - vsd / 2.0) * vsd * (1.0 + LAMBDA * vsd);
        gds = Beta * (1.0 + LAMBDA * vsd) * (vsg - n_vt - vsd) + Beta * LAMBDA * vsd * (vsg - n_vt - vsd / 2.0);
        gm = Beta * vsd * (1.0 + LAMBDA * vsd);
        gmb = gm * gamma / (2.0 * sqrt(phi - vsb));

       CGS = 2.0/3.0 * mCox * (1.0 - pow(vsg - vsd - n_vt,2)/pow(2*(vsg - n_vt) - vsd,2)) + CGSO * W;
       CGD = 2.0/3.0 * mCox * (1.0 - pow(vsg - n_vt,2)/pow(2*(vsg - n_vt) - vsd,2)) + CGDO * W;
       CGB = 0 + CGBO * L;
    }
    else if ((vds <= (vgs - vt)) && (vgs <= vt))
    { // # the transistor is in saturation
        id = (Beta / 2.0) * pow((vsg - n_vt), 2) * (1 + LAMBDA * vsd);
        gds = (Beta / 2.0) * LAMBDA * pow((vsg - n_vt), 2);
        gm = Beta * (1.0 + LAMBDA * vsd) * (vsg - n_vt);
        gmb = gm * gamma / (2.0 * sqrt(phi - vsb));

        CGS = 2.0/3.0 * mCox + CGSO * W;
        CGD = 0.0 + CGDO * W;
        CGB = 0.0 + CGBO * L;
    }
    else
    { // # the transistor is in cutoff
        id = 0.0; 
        gds = 0.0;
        gm = 0.0;
        gmb = 0.0;

        // CGS = 0 + CGSO * W;
        // CGD = 0 + CGDO * W;
        // CGB = mCox + CGBO * L;

        if(vsg - n_vt <= -phi){
            CGS = 0 + CGSO * W;
            CGD = 0 + CGDO * W;
            CGB = mCox + CGBO * L;
        }
        else if(vsg - n_vt > -phi && vsg - n_vt <= -phi/2){
            CGS = 0 + CGSO * W;
            CGD = 0 + CGDO * W;
            CGB = -mCox * ((vsg - n_vt)/phi) + CGBO * L;
        }
        else{
            CGS = 2.0/3.0 * mCox + 4.0/3.0 * mCox * (vsg - n_vt)/phi + CGSO * W;
            CGD = 0.0 + CGDO * W;
            CGB = -mCox * ((vsg - n_vt)/phi) + CGBO * L;
        }
    }

    
        // C_assigner_BE(T_nodes - (4 * number) + 2, T_nodes - (4 * number) + 3, CBD, h, LHS, RHS, solution, mode); // # Capacitor CBD
        // C_assigner_BE(T_nodes - (4 * number) + 4, T_nodes - (4 * number) + 3, CBS, h, LHS, RHS, solution, mode); // # Capacitor CBS
        // C_assigner_BE(T_nodes - (4 * number) + 1, T_nodes - (4 * number) + 4, CGS, h, LHS, RHS, solution, mode); // # Capacitor GSO
        // C_assigner_BE(T_nodes - (4 * number) + 1, T_nodes - (4 * number) + 3, CGB, h, LHS, RHS, solution, mode); // # Capacitor GBO
        // C_assigner_BE(T_nodes - (4 * number) + 4, T_nodes - (4 * number) + 3, CGD, h, LHS, RHS, solution, mode); // # Capacitor GDO

        // C_assigner_BE(node_vd, node_vg, CGD, h, LHS, RHS, solution, mode); // # Capacitor CGD
        // C_assigner_BE(node_vs, node_vg, CGS, h, LHS, RHS, solution, mode); // # Capacitor CGS
        // C_assigner_BE(node_vb, node_vg, CGB, h, LHS, RHS, solution, mode); // # Capacitor CGB
        // C_assigner_BE(node_vd, node_vb, CBD, h, LHS, RHS, solution, mode); // # Capacitor CBD
        // C_assigner_BE(node_vs, node_vb, CBS, h, LHS, RHS, solution, mode); // # Capacitor CBS
    

    // gds cannot be zero.
    if(gds == 0){
        gds = 1.0e-12;
    }
    // if(id == 0){
    //     id = 1.0e-12;
    // }
    if(gm == 0){
        gm = 1.0e-12;
    }

    I_DSeq = id - gds * vsd - gm * vsg;



        Is_assigner(T_nodes - (3 * number) + 3, T_nodes - (3 * number) + 1, I_DSeq, RHS);
        // VCCS_assigner(node_vs, node_vd, node_vs, node_vb, gmb, LHS); // assigning gmb
        R_assigner(T_nodes - (3 * number) + 3, T_nodes - (3 * number) + 1, gds, LHS);           // # assigning gds
        VCCS_assigner(T_nodes - (3 * number) + 3, T_nodes - (3 * number) + 1, T_nodes - (3 * number) + 3, T_nodes - (3 * number) + 2, gm, LHS);  // # assigning gm

        // std::cout << "gds_PMOS: " << gds << std::endl;  
        // std::cout << "gm_PMOS: " << gm << std::endl;
        // std::cout << "id_eq_PMOS: " << I_DSeq << std::endl;
        // std::cout << "id_PMOS: " << id << std::endl;

        DEBUG_PRINT("PMOS gds: " << gds);
        DEBUG_PRINT("PMOS gm: " << gm);
        DEBUG_PRINT("PMOS id_eq: " << I_DSeq);
        DEBUG_PRINT("PMOS id: " << id);



    // if(mode == 0 && NR_iteration_counter == 0){

    //     Capacitor c1{0, "M"+std::to_string(number)+".1", node_vd, node_vg, CGD}; // M1.1
    //     Capacitor c2{0, "M"+std::to_string(number)+".2", node_vs, node_vg, CGS}; // M1.2
    //     Capacitor c3{0, "M"+std::to_string(number)+".3", node_vb, node_vg, CGB}; // M1.3
    //     Capacitor c4{0, "M"+std::to_string(number)+".4", node_vd, node_vb, CBD}; // M1.4
    //     Capacitor c5{0, "M"+std::to_string(number)+".5", node_vs, node_vb, CBS}; // M1.5

    //     C_list.push_back(c1);
    //     C_list.push_back(c2);
    //     C_list.push_back(c3);
    //     C_list.push_back(c4);
    //     C_list.push_back(c5);
    // }
    // else{

    //     for(auto& Cap: C_list){                                      // Update the Capacitance values
    //         if(Cap.name == "M"+std::to_string(number)+".1"){
    //             Cap.value = CGD;
    //         }
    //         else if(Cap.name == "M"+std::to_string(number)+".2"){
    //             Cap.value = CGS;
    //         }
    //         else if(Cap.name == "M"+std::to_string(number)+".3"){
    //             Cap.value = CGB;
    //         }
    //         else if(Cap.name == "M"+std::to_string(number)+".4"){
    //             Cap.value = CBD;
    //         }
    //         else if(Cap.name == "M"+std::to_string(number)+".5"){
    //             Cap.value = CBS;
    //         }

    //     }
    // }
    

    return id;
}



// Assigning the stamp matrices for dynamic and non-linear components and update the LHS and RHS matrices
std::pair<arma::sp_mat,arma::mat> DynamicNonLinear(const CKTcircuit &ckt, double h, arma::mat pre_solution, int mode,const double time_trans, std::vector<Capacitor> &C_list, int NR_iteration_counter){
   
    arma::sp_mat LHS = ckt.cktdematrix->get_init_LHS();
    arma::mat RHS = ckt.cktdematrix->get_init_RHS();

    for (const auto& element : ckt.CKTelements) {
        std::visit([&](auto&& arg) {
            if constexpr(std::is_same_v<std::decay_t<decltype(arg)>, Capacitor>){
                // Linear Capacitor
                C_assigner_BE(arg.nodePos, arg.nodeNeg, arg.value, h, LHS, RHS, vec_trans.back().solution, mode);

            }
            else if constexpr(std::is_same_v<std::decay_t<decltype(arg)>, Pulsevoltage>){
                
                double val_pulse = V_pulse_value(arg.V1, arg.V2, time_trans, arg.td, arg.tr, arg.tf, arg.pw, arg.per);
                
                RHS.row(arg.RHS_locate-1).col(0) += val_pulse;

            }
            else if constexpr(std::is_same_v<std::decay_t<decltype(arg)>, NMOS>){
                NMOS_assigner(arg.id, arg.node_vd, arg.node_vg, arg.node_vs, arg.node_vb, arg.W, arg.L, h, pre_solution, ckt.T_nodes, LHS, RHS, mode, C_list, NR_iteration_counter);
            }
            else if constexpr(std::is_same_v<std::decay_t<decltype(arg)>, Diode>){
                Diode_assigner(arg.nodePos, arg.nodeNeg, arg.Is, arg.VT, LHS, RHS, pre_solution, mode);

            }
            else if constexpr(std::is_same_v<std::decay_t<decltype(arg)>, PMOS>){
                PMOS_assigner(arg.id, arg.node_vd, arg.node_vg, arg.node_vs, arg.node_vb, arg.W, arg.L, h, pre_solution, ckt.T_nodes, LHS, RHS, mode, C_list, NR_iteration_counter);
                
            }

        }, element.element);
    }


    
    return {LHS, RHS};
}

bool isConverge(const std::deque<arma::mat> &NR_solutions, const CKTcircuit &ckt)
{   
    // Validate inputs
    if (NR_solutions.size() < 3) {
        throw std::runtime_error("Not enough solutions in isConverge function for convergence check.");
    }

    arma::mat pre_solution = NR_solutions.front();
    arma::mat current_solution = NR_solutions.at(1);
    arma::mat next_solution = NR_solutions.back();

    if(pre_solution.n_rows != current_solution.n_rows || pre_solution.n_rows != next_solution.n_rows || current_solution.n_rows != next_solution.n_rows){
        std::cerr << "The size of pre_solution, current_solution, next_solution is not the same in isConverge function." << std::endl;
        exit(1);
    }

    arma::mat pre_voltages = pre_solution.submat(0, 0, ckt.T_nodes-1, 0);
    arma::mat current_voltages = current_solution.submat(0, 0, ckt.T_nodes-1, 0);
    arma::mat next_voltages = next_solution.submat(0, 0, ckt.T_nodes-1, 0);

    arma::mat pre_current = pre_solution.submat(ckt.T_nodes, 0, pre_solution.n_rows - 1, 0);
    arma::mat current_current = current_solution.submat(ckt.T_nodes, 0, current_solution.n_rows - 1, 0);
    arma::mat next_current = next_solution.submat(ckt.T_nodes, 0, next_solution.n_rows - 1, 0);


    // |v(k+1)-v(k)| <= RELTOL*max(|v(k+1)|,|v(k)|) + VNTOL

        arma::mat Vmax_next_current = arma::max(arma::abs(next_voltages), arma::abs(current_voltages));
        arma::mat Vdelta_next_current = arma::abs(next_voltages - current_voltages);
        arma::vec tolerance_1 = RELTOL * Vmax_next_current + VNTOL;
        // tolerance_1.print("tolerance_1 =");

        // If |v(k+1)-v(k)| is bigger arma::any will return true and the function will return false.
        if (arma::any(Vdelta_next_current > tolerance_1))
        {
            return false;
        }
    
 
    // |i(k+1)-i(k)| <= RELTOL*max(|i(k+1)|,|i(k)|) + ABSTOL

        arma::mat Imax_next_current = arma::max(arma::abs(next_current), arma::abs(current_current));
        arma::mat Idelta_next_current = arma::abs(next_current - current_current);
        arma::vec tolerance_2 = RELTOL * Imax_next_current + ABSTOL;
        // tolerance_2.print("tolerance_2 =");

        if (arma::any(Idelta_next_current > tolerance_2))
        {
            return false;
        }
    

    // |v(k) - v(k-1)| <= RELTOL*max(|v(k)|,|v(k-1)|) + VNTOL

        arma::mat Vmax_current_pre = arma::max(arma::abs(current_voltages), arma::abs(pre_voltages));
        arma::mat Vdelta_current_pre = arma::abs(current_voltages - pre_voltages);
        arma::vec tolerance_3 = RELTOL * Vmax_current_pre + VNTOL;

        if (arma::any(Vdelta_current_pre > tolerance_3))
        {
            return false;
        }
    

    // |i(k) - i(k-1)| <= RELTOL*max(|i(k)|,|i(k-1)|) + ABSTOL

        arma::mat Imax_current_pre = arma::max(arma::abs(current_current), arma::abs(pre_current));
        arma::mat Idelta_current_pre = arma::abs(current_current - pre_current);
        arma::vec tolerance_4 = RELTOL * Imax_current_pre + ABSTOL;

        if (arma::any(Idelta_current_pre > tolerance_4))
        {
            return false;
        }
    

    // // |v(k+1) - v(k-1)| <= √|v(k) - v(k-1)|^2 + |v(k+1) - v(k)|^2
    
    //     arma::mat Vdelta_next_pre = arma::abs(next_voltages - pre_voltages);
    //     arma::mat V_Perpendicular = arma::sqrt(arma::pow(arma::abs(Vdelta_current_pre),2) + arma::pow(arma::abs(Vdelta_next_current),2));
        
    //     if (arma::any(arma::vectorise(Vdelta_next_pre > V_Perpendicular)) )
    //     {
    //         return false;
    //     }
    

    // // |i(k+1) - i(k-1)| <= √|i(k) - i(k-1)|^2 + |i(k+1) - i(k)|^2

    //     arma::mat Idelta_next_pre = arma::abs(next_current - pre_current);
    //     arma::mat I_Perpendicular = arma::sqrt(arma::pow(arma::abs(Idelta_current_pre),2) + arma::pow(arma::abs(Idelta_next_current),2));

    //     if (arma::any(arma::vectorise(Idelta_next_pre > I_Perpendicular)))
    //     {
    //         return false;
    //     }
    

    return true;
}


// Newton Raphson system solver for non-linear and dynamic elements
arma::mat NewtonRaphson_system(const CKTcircuit &ckt, const arma::mat &pre_solution, const double &h, const int &mode, const double time_trans, std::vector<Capacitor> &C_list){
    
    // std::cout << "Enter NewtonRaphson system" << std::endl;
    DEBUG_PRINT("Enter NewtonRaphson system");

    int NR_iteration_counter = 0;
    bool isconverge = false;
    arma::mat solution = pre_solution;

    std::deque<arma::mat> NR_solutions;
    NR_solutions.push_back(pre_solution);  // save the previous solution

    do{ auto t1 = std::chrono::high_resolution_clock::now();
        std::pair<arma::sp_mat, arma::mat> matrices = DynamicNonLinear(ckt,h, solution, mode, time_trans, C_list, NR_iteration_counter);
        auto t2 = std::chrono::high_resolution_clock::now();
        const arma::sp_mat &LHS = matrices.first;
        const arma::mat &RHS = matrices.second;

        // double rcond_LHS = rcond(LHS);

        // Solve Ax = b
        // J(v) * x(k+1) = [J(v)]x(k) - f(x(k))
        auto t3 = std::chrono::high_resolution_clock::now();
        // solution = arma::solve(LHS, RHS);
        solution = arma::spsolve(LHS, RHS, "superlu");
        auto t4 = std::chrono::high_resolution_clock::now();
        NR_iteration_counter += 1;


        //update the NR_solutions(current_solution, pre_solution, next_solution) 
        NR_solutions.push_back(solution);
        if (NR_solutions.size() > 3)
        {
            NR_solutions.pop_front();
        }

        // Check for convergence only if you have 3 solutions
        if (NR_solutions.size() == 3) {
            isconverge = isConverge(NR_solutions, ckt);
        }

        if (NR_iteration_counter > 100)
        {
            std::cerr << "Not Converge at 100 iterations" << std::endl;
            exit(1);
            break;
        }
        // NR_ITE = NR_iteration_counter; // Global variable for the number of NR iterations
        auto t5 = std::chrono::high_resolution_clock::now();

        auto duration1 = std::chrono::duration<double, std::milli>( t2 - t1 );
        auto duration2 = std::chrono::duration<double, std::milli>( t4 - t3 );
        auto duration3 = std::chrono::duration<double, std::milli>( t5 - t1 );
        std::cout << "Time for DynamicNonLinear: " << duration1.count() << " ms" << std::endl;
        std::cout << "Time for spsolve: " << duration2.count() << " ms" << std::endl;
        std::cout << "Total Time: " << duration3.count() << " ms" << std::endl;
    }while(!isconverge);


    return solution;
}

// Overloading the NewtonRaphson_system function for the multiple time steps
// This version is used in multi_solution_solver
// arma::mat NewtonRaphson_system(const CKTcircuit &ckt, const arma::mat &pre_solution, const double &h, const int &mode, const double time_trans, 
//     std::vector<Capacitor> &C_list, bool &ITE_reach){
    
//     // std::cout << "Enter NewtonRaphson system" << std::endl;
//     DEBUG_PRINT("Enter NewtonRaphson system");

//     auto startNR = std::chrono::high_resolution_clock::now();

//     int NR_iteration_counter = 0;
//     bool isconverge = false;
//     arma::mat solution = pre_solution;

//     std::deque<arma::mat> NR_solutions;
//     NR_solutions.push_back(pre_solution);  // save the previous solution

//     do{
//         NR_iteration_counter += 1;

//         if(NR_iteration_counter > 100){
//             ITE_reach = true;
//             break;
//         }

//         std::pair<arma::mat, arma::mat> matrices = DynamicNonLinear(ckt,h, solution, mode, time_trans, C_list, NR_iteration_counter);
//         const arma::mat &LHS = matrices.first;
//         const arma::mat &RHS = matrices.second;

//         // Solve Ax = b
//         // J(v) * x(k+1) = [J(v)]x(k) - f(x(k))

//         solution = arma::solve(LHS, RHS);

//         //update the NR_solutions(current_solution, pre_solution, next_solution) 
//         NR_solutions.push_back(solution);
//         if (NR_solutions.size() > 3)
//         {
//             NR_solutions.pop_front();
//         }

//         // Check for convergence only if you have 3 solutions
//         if (NR_solutions.size() == 3) {
//             isconverge = isConverge(NR_solutions, ckt);
//         }

//     }while(!isconverge);

//     return solution;
// }

std::pair<arma::mat,arma::mat> get_currents_voltages(const std::vector<Capacitor> &pre_C_list, const double h, const arma::mat &solution , const arma::mat &pre_solution){

    // h/C * i = u(k+1) - u(k)

    int G_rows = pre_C_list.size();
    int G_cols = G_rows;
    double vol{}, pre_vol{}, delta_vol{};               // Delta voltage across the capacitor u(k+1) - u(k)

    arma::mat current_matrix = arma::zeros(G_rows, 1);  // Currents matrix
    arma::mat delta_v = arma::zeros(G_rows, 1);         // Delta voltage matrix
    arma::vec G_vec(pre_C_list.size());                 //Initialize arma::vec with the size of C_list
    arma::mat volt = arma::zeros(G_rows, 1);

    for(size_t i = 0; i < pre_C_list.size(); ++i){

        G_vec(i) = h / pre_C_list.at(i).value;

        if(pre_C_list.at(i).nodePos == 0){

            vol = solution(pre_C_list.at(i).nodeNeg - 1, 0);
            pre_vol = pre_solution(pre_C_list.at(i).nodeNeg - 1, 0);
            delta_vol = vol - pre_vol;          //u(k+1) - u(k)
        }
        else if(pre_C_list.at(i).nodeNeg == 0){

            vol = solution(pre_C_list.at(i).nodePos - 1, 0);
            pre_vol = pre_solution(pre_C_list.at(i).nodePos - 1, 0);
            delta_vol = vol - pre_vol;          //u(k+1) - u(k)
        }
        else{

            vol = solution(pre_C_list.at(i).nodePos - 1, 0) - solution(pre_C_list.at(i).nodeNeg - 1, 0);
            pre_vol = pre_solution(pre_C_list.at(i).nodePos - 1, 0) - pre_solution(pre_C_list.at(i).nodeNeg - 1, 0);
            delta_vol = vol - pre_vol;          //u(k+1) - u(k)
        }

        delta_v(i,0) = delta_vol;              // It may be negative!
        volt(i,0) = vol;
        // std::cout << "c: " << pre_C_list.at(i).value << std::endl;
        // std::cout << "C/h: " << pre_C_list.at(i).value/h << std::endl;
        // std::cout << "delta_vol: " << delta_vol << std::endl;
        // std::cout << "C/h * delta_vol" << pre_C_list.at(i).value/h * delta_vol << std::endl;
    }

    arma::mat G_matrix = arma::diagmat(G_vec);  // Diagonal matrix
    // G_matrix.print("G_matrix =");
    // delta_v.print("delta_v =");

    ARMA_PRINT(G_matrix, "G_matrix = ");
    ARMA_PRINT(delta_v, "delta_v = ");
    

    current_matrix = arma::solve(G_matrix, delta_v);

    return {current_matrix, volt};
}


// New timestep options
void timestep_options(double & temp_h, double & next_h_up, double & next_h_down, const Transient &trans, bool & TMAX_reach, bool & TMIN_reach) {
    if(temp_h >= trans.h_MAX) {
        temp_h = trans.h_MAX;
        TMAX_reach = true;
        
    }
    if(temp_h <= trans.h_MIN) {
        temp_h = trans.h_MIN;
        TMIN_reach = true;
        
    }

    next_h_up = temp_h * 2;
    next_h_down = temp_h / 2;

    if(next_h_up > trans.h_MAX) {
        next_h_up = trans.h_MAX;
        // TMAX_reach = true;
        
    }
    if(next_h_up < trans.h_MIN) {
        next_h_down = trans.h_MIN;
        // TMAX_reach = true;
        
    }
    if(next_h_down > trans.h_MAX) {
        next_h_up = trans.h_MAX;
        // TMAX_reach = true;
        
    }
    if(next_h_down < trans.h_MIN) {
        next_h_down = trans.h_MIN;
        // TMIN_reach = true;
        
    }

}

int h_reach_new(const Truncation_error LTE, int &last_decision, const multi_timestep &multi_h){

        /* index_h:
            0 = temp_h / 4
            1 = down
            2 = mid
            3 = up
        */ 
    int reach{};

    // if(multi_h.TMAX_reach){
    //     bool mid = multi_h.h_mid <= TRTOL * LTE.LTE_BE_mid_h.min();
    //     bool down = multi_h.h_down <= TRTOL * LTE.LTE_BE_down_h.min();

    //     if(down == false){
    //         last_decision = 2; 
    //         reach = 0;          // temp_h / 4
    //         return reach;
    //     }

    //     if(mid == true && down == true){

    //         reach = 2;          // mid
    //         return reach;
    //     }

    //     if(mid == false && down == false){

    //         last_decision = 2;
    //         reach = 0;          // temp_h / 4
    //         return reach;

    //     }

    //     if(mid == false && down == true){
    //         last_decision = 0;
    //         reach = 1;          // down
    //         return reach;
    //     }
        
    //     std::cerr << "h_reach_new function error" << std::endl;
    //     exit(1);
    
    // }

    // else if(multi_h.TMIN_reach){
    //     bool mid = multi_h.h_mid <= TRTOL * LTE.LTE_BE_mid_h.min();
    //     bool up = multi_h.h_up <= TRTOL * LTE.LTE_BE_up_h.min();
   
    //     if(mid == false){
    //         last_decision = 2;   // temp_h / 4
    //         reach = 0;
    //         return reach;
    //     }

    //     if(mid == true && up == true){

    //         reach = 3;          // up
    //         return reach;
    //     }

    //     if(mid == false && up == false){

    //         std::cerr << "Time step is too small (Reach minimum!)" << std::endl;
    //         exit(1);

    //     }

    //     if(mid == true && up == false){
    //         last_decision = 0;
    //         reach = 2;          // mid
    //         return reach;
    //     }
        
    //     std::cerr << "h_reach_new function error" << std::endl;
    //     exit(1);
    // }

    // else{
        // bool mid;
        // bool up;
        // bool down;

        // if(multi_h.ITE_reach_mid){
        //     mid = false;
        // }
        // else{
        //     mid = multi_h.h_mid <= TRTOL * LTE.LTE_BE_mid_h.min();
        // }
        // if(multi_h.ITE_reach_up){
        //     up = false;
        // }
        // else{
        //     up = multi_h.h_up <= TRTOL * LTE.LTE_BE_up_h.min();
        // }
        // if(multi_h.ITE_reach_down){
        //     down = false;
        // }
        // else{
        //     down = multi_h.h_down <= TRTOL * LTE.LTE_BE_down_h.min();
        // }
        bool mid = multi_h.h_mid <= TRTOL * LTE.LTE_BE_mid_h.min();
        bool up = multi_h.h_up <= TRTOL * LTE.LTE_BE_up_h.min();
        bool down = multi_h.h_down <= TRTOL * LTE.LTE_BE_down_h.min();

        
        if(down == false){
            last_decision = 2; // temp_h / 8
            reach = 0;
            return reach;
        }

        if(mid == true && up == true && down == true){

            reach = 3;  // up
            return reach;
        }

        if(mid == false && up == false && down == false){

            last_decision = 2;
            reach = 0;  // temp_h / 8
            return reach;

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
        
        std::cerr << "h_reach_new function error" << std::endl;
        exit(1);
    // }
    

}



multi_timestep multi_solution_solver(const double &h, const Transient &trans, const CKTcircuit &ckt){
    
    multi_timestep multi_h;

    multi_h.TMAX_reach = false;
    multi_h.TMIN_reach = false;
    multi_h.ITE_reach_mid = false;
    multi_h.ITE_reach_up = false;
    multi_h.ITE_reach_down = false;

    multi_h.h_mid = h;

    timestep_options(multi_h.h_mid, multi_h.h_up, multi_h.h_down, trans, multi_h.TMAX_reach, multi_h.TMIN_reach);  // calculate the next time step, up mid down

    multi_h.t_mid = trans.time_trans + multi_h.h_mid;
    multi_h.t_up = trans.time_trans + multi_h.h_up;
    multi_h.t_down = trans.time_trans + multi_h.h_down;

    // Double check the time steps
    // if(multi_h.t_mid == multi_h.t_up){
    //     if(multi_h.TMAX_reach == false){
    //         std::cerr << "The TMAX_reach is false in multi_solution_solver function" << std::endl;
    //         exit(1);
    //     }
    // }

    // if(multi_h.t_mid == multi_h.t_down){
    //     if(multi_h.TMIN_reach == false){
    //         std::cerr << "The TMIN_reach is false in multi_solution_solver function" << std::endl;
    //         exit(1);
    //     }
    // }

    // std::cout << "The time steps are:" << multi_h.h_mid << " " << multi_h.h_up << " " << multi_h.h_down << std::endl;
    DEBUG_PRINT("The time steps are: " << multi_h.h_mid << " " << multi_h.h_up << " " << multi_h.h_down);
    
    // Mid time step
    pool.detach_task(
        [&trans, &ckt, &multi_h]
        {   
            begin_mid.push_back(std::chrono::high_resolution_clock::now());
            std::vector<Capacitor> C_list_mid = trans.C_list;
            arma::mat mid_solution = NewtonRaphson_system(ckt, vec_trans.back().solution, multi_h.h_mid, 1, multi_h.t_mid, C_list_mid);
            ARMA_PRINT(mid_solution, "mid_solution =");

                multi_h.solution_mid = mid_solution;
                std::pair<arma::mat, arma::mat> mid_currents_voltages = get_currents_voltages(C_list_mid, multi_h.h_mid, mid_solution, vec_trans.back().solution);
                multi_h.C_current_mid = mid_currents_voltages.first;
                multi_h.C_voltage_mid = mid_currents_voltages.second;
                arma::mat mid_Capacitance = get_capacitance(C_list_mid);       // Update the capacitance values because the MOSFETs
                multi_h.C_charge_mid = mid_Capacitance % multi_h.C_voltage_mid;  // element-wise multiplication of two objects (Schur product)
                multi_h.Capacitance_mid = mid_Capacitance;
                multi_h.C_list_mid = C_list_mid;

            end_mid.push_back(std::chrono::high_resolution_clock::now());
            
            // sync_out.print("Time mid_solution_solver = ", tmr.current_ms(),"\n");
        }
    );
    
    // Up time step
    // if(multi_h.TMAX_reach == false){
        pool.detach_task(
            [&trans, &ckt, &multi_h]
            {   
                begin_up.push_back(std::chrono::high_resolution_clock::now());
                std::vector<Capacitor> C_list_up = trans.C_list;
                arma::mat up_solution = NewtonRaphson_system(ckt, vec_trans.back().solution, multi_h.h_up, 1, multi_h.t_up, C_list_up);
                // up_solution.print("up_solution =");
                ARMA_PRINT(up_solution, "up_solution =");

                    multi_h.solution_up = up_solution;
                    std::pair<arma::mat, arma::mat> up_currents_voltages = get_currents_voltages(C_list_up, multi_h.h_up, up_solution, vec_trans.back().solution);
                    multi_h.C_current_up = up_currents_voltages.first;
                    multi_h.C_voltage_up = up_currents_voltages.second;
                    arma::mat up_Capacitance = get_capacitance(C_list_up);       // Update the capacitance values because the MOSFETs
                    multi_h.C_charge_up = up_Capacitance % multi_h.C_voltage_up;  // element-wise multiplication of two objects (Schur product)
                    multi_h.Capacitance_up = up_Capacitance;
                    multi_h.C_list_up = C_list_up;

                end_up.push_back(std::chrono::high_resolution_clock::now());
                // sync_out.print("Time up_solution_solver = ", tmr.current_ms(),"\n");
            }
        );
    // }

    // Down time step
    // if(multi_h.TMIN_reach == false){
        pool.detach_task(
            [&trans, &ckt, &multi_h]
            {   
                begin_down.push_back(std::chrono::high_resolution_clock::now());
                std::vector<Capacitor> C_list_down = trans.C_list;
                arma::mat down_solution = NewtonRaphson_system(ckt, vec_trans.back().solution, multi_h.h_down, 1, multi_h.t_down, C_list_down);
                // down_solution.print("down_solution =");
                ARMA_PRINT(down_solution, "down_solution =");


                    multi_h.solution_down = down_solution;
                    std::pair<arma::mat, arma::mat> down_currents_voltages = get_currents_voltages(C_list_down, multi_h.h_down, down_solution, vec_trans.back().solution);
                    multi_h.C_current_down = down_currents_voltages.first;
                    multi_h.C_voltage_down = down_currents_voltages.second;
                    arma::mat down_Capacitance = get_capacitance(C_list_down);       // Update the capacitance values because the MOSFETs
                    multi_h.C_charge_down = down_Capacitance % multi_h.C_voltage_down;  // element-wise multiplication of two objects (Schur product)
                    multi_h.Capacitance_down = down_Capacitance;
                    multi_h.C_list_down = C_list_down;
                
                end_down.push_back(std::chrono::high_resolution_clock::now());
            }
        );
    // }

    pool.wait();

    if(multi_h.h_mid == trans.h_MIN && multi_h.ITE_reach_mid && multi_h.ITE_reach_up && multi_h.ITE_reach_down){
        std::cerr << "The time step is too small (Reach minimum!)" << std::endl;
        exit(1);
    }
            

    // std::vector<Capacitor> C_list_mid = trans.C_list;
    // arma::mat mid_solution = NewtonRaphson_system(ckt, vec_trans.back().solution, multi_h.h_mid, 1, multi_h.t_mid, C_list_mid);
    // // mid_solution.print("mid_solution =");
    // ARMA_PRINT(mid_solution, "mid_solution =");
    // multi_h.solution_mid = mid_solution;
    // auto mid_currents_voltages = get_currents_voltages(C_list_mid, multi_h.h_mid, mid_solution, vec_trans.back().solution);
    // multi_h.C_current_mid = mid_currents_voltages.first;
    // multi_h.C_voltage_mid = mid_currents_voltages.second;
    // arma::mat mid_Capacitance = get_capacitance(C_list_mid);       // Update the capacitance values because the MOSFETs
    // multi_h.C_charge_mid = mid_Capacitance % multi_h.C_voltage_mid;  // element-wise multiplication of two objects (Schur product)
    // // multi_h.LHS_mid = mid_matrixs.first;
    // // multi_h.RHS_mid = mid_matrixs.second;
    // multi_h.Capacitance_mid = mid_Capacitance;
    // multi_h.C_list_mid = C_list_mid;


    // std::vector<Capacitor> C_list_up = trans.C_list;
    // arma::mat up_solution = NewtonRaphson_system(ckt, vec_trans.back().solution, multi_h.h_up, 1, multi_h.t_up, C_list_up);
    // // up_solution.print("up_solution =");
    // ARMA_PRINT(up_solution, "up_solution =");
    // multi_h.solution_up = up_solution;
    // auto up_currents_voltages = get_currents_voltages(C_list_up, multi_h.h_up, up_solution, vec_trans.back().solution);
    // multi_h.C_current_up = up_currents_voltages.first;
    // multi_h.C_voltage_up = up_currents_voltages.second;
    // arma::mat up_Capacitance = get_capacitance(C_list_up);       // Update the capacitance values because the MOSFETs
    // multi_h.C_charge_up = up_Capacitance % multi_h.C_voltage_up;  // element-wise multiplication of two objects (Schur product)
    // // multi_h.LHS_up = up_matrixs.first;
    // // multi_h.RHS_up = up_matrixs.second;
    // multi_h.Capacitance_up = up_Capacitance;
    // multi_h.C_list_up = C_list_up;


    // std::vector<Capacitor> C_list_down = trans.C_list;
    // arma::mat down_solution = NewtonRaphson_system(ckt, vec_trans.back().solution, multi_h.h_down, 1, multi_h.t_down, C_list_down);
    // // down_solution.print("down_solution =");
    // ARMA_PRINT(down_solution, "down_solution =");
    // multi_h.solution_down = down_solution;
    // auto down_currents_voltages = get_currents_voltages(C_list_down, multi_h.h_down, down_solution, vec_trans.back().solution);
    // multi_h.C_current_down = down_currents_voltages.first;
    // multi_h.C_voltage_down = down_currents_voltages.second;
    // arma::mat down_Capacitance = get_capacitance(C_list_down);       // Update the capacitance values because the MOSFETs
    // multi_h.C_charge_down = down_Capacitance % multi_h.C_voltage_down;  // element-wise multiplication of two objects (Schur product)
    // // multi_h.LHS_down = down_matrixs.first;
    // // multi_h.RHS_down = down_matrixs.second;
    // multi_h.Capacitance_down = down_Capacitance;
    // multi_h.C_list_down = C_list_down;


    return multi_h;
}

arma::mat multi_next_h(Transient &trans, const CKTcircuit &ckt){

    // std::cout << "Enter multi_next_h" << std::endl;
    DEBUG_PRINT("Enter multi_next_h");

    arma::mat solution;

    double temp_h = vec_trans.back().h;
    double lats_step = vec_trans.back().h;
    arma::mat pre_current = vec_trans.back().C_current;
    arma::mat pre_charge = vec_trans.back().C_charge;

    int index_h{};
    int pre_decision{};
    multi_timestep multi_h;
    multi_timestep pre_multi_h;
    Truncation_error LTE;

    do{
        timer.start();

        multi_h = multi_solution_solver(temp_h, trans, ckt);

        timer.stop();
        timer.total();

        ARMA_PRINT(multi_h.Capacitance_mid, "Capacitance_mid = ");
        ARMA_PRINT(multi_h.Capacitance_up, "Capacitance_up = ");
        ARMA_PRINT(multi_h.Capacitance_down, "Capacitance_down = ");

        ARMA_PRINT(multi_h.C_current_mid, "C_current_mid = ");
        ARMA_PRINT(multi_h.C_current_up, "C_current_up = ");
        ARMA_PRINT(multi_h.C_current_down, "C_current_down = ");

        ARMA_PRINT(multi_h.C_charge_mid, "C_charge_mid = ");
        ARMA_PRINT(multi_h.C_charge_up, "C_charge_up = ");
        ARMA_PRINT(multi_h.C_charge_down, "C_charge_down = ");

        // if(multi_h.TMAX_reach){
        //     arma::mat mid_max_current = arma::max(arma::abs(multi_h.C_current_mid), arma::abs(pre_current));
        //     arma::mat down_max_current = arma::max(arma::abs(multi_h.C_current_down), arma::abs(pre_current));

        //     arma::mat mid_max_charge = arma::max(arma::abs(multi_h.C_charge_mid), arma::abs(pre_charge));
        //     arma::mat down_max_charge = arma::max(arma::abs(multi_h.C_charge_down), arma::abs(pre_charge));

        //     // LTE current = ABSTOL + RELTOL * max(|ik+1|,|ik|) in SPICE book
        //     arma::mat mid_LTE_current = ABSTOL + RELTOL * mid_max_current;
        //     arma::mat down_LTE_current = ABSTOL + RELTOL * down_max_current;

        //     // LTE charge = RELTOL * max(|qk+1|,|qk|, chgtol) / h(k) in SPICE book
        //     arma::mat CHGTOL_MAT = arma::mat(mid_max_charge.n_rows, mid_max_charge.n_cols).fill(CHGTOL);
        //     arma::mat mid_LTE_charge = RELTOL * arma::max(mid_max_charge, CHGTOL_MAT) / multi_h.h_mid;
        //     arma::mat down_LTE_charge = RELTOL * arma::max(down_max_charge, CHGTOL_MAT) / multi_h.h_down;

        //     // LTE bound = max(LTE_current, LTE_charge)
        //     LTE.LTE_bound_mid = arma::max(mid_LTE_current, mid_LTE_charge);
        //     LTE.LTE_bound_down = arma::max(down_LTE_current, down_LTE_charge);

        //     arma::mat BEcur_diff_mid = multi_h.C_current_mid - pre_current;
        //     arma::mat BEcur_diff_down = multi_h.C_current_down - pre_current;

        //     // tol_h = 2C / |i(k+1) - i(k)| * LTE_bound 
        //     LTE.LTE_BE_mid_h = (2 * multi_h.Capacitance_mid) / arma::abs(BEcur_diff_mid) % LTE.LTE_bound_mid;
        //     LTE.LTE_BE_down_h = (2 * multi_h.Capacitance_down) / arma::abs(BEcur_diff_down) % LTE.LTE_bound_down;
        // }

        // else if(multi_h.TMIN_reach){
        //     arma::mat mid_max_current = arma::max(arma::abs(multi_h.C_current_mid), arma::abs(pre_current));
        //     arma::mat up_max_current = arma::max(arma::abs(multi_h.C_current_up), arma::abs(pre_current));

        //     arma::mat mid_max_charge = arma::max(arma::abs(multi_h.C_charge_mid), arma::abs(pre_charge));
        //     arma::mat up_max_charge = arma::max(arma::abs(multi_h.C_charge_up), arma::abs(pre_charge));

        //     // LTE current = ABSTOL + RELTOL * max(|ik+1|,|ik|) in SPICE book
        //     arma::mat mid_LTE_current = ABSTOL + RELTOL * mid_max_current;
        //     arma::mat up_LTE_current = ABSTOL + RELTOL * up_max_current;

        //     // LTE charge = RELTOL * max(|qk+1|,|qk|, chgtol) / h(k) in SPICE book
        //     arma::mat CHGTOL_MAT = arma::mat(mid_max_charge.n_rows, mid_max_charge.n_cols).fill(CHGTOL);
        //     arma::mat mid_LTE_charge = RELTOL * arma::max(mid_max_charge, CHGTOL_MAT) / multi_h.h_mid;
        //     arma::mat up_LTE_charge = RELTOL * arma::max(up_max_charge, CHGTOL_MAT) / multi_h.h_up;

        //     // LTE bound = max(LTE_current, LTE_charge)
        //     LTE.LTE_bound_mid = arma::max(mid_LTE_current, mid_LTE_charge);
        //     LTE.LTE_bound_up = arma::max(up_LTE_current, up_LTE_charge);

        //     arma::mat BEcur_diff_mid = multi_h.C_current_mid - pre_current;
        //     arma::mat BEcur_diff_up = multi_h.C_current_up - pre_current;

        //     // tol_h = 2C / |i(k+1) - i(k)| * LTE_bound 
        //     LTE.LTE_BE_mid_h = (2 * multi_h.Capacitance_mid) / arma::abs(BEcur_diff_mid) % LTE.LTE_bound_mid;
        //     LTE.LTE_BE_up_h = (2 * multi_h.Capacitance_up) / arma::abs(BEcur_diff_up) % LTE.LTE_bound_up;
        // }

        // else{
        // if(multi_h.ITE_reach_mid == false){
        //     arma::mat mid_max_current = arma::max(arma::abs(multi_h.C_current_mid), arma::abs(pre_current));
        //     arma::mat mid_max_charge = arma::max(arma::abs(multi_h.C_charge_mid), arma::abs(pre_charge));
        //     arma::mat mid_LTE_current = ABSTOL + RELTOL * mid_max_current;

        //     arma::mat CHGTOL_MAT = arma::mat(mid_max_charge.n_rows, mid_max_charge.n_cols).fill(CHGTOL);
        //     arma::mat mid_LTE_charge = RELTOL * arma::max(mid_max_charge, CHGTOL_MAT) / multi_h.h_mid;
        //     LTE.LTE_bound_mid = arma::max(mid_LTE_current, mid_LTE_charge);
        //     arma::mat BEcur_diff_mid = multi_h.C_current_mid - pre_current;
        //     LTE.LTE_BE_mid_h = (2 * multi_h.Capacitance_mid) / arma::abs(BEcur_diff_mid) % LTE.LTE_bound_mid;
        // }
        // if(multi_h.ITE_reach_up == false){
        //     arma::mat up_max_current = arma::max(arma::abs(multi_h.C_current_up), arma::abs(pre_current));
        //     arma::mat up_max_charge = arma::max(arma::abs(multi_h.C_charge_up), arma::abs(pre_charge));
        //     arma::mat up_LTE_current = ABSTOL + RELTOL * up_max_current;
        //     arma::mat CHGTOL_MAT = arma::mat(up_max_charge.n_rows, up_max_charge.n_cols).fill(CHGTOL);
        //     arma::mat up_LTE_charge = RELTOL * arma::max(up_max_charge, CHGTOL_MAT) / multi_h.h_up;
        //     LTE.LTE_bound_up = arma::max(up_LTE_current, up_LTE_charge);
        //     arma::mat BEcur_diff_up = multi_h.C_current_up - pre_current;
        //     LTE.LTE_BE_up_h = (2 * multi_h.Capacitance_up) / arma::abs(BEcur_diff_up) % LTE.LTE_bound_up;
        // }
        // if(multi_h.ITE_reach_down == false){
        //     arma::mat down_max_current = arma::max(arma::abs(multi_h.C_current_down), arma::abs(pre_current));
        //     arma::mat down_max_charge = arma::max(arma::abs(multi_h.C_charge_down), arma::abs(pre_charge));
        //     arma::mat down_LTE_current = ABSTOL + RELTOL * down_max_current;
        //     arma::mat CHGTOL_MAT = arma::mat(down_max_charge.n_rows, down_max_charge.n_cols).fill(CHGTOL);
        //     arma::mat down_LTE_charge = RELTOL * arma::max(down_max_charge, CHGTOL_MAT) / multi_h.h_down;
        //     LTE.LTE_bound_down = arma::max(down_LTE_current, down_LTE_charge);
        //     arma::mat BEcur_diff_down = multi_h.C_current_down - pre_current;
        //     LTE.LTE_BE_down_h = (2 * multi_h.Capacitance_down) / arma::abs(BEcur_diff_down) % LTE.LTE_bound_down;
        // }

            arma::mat mid_max_current = arma::max(arma::abs(multi_h.C_current_mid), arma::abs(pre_current));
            arma::mat up_max_current = arma::max(arma::abs(multi_h.C_current_up), arma::abs(pre_current));
            arma::mat down_max_current = arma::max(arma::abs(multi_h.C_current_down), arma::abs(pre_current));

            arma::mat mid_max_charge = arma::max(arma::abs(multi_h.C_charge_mid), arma::abs(pre_charge));
            arma::mat up_max_charge = arma::max(arma::abs(multi_h.C_charge_up), arma::abs(pre_charge));
            arma::mat down_max_charge = arma::max(arma::abs(multi_h.C_charge_down), arma::abs(pre_charge));

            // LTE current = ABSTOL + RELTOL * max(|ik+1|,|ik|) in SPICE book
            arma::mat mid_LTE_current = ABSTOL + RELTOL * mid_max_current;
            arma::mat up_LTE_current = ABSTOL + RELTOL * up_max_current;
            arma::mat down_LTE_current = ABSTOL + RELTOL * down_max_current;

            // LTE charge = RELTOL * max(|qk+1|,|qk|, chgtol) / h(k) in SPICE book
            arma::mat CHGTOL_MAT = arma::mat(mid_max_charge.n_rows, mid_max_charge.n_cols).fill(CHGTOL);
            arma::mat mid_LTE_charge = RELTOL * arma::max(mid_max_charge, CHGTOL_MAT) / multi_h.h_mid;
            arma::mat up_LTE_charge = RELTOL * arma::max(up_max_charge, CHGTOL_MAT) / multi_h.h_up;
            arma::mat down_LTE_charge = RELTOL * arma::max(down_max_charge, CHGTOL_MAT) / multi_h.h_down;

            // LTE bound = max(LTE_current, LTE_charge)
            LTE.LTE_bound_mid = arma::max(mid_LTE_current, mid_LTE_charge);
            LTE.LTE_bound_up = arma::max(up_LTE_current, up_LTE_charge);
            LTE.LTE_bound_down = arma::max(down_LTE_current, down_LTE_charge);

            arma::mat BEcur_diff_mid = multi_h.C_current_mid - pre_current;
            arma::mat BEcur_diff_up = multi_h.C_current_up - pre_current;
            arma::mat BEcur_diff_down = multi_h.C_current_down - pre_current;

            // tol_h = 2C / |i(k+1) - i(k)| * LTE_bound 
            LTE.LTE_BE_mid_h = (2 * multi_h.Capacitance_mid) / arma::abs(BEcur_diff_mid) % LTE.LTE_bound_mid;
            LTE.LTE_BE_up_h = (2 * multi_h.Capacitance_up) / arma::abs(BEcur_diff_up) % LTE.LTE_bound_up;
            LTE.LTE_BE_down_h = (2 * multi_h.Capacitance_down) / arma::abs(BEcur_diff_down) % LTE.LTE_bound_down;
        // }


        index_h = h_reach_new(LTE, pre_decision, multi_h);

        /* index_h:
            0 = temp_h / 4
            1 = down
            2 = mid
            3 = up
        */ 
       switch(index_h){
        case 0:
            temp_h = temp_h / 4.0;
            pre_multi_h = multi_h;
            break;
        case 1:
            trans.h = multi_h.h_down;
            trans.LHS = multi_h.LHS_down;
            trans.RHS = multi_h.RHS_down;
            trans.solution = multi_h.solution_down;
            trans.C_list = multi_h.C_list_down;
            trans.Capacitance = multi_h.Capacitance_down;
            solution = multi_h.solution_down;

            trans.C_current = multi_h.C_current_down;
            trans.C_voltage = multi_h.C_voltage_down;
            trans.C_charge = multi_h.C_charge_down;
            break;
        case 2:
            trans.h = multi_h.h_mid;
            trans.LHS = multi_h.LHS_mid;
            trans.RHS = multi_h.RHS_mid;
            trans.solution = multi_h.solution_mid;
            trans.C_list = multi_h.C_list_mid;
            trans.Capacitance = multi_h.Capacitance_mid;
            solution = multi_h.solution_mid;

            trans.C_current = multi_h.C_current_mid;
            trans.C_voltage = multi_h.C_voltage_mid;
            trans.C_charge = multi_h.C_charge_mid;
            break;
        case 3:
            trans.h = multi_h.h_up;
            trans.LHS = multi_h.LHS_up;
            trans.RHS = multi_h.RHS_up;
            trans.solution = multi_h.solution_up;
            trans.C_list = multi_h.C_list_up;
            trans.Capacitance = multi_h.Capacitance_up;
            solution = multi_h.solution_up;

            trans.C_current = multi_h.C_current_up;
            trans.C_voltage = multi_h.C_voltage_up;
            trans.C_charge = multi_h.C_charge_up;
            break;

        default:
            std::cerr << "Error in multi_next_h function" << std::endl;
            exit(1);
       }


    }while(index_h == 0);


    return solution;
}

void save_csv(const CKTcircuit &ckt){
    std::ofstream file("final_solution.csv");

    // Write the header
    file << "Time, Time Step";
    for (size_t j = 0; j < ckt.external_nodes; ++j) {
        file << ", Voltage " << (j + 1);
    }
    for (size_t z = ckt.no_of_V_sources; z > 0; --z) {
        file << ", Current " << (ckt.no_of_V_sources - z + 1);
    }
    // file << ", Current 2";
    file << std::endl;

    // Write the data
    for(size_t i = 0; i < vec_trans.size(); ++i){
        file << std::scientific << std::setprecision(20); // Set precision to 20 decimal places
        file << vec_trans.at(i).time_trans << ", " << vec_trans.at(i).h ; // Time and Timestep

        for(size_t j = 0; j < ckt.external_nodes; ++j){
            file << ", " << vec_trans.at(i).solution(j,0); // Voltages
        }

        for(size_t z = ckt.no_of_V_sources; z > 0; --z){
            file << ", " << vec_trans.at(i).solution(ckt.cktdematrix-> n_rows - z,0); // Currents
        }
        // file << ", " << vec_trans.at(i).C_current(0,0);
        file << std::endl;
    }

    file.close();

}

void save_threads_time(const std::chrono::time_point<std::chrono::high_resolution_clock> &t1, const std::chrono::time_point<std::chrono::high_resolution_clock> &t2){
    std::ofstream file("Threads time.csv");
    if(begin_mid.size() != begin_down.size() || begin_mid.size() != begin_up.size() || begin_down.size() != begin_up.size()){
        std::cerr << "The size of the vectors are not equal" << std::endl;
        exit(1);
    }
    if(end_mid.size() != end_down.size() || end_mid.size() != end_up.size() || end_down.size() != end_up.size()){
        std::cerr << "The size of the vectors are not equal" << std::endl;
        exit(1);
    }

    std::vector<std::chrono::duration<double, std::milli>> bg_mid(begin_mid.size()), bg_up(begin_mid.size()), bg_down(begin_mid.size()), 
        ed_mid(begin_mid.size()), ed_up(begin_mid.size()), ed_down(begin_mid.size());
    
    for(int i = 0; i < begin_mid.size(); ++i){
        bg_mid.at(i) = begin_mid.at(i) - t1;
        bg_up.at(i) = begin_up.at(i) - t1;
        bg_down.at(i) = begin_down.at(i) - t1;

        ed_mid.at(i) = end_mid.at(i) - t1;
        ed_up.at(i) = end_up.at(i) - t1;
        ed_down.at(i) = end_down.at(i) - t1;
    }

    file << "Begin Mid, Begin Up, Begin Down, End Mid, End Up, End Down" << std::endl;
    for(int i = 0; i < begin_mid.size(); ++i){
        file << bg_mid.at(i).count() << ", " << bg_up.at(i).count() << ", " << bg_down.at(i).count() << ", " 
        << ed_mid.at(i).count() << ", " << ed_up.at(i).count() << ", " << ed_down.at(i).count() << std::endl;
    }

    std::ofstream file2("simulation end time.csv");
    file2 << "Simulation end time" << std::endl;
    std::chrono::duration<double, std::milli> total_time = (t2 - t1);
    file2 << total_time.count() << std::endl;
}

void dummy_task()
{
    // std::cout << "dummy task" << std::endl;
    std::this_thread::sleep_for(std::chrono::microseconds(100));
}

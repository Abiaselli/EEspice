#define ARMA_DONT_USE_WRAPPER
#define ARMA_USE_MKL_ALLOC
#define ARMA_USE_SUPERLU

#include <omp.h>
#include "model_setup.hpp"
#include "Transient.hpp"
#include "DC.hpp"
#include "AC.hpp"
#include "saveCSV.hpp"
#include "batch.hpp"
#include "simulation_exceptions.hpp"

// Main function for the circuit simulation
int main(int argc, const char **argv)
{
    // Check if a netlist file is provided as command line argument
    if (argc < 2) {
        std::cerr << "Usage: ./eespice <netlist_file>" << std::endl;
        return 1;
    }

    try {
        XB_Timer total;
        total.start();
        setDebugMode(false); // Set the debug mode to false or true

        // Parse netlist file
        Modelmap modmap;
        Circuitmap cktmap;
        CircuitParser parser(argv[1]);
        parser_netlist(parser, cktmap, modmap);

        // Model setup using the temperature
        modelSetup(modmap, nomTemp);

        if(parser.is_batch){
            std::cout << "Starting batch simulation..." << std::endl;
            batchMode = true; // Set the global batch mode to true
            auto batch_results = batch::run_batch_simulation(cktmap, parser, modmap);
            std::cout << "Batch simulation finished. Saving " << batch_results.size() << " results." << std::endl;
            batch::save_csv_batch(batch_results);
        }
        else{           
            // CKT circuit setup
            CKTcircuit ckt;
            {
                ScopedTimer setupTimer(ckt.sim_stats.simTime.setup_time); // Time the CKT setup
                ckt.map = cktmap; // Assign the circuit map to the CKTcircuit
                auto MatrixPtr = std::make_shared<Matrix>();  // Create the Matrix as a shared pointer.
                MatrixPtr->use_sparse = parser.use_sparse;  // Set sparse matrix flag based on netlist directive
                CKTsetup(ckt, parser, MatrixPtr, modmap); // Pass the parser to the ckt and the initialise LHS and RHS matrices
                CKTload(ckt);
                ckt.cktmatrix->set_initmatrix(); // Set the initial LHS and RHS matrices
                ckt.sim_stats.MNA_Matrix_size = ckt.cktmatrix->LHS.rows();    // Store the size of MNA matrix to the simulation statistics

                // Check for multithreading from parser
                if (ckt.CKTmultithreaded){
                    std::cout << "Multithreading enabled for BSIM4 transistors with " << ckt.num_threads << " threads." << std::endl;
                }
            }

            if(parser.is_op){
                bool non_linear = CKTisNonLinear(ckt.CKTelements);
                OPResult op_result = OP_ops(ckt, modmap, non_linear);
                printOperatingPointWithNames(op_result.solution, ckt.map);
                save_txt_op("op_solution.txt", op_result, ckt.map);
                std::cout << "Operating point simulation completed." << std::endl;
            }
            else if(parser.is_transient){
                TransientSimulator trans_sim = Transsetup(parser, ckt);
                std::vector<Transient> vec_trans_result = Transient_ops(ckt, trans_sim, modmap);
                save_csv("tran_solution.csv", ckt, vec_trans_result, ckt.map);
                std::cout << "Transient simulation completed." << std::endl;

            }
            else if(parser.is_dc){
                dc::DCSimulator dcSim = dc::DCsetup(parser, ckt);
                std::vector<dc::DCResult> vec_dc_result = dc::DC_ops(ckt, dcSim, modmap);
                save_csv_dc("dc_solution.csv", ckt, vec_dc_result, ckt.map);
                std::cout << "DC simulation completed." << std::endl;
            }
            else if(parser.is_ac){
                if(!CKTcheckAC(ckt.CKTelements)){throw SimulationException("Error: No AC sources found!", "CKTcheckAC");}
                CKTloadAC(ckt);
                ckt.cktmatrix->set_init_cxmatrix(); // Set the initial complex LHS and RHS matrices for AC analysis
                ac::ACsimulator acSim = ac::ACsetup(parser, ckt);
                std::vector<ac::ACResult> vec_ac_result = ac::AC_ops(ckt, acSim, modmap);
                save_csv_ac("ac_solution.csv", ckt, vec_ac_result, ckt.map, acSim.type);
                std::cout << "AC simulation completed." << std::endl;
            }
            else{
                throw SimulationException("Error: No simulation type specified!", "main");
            }

            if(parser.acct){
                total.stop();
                ckt.sim_stats.simTime.total_time = total;
                ckt.sim_stats.simTime.parse_time = parser.parseTimer;
                ckt.sim_stats.printStatistics();
            }
        }

    }catch (const SimulationException& e) {
        std::cerr << "Simulation failed: " << e.what() << std::endl;
        return 1; // Exit gracefully instead of exit(1)
    } catch (const std::exception& e) {
        std::cerr << "Unexpected error: " << e.what() << std::endl;
        return 1;
    }


    return 0;
}

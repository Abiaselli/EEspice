#define ARMA_DONT_USE_WRAPPER
#define ARMA_USE_MKL_ALLOC
// #define ARMA_USE_SUPERLU

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
        auto batch_results = batch::run_batch_simulation(cktmap, parser, modmap);
        std::cout << "Batch simulation finished. Saving " << batch_results.size() << " results." << std::endl;
        batch::save_csv_batch(batch_results);
    }
    else{
        try {
            // CKT circuit setup
            CKTcircuit ckt;
            ckt.map = cktmap; // Assign the circuit map to the CKTcircuit
            auto denseMatrixPtr = std::make_shared<DenseMatrix>();  // Create the DenseMatrix as a shared pointer.
            CKTsetup(ckt, parser, denseMatrixPtr, modmap); // Pass the parser to the ckt and the initialise LHS and RHS matrices
            CKTload(ckt);
            ckt.cktdematrix->set_initmatrix(); // Set the initial LHS and RHS matrices

            if(parser.is_op){
                bool non_linear = false;
                if(!ckt.CKTelements.nmos.empty() || !ckt.CKTelements.pmos.empty() || !ckt.CKTelements.diodes.empty()){
                    non_linear = true;
                }
                OPResult op_result = OP_ops(ckt, modmap, non_linear);
                printOperatingPoint(op_result.solution, ckt);
            }
            if(parser.is_transient){
                TransientSimulator trans_sim = Transsetup(parser, ckt);
                std::vector<Transient> vec_trans_result = Transient_ops(ckt, trans_sim, modmap);
                save_csv("tran_solution.csv", ckt, vec_trans_result, ckt.map);

            }
            if(parser.is_dc){
                dc::DCSimulator dcSim = dc::DCsetup(parser, ckt);
                std::vector<dc::DCResult> vec_dc_result = dc::DC_ops(ckt, dcSim, modmap);
                save_csv_dc("dc_solution.csv", ckt, vec_dc_result, ckt.map);
            }
            if(parser.is_ac){
                CKTloadAC(ckt);
                ckt.cktdematrix->set_init_cxmatrix(); // Set the initial complex LHS and RHS matrices for AC analysis
                AC::ACsimulator acSim = AC::ACsetup(parser, ckt);
                std::vector<AC::AC> vec_ac_result = AC::AC_ops(ckt, acSim, modmap);
            }
        } catch (const SimulationException& e) {
            std::cerr << "Simulation failed: " << e.what() << std::endl;
            return 1; // Exit gracefully instead of exit(1)
        } catch (const std::exception& e) {
            std::cerr << "Unexpected error: " << e.what() << std::endl;
            return 1;
        }
    }

    /*-----------------------------------------------------------------------------------------------------------*/
    // SAVING THE SOLUTION AND TIME MATRICES INTO CSV FILES

    // std::chrono::duration<double, std::milli> time_span = (t3 - t1);
    // std::chrono::duration<double, std::milli> analysis_time = (tstop_trans - t1);
    // std::cout << "Total analysis time:" << (analysis_time).count()  << "ms\n";
    // std::cout << "Total time:" << time_span.count() << "ms\n";
    // std::cout << "Total NR iteration:" << total_NR_iteration << std::endl;
    // std::cout << "Total timepoint:" << total_timepoint << std::endl;

    return 0;
    /*-----------------------------------------------------------------------------------------------------------*/
}

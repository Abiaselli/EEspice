#define ARMA_DONT_USE_WRAPPER
#define ARMA_USE_MKL_ALLOC
// #define ARMA_USE_SUPERLU

#include "model_setup.hpp"
#include "Transient.hpp"
#include "DC.hpp"
#include "AC.hpp"
#include "saveCSV.hpp"
#include "batch.hpp"

// Main function for the circuit simulation
int main(int argc, const char **argv)
{

    setDebugMode(false); // Set the debug mode to false or true

    // for(int i = 0; i < 3; i++){
    //     pool.detach_task(dummy_task);
    // }
    // pool.wait();

    auto t1 = std::chrono::high_resolution_clock::now(); // Start time
    // Parse netlist file
    Modelmap modmap;
    Circuitmap cktmap;
    CircuitParser parser("Netlist/batch.cir");
    parser_netlist(parser, cktmap, modmap);

    // Model setup using the temperature
    modelSetup(modmap, nomTemp);

    if(parser.is_batch){
        std::cout << "Starting batch simulation..." << std::endl;
        auto batch_results = batch::run_batch_simulation(cktmap, parser, modmap);
        std::cout << "Batch simulation finished. Saving " << batch_results.size() << " results." << std::endl;
    }
    else{
        // CKT circuit setup
        CKTcircuit ckt;
        ckt.map = cktmap; // Assign the circuit map to the CKTcircuit
        auto denseMatrixPtr = std::make_shared<DenseMatrix>();  // Create the DenseMatrix as a shared pointer.
        CKTsetup(ckt, parser, denseMatrixPtr, modmap); // Pass the parser to the ckt and the initialise LHS and RHS matrices
        CKTload(ckt);
        ckt.cktdematrix->set_initmatrix(); // Set the initial LHS and RHS matrices

        if(parser.is_transient){
            TransientSimulator trans_sim = Transsetup(parser, ckt);
            std::vector<Transient> vec_trans_result = Transient_ops(ckt, trans_sim, modmap);
            save_csv(ckt, vec_trans_result, ckt.map);

        }
        if(parser.is_dc){
            DCSimulator dcSim = dc::DCsetup(parser, ckt);
            std::vector<DC> vec_dc_result = dc::DC_ops(ckt, dcSim, modmap);
            save_csv_dc(ckt, vec_dc_result, ckt.map);
        }
        if(parser.is_ac){
            CKTloadAC(ckt);
            ckt.cktdematrix->set_init_cxmatrix(); // Set the initial complex LHS and RHS matrices for AC analysis
            AC::ACsimulator acSim = AC::ACsetup(parser, ckt);
            std::vector<AC::AC> vec_ac_result = AC::AC_ops(ckt, acSim, modmap);
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

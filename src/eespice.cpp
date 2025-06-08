#define ARMA_DONT_USE_WRAPPER
#define ARMA_USE_MKL_ALLOC
// #define ARMA_USE_SUPERLU

#include "model_setup.hpp"
#include "Transient.hpp"
#include "DC.hpp"
#include "saveCSV.hpp"

// Main function for the circuit simulation
int main(int argc, const char **argv)
{

    setDebugMode(false); // Set the debug mode to false or true

    // for(int i = 0; i < 3; i++){
    //     pool.detach_task(dummy_task);
    // }
    // pool.wait();

    auto t1 = std::chrono::high_resolution_clock::now(); // Start time

    CKTcircuit ckt;
    auto denseMatrixPtr = std::make_shared<DenseMatrix>();  // Create the DenseMatrix as a shared pointer.

    Modelmap modmap;

    CircuitParser parser("Netlist/Ring.cir");
    parser_netlist(parser, ckt.map, modmap);

    // Model setup using the temperature
    modelSetup(modmap, nomTemp);

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

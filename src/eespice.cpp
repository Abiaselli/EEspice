#define ARMA_DONT_USE_WRAPPER
#define ARMA_USE_MKL_ALLOC
// #define ARMA_USE_SUPERLU

#include "Transient.hpp"
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
    Circuitmap map;

    CircuitParser parser("Netlist/Ring.cir");
    parser_netlist(parser, map);

    CKTsetup(ckt, parser, denseMatrixPtr); // Pass the parser to the ckt and the initialise LHS and RHS matrices

    CKTload(ckt);
    ckt.cktdematrix->set_initmatrix(); // Set the initial LHS and RHS matrices

    TransientSimulator trans_sim = Transsetup(parser, ckt);

    auto vec_trans_result = Transient_ops(ckt, trans_sim);

    auto tstop_trans = std::chrono::high_resolution_clock::now();

    /*-----------------------------------------------------------------------------------------------------------*/
    // SAVING THE SOLUTION AND TIME MATRICES INTO CSV FILES
    auto t2 = std::chrono::high_resolution_clock::now(); // End time

    save_csv(ckt, vec_trans_result, map);

    auto t3 = std::chrono::high_resolution_clock::now(); // End time

    std::chrono::duration<double, std::milli> time_span = (t3 - t1);
    std::chrono::duration<double, std::milli> analysis_time = (tstop_trans - t1);
    std::cout << "Total analysis time:" << (analysis_time).count()  << "ms\n";
    std::cout << "Total time:" << time_span.count() << "ms\n";
    std::cout << "Total NR iteration:" << total_NR_iteration << std::endl;
    std::cout << "Total timepoint:" << total_timepoint << std::endl;

    return 0;
    /*-----------------------------------------------------------------------------------------------------------*/
}

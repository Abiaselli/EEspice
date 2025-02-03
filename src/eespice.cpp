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
    DenseMatrix dematrix;
    Circuitmap map;

    CircuitParser parser("Netlist/Ring.cir");
    parser_netlist(parser, map);

    CKTsetup(ckt, parser, dematrix); // Pass the parser to the ckt and the initialise LHS and RHS matrices
    ckt.setcktmatrix(dematrix);

    CKTload(ckt);
    ckt.cktdematrix->set_initmatrix(); // Set the initial LHS and RHS matrices

    TransientSimulator trans_sim;
    Transsetup(trans_sim.trans_config, parser, ckt);

    auto vec_trans_result = Transient_ops(ckt, dematrix, trans_sim);

    auto tstop_trans = std::chrono::high_resolution_clock::now();

    /* Getting number of milliseconds as a double. */
    // std::chrono::duration<double, std::milli> OP_time = (tstop_op - tstart_op);
    // std::chrono::duration<double, std::milli> trans_time = (tstop_trans - tstart_trans);

    // std::cout << "DC OP time:" << OP_time.count() << "ms\n";
    // std::cout << "Transient time:" << trans_time.count() << "ms\n";

    /*-----------------------------------------------------------------------------------------------------------*/
    // SAVING THE SOLUTION AND TIME MATRICES INTO CSV FILES
    auto t2 = std::chrono::high_resolution_clock::now(); // End time

    save_csv(ckt, vec_trans_result, map);

    auto t3 = std::chrono::high_resolution_clock::now(); // End time

    std::chrono::duration<double, std::milli> time_span = (t3 - t1);
    std::chrono::duration<double, std::milli> analysis_time = (tstop_trans - t1);
    std::cout << "Total analysis time:" << (analysis_time).count()  << "ms\n";
    std::cout << "Total time:" << time_span.count() << "ms\n";
    // std::cout << "The Total time for multi-solver is: " << timer.total_ms() << "ms\n";

    // save_threads_time(t1, tstop_trans);

    return 0;
    /*-----------------------------------------------------------------------------------------------------------*/
}

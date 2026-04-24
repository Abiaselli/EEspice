#include <omp.h>
#include <filesystem>
#include "model_setup.hpp"
#include "Transient.hpp"
#include "DC.hpp"
#include "AC.hpp"
#include "saveCSV.hpp"
#include "saveRAW.hpp"
#include "batch.hpp"
#include "simulation_exceptions.hpp"

namespace {
enum class OutputKind { CSV, RAW };

// Extension rules: .csv / .txt -> CSV (text/legacy), everything else
// (including .raw, no extension, unknown) -> RAW (ngspice-compatible).
OutputKind resolve_output_kind(const std::string &path) {
    namespace fs = std::filesystem;
    std::string ext = fs::path(path).extension().string();
    for (auto &c : ext) c = static_cast<char>(std::tolower(c));
    if (ext == ".csv") return OutputKind::CSV;
    if (ext == ".txt") return OutputKind::CSV;  // text/CSV branch
    return OutputKind::RAW;
}
} // anonymous namespace

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
            // Batch default flips to a single multi-plot .raw file. If
            // .output points at something with a .csv extension (or a bare
            // directory name like "batch_results"), fall back to the
            // directory-of-CSVs layout. OutputKind on an empty extension
            // returns RAW, so the no-.output case picks up the new default.
            const std::string batch_default_raw = "batch_results.raw";
            const std::string batch_path = parser.output_path.empty()
                                               ? batch_default_raw
                                               : parser.output_path;
            if (resolve_output_kind(batch_path) == OutputKind::RAW) {
                save_raw_batch(batch_results, batch_path);
            } else {
                batch::save_csv_batch(batch_results, parser.output_path);
            }
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
                CKTdiscoverPattern(ckt); // Discover pattern, lock it, and save baselines for efficient NR iteration resets
                ckt.sim_stats.MNA_Matrix_size = ckt.cktmatrix->LHS.rows();    // Store the size of MNA matrix to the simulation statistics

                // Check for multithreading from parser
                if (ckt.CKTmultithreaded){
                    std::cout << "Multithreading enabled for BSIM4 transistors with " << ckt.num_threads << " threads." << std::endl;
                }
            }

            auto resolve_output = [&](const std::string &default_name) -> std::string {
                std::string path = parser.output_path.empty() ? default_name : parser.output_path;
                std::filesystem::path parent = std::filesystem::path(path).parent_path();
                if (!parent.empty()) {
                    std::filesystem::create_directories(parent);
                }
                return path;
            };

            if(parser.is_op){
                bool non_linear = CKTisNonLinear(ckt.CKTelements);
                OPResult op_result = OP_ops(ckt, modmap, non_linear);
                printOperatingPointWithNames(op_result.solution, ckt.map);
                std::string path = resolve_output("op_solution.raw");
                if (resolve_output_kind(path) == OutputKind::RAW) {
                    save_raw_op(path, ckt, op_result, ckt.map);
                } else {
                    // .csv / .txt both route to the existing human-readable text report.
                    save_txt_op(path, op_result, ckt.map);
                }
                std::cout << "Operating point simulation completed." << std::endl;
            }
            else if(parser.is_transient){
                TransientSimulator trans_sim = Transsetup(parser, ckt);
                std::vector<Transient> vec_trans_result = Transient_ops(ckt, trans_sim, modmap);
                std::string path = resolve_output("tran_solution.raw");
                if (resolve_output_kind(path) == OutputKind::RAW) {
                    save_raw_tran(path, ckt, vec_trans_result, ckt.map);
                } else {
                    save_csv(path, ckt, vec_trans_result, ckt.map);
                }
                std::cout << "Transient simulation completed." << std::endl;

            }
            else if(parser.is_dc){
                dc::DCSimulator dcSim = dc::DCsetup(parser, ckt);
                std::vector<dc::DCResult> vec_dc_result = dc::DC_ops(ckt, dcSim, modmap);
                std::string path = resolve_output("dc_solution.raw");
                if (resolve_output_kind(path) == OutputKind::RAW) {
                    save_raw_dc(path, ckt, vec_dc_result, ckt.map);
                } else {
                    save_csv_dc(path, ckt, vec_dc_result, ckt.map);
                }
                std::cout << "DC simulation completed." << std::endl;
            }
            else if(parser.is_ac){
                if(!CKTcheckAC(ckt.CKTelements)){throw SimulationException("Error: No AC sources found!", "CKTcheckAC");}
                CKTloadAC(ckt);
                ckt.cktmatrix->set_init_cxmatrix(); // Set the initial complex LHS and RHS matrices for AC analysis
                ac::ACsimulator acSim = ac::ACsetup(parser, ckt);
                std::vector<ac::ACResult> vec_ac_result = ac::AC_ops(ckt, acSim, modmap);
                std::string path = resolve_output("ac_solution.raw");
                if (resolve_output_kind(path) == OutputKind::RAW) {
                    save_raw_ac(path, ckt, vec_ac_result, ckt.map);
                } else {
                    save_csv_ac(path, ckt, vec_ac_result, ckt.map, acSim.type);
                }
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

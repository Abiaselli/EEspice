#pragma once

// Binary ngspice-compatible .raw output writer for EESpice.
// RAW schema emitted by this writer:
//
//   A file contains one or more plot blocks. Batch mode writes a multi-plot
//   file by appending complete plot blocks back-to-back.
//
//   Per plot:
//     Title: <text>
//     Date: <UTC timestamp>
//     Command: ngspice-compatible EESpice, <build-or-batch-metadata>
//     Plotname: op<N> | tran<N> | dc<N> | ac<N>
//     Flags: real | complex
//     No. Variables: <nvars>
//     No. Points: <npoints>
//     Variables:
//       <index>\t<name>\t<type>
//       ...
//     Binary:
//       payload
//
//   Payload is point-major: for each point p, variables are written in the
//   exact Variables order. Real plots store one native binary double per
//   variable. Complex plots store two native binary doubles per variable:
//   (real, imag). Real-valued scale variables in complex plots, such as AC
//   frequency, are encoded as (value, 0.0). EESpice's supported Linux/x86_64
//   output is little-endian, although the C++ writer emits native-endian
//   doubles.
//
//   Variable order by analysis:
//     OP:   v(<node>)..., i(<branch>)... with no explicit scale variable.
//     TRAN: time, v(<node>)..., i(<branch>)...
//     DC:   <sweep-source>, v(<node>)..., i(<branch>)...
//     AC:   frequency, v(<node>)..., i(<branch>)... as complex values.
//
// See helper/decode_raw.py for a small Python stdlib decoder example.

#include <complex>
#include <ctime>
#include <cstdint>
#include <fstream>
#include <functional>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

#include "AC_calcs.hpp"
#include "CKT.hpp"
#include "DC.hpp"
#include "OP_calcs.hpp"
#include "Transient.hpp"
#include "saveCSV.hpp"
#include "simulation_exceptions.hpp"

namespace {

struct RawVarDesc {
    std::string name;
    std::string type;
};

// WHY this depends on saveCSV.hpp: we reuse `buildNodeIndexToNameMap` to
// enumerate external nodes in index order. There are exactly two writer
// backends (CSV, RAW) and no third is planned, so a shared abstraction
// would be speculative.
//
// WHY helpers live in an anonymous namespace in a header: in this project
// only eespice.cpp includes saveRAW.hpp, so the per-TU copy footprint is
// one copy. If a second TU ever includes this header, split into
// saveRAW.hpp / saveRAW.cpp.

// Writes the header block up through the "No. Points" line, matching
// ngspice rawfile.c:115-122. Flags line has NO trailing whitespace.
// Date format pinned to "%a %b %e %H:%M:%S %Y" (UTC) for cross-machine
// determinism (see plan Step 2 / Risk #10).
void write_raw_header_block(std::ofstream &file,
                            const std::string &title,
                            const std::string &plotname,
                            bool is_complex,
                            int nvars,
                            int npoints,
                            const std::string &command_suffix = "")
{
    std::time_t t = std::time(nullptr);
    std::tm *gm = std::gmtime(&t);
    std::ostringstream date_ss;
    date_ss << std::put_time(gm, "%a %b %e %H:%M:%S %Y");

    file << "Title: " << title << "\n";
    file << "Date: " << date_ss.str() << "\n";
    // Command line: identifies the producing simulator. spicelib detects the
    // raw-file dialect by matching "ngspice" / "ltspice" substrings here — so
    // we advertise "ngspice-compatible" explicitly, otherwise spicelib throws
    // SpiceReadException with "file dialect is not specified". The default
    // suffix is "Build __DATE__ __TIME__" pinning the build identity; batch
    // runs override it with per-plot config metadata (e.g. "Batch: m1.W=...")
    // so multi-plot RAW readers can map each plot back to its parameter values.
    file << "Command: ngspice-compatible EESpice, ";
    if (command_suffix.empty()) {
        file << "Build " << __DATE__ << " " << __TIME__;
    } else {
        file << command_suffix;
    }
    file << "\n";
    file << "Plotname: " << plotname << "\n";
    file << "Flags: " << (is_complex ? "complex" : "real") << "\n";
    file << "No. Variables: " << nvars << "\n";
    file << "No. Points: " << npoints << "\n";
}

// Writes the Variables section, tab-separated `index\tname\ttype` per
// rawfile.c:163-221, then the "Binary:\n" delimiter that precedes the
// binary payload.
void write_raw_variables_block(std::ofstream &file,
                               const std::vector<RawVarDesc> &vars)
{
    file << "Variables:\n";
    for (size_t i = 0; i < vars.size(); ++i) {
        file << "\t" << i << "\t" << vars[i].name << "\t" << vars[i].type << "\n";
    }
    file << "Binary:\n";
}

// Row-major binary payload: outer loop over points, inner over variables
// (rawfile.c:223-258). Each value is a native-endian little-endian double.
void write_raw_binary_real(std::ofstream &file,
                           int n_points, int n_vars,
                           std::function<double(int point, int var)> getter)
{
    for (int p = 0; p < n_points; ++p) {
        for (int v = 0; v < n_vars; ++v) {
            double value = getter(p, v);
            file.write(reinterpret_cast<const char *>(&value), sizeof(double));
        }
    }
}

// For complex plots, every variable is written as (real, imag) pair —
// real-valued scales (frequency) are coerced to (v, 0.0) per rawfile.c:234-238.
// `is_complex_var[v]` gates whether to take imag() from the callback's value
// or zero it out.
void write_raw_binary_complex(std::ofstream &file,
                              int n_points, int n_vars,
                              std::function<std::complex<double>(int point, int var)> getter,
                              const std::vector<bool> &is_complex_var)
{
    for (int p = 0; p < n_points; ++p) {
        for (int v = 0; v < n_vars; ++v) {
            std::complex<double> z = getter(p, v);
            double re = z.real();
            double im = is_complex_var[v] ? z.imag() : 0.0;
            file.write(reinterpret_cast<const char *>(&re), sizeof(double));
            file.write(reinterpret_cast<const char *>(&im), sizeof(double));
        }
    }
}

// Formats a BatchRunResult's CircuitConfig as the trailing part of the RAW
// Command: line, e.g. "Batch: m1.W=1.000000e-04, m1.L=9.000000e-07".
// Returns "" for an empty config so write_raw_header_block falls back to the
// default Build timestamp. Precision matches saveCSV.hpp:305 so CSV and RAW
// agree on parameter values. CircuitConfig is std::map (batch.hpp:22) so the
// iteration order is alphabetical by key — deterministic across runs.
std::string make_batch_command_suffix(const batch::CircuitConfig &config)
{
    if (config.empty()) return "";
    std::ostringstream ss;
    ss << "Batch: ";
    bool first = true;
    for (const auto &[name, value] : config) {
        if (!first) ss << ", ";
        ss << name << "="
           << std::scientific << std::setprecision(6) << value;
        first = false;
    }
    return ss.str();
}

// Opens `path` for binary output, surfacing filesystem failure as the
// project-standard SimulationException (mirrors saveCSV.hpp:67-71 pattern).
// 5 callers: each `save_raw_*(filename, ...)` overload plus save_raw_batch.
std::ofstream open_raw_file(const std::string &path)
{
    std::ofstream file(path, std::ios::binary);
    if (!file.is_open()) {
        throw SimulationException(
            "Error: Could not open file '" + path + "' for writing.",
            "OUTPUT_FILE_OPEN_FAILED");
    }
    return file;
}

} // anonymous namespace


// -----------------------------------------------------------------------------
// OP writer
//
// Plotname: op1, Flags: real, No. Points: 1. No separate scale variable:
// variable 0 is the first node voltage (ngspice rawfile.c:150-160 treats
// variable 0 as scale when pl_scale is absent; spicelib tolerates this).
// -----------------------------------------------------------------------------
void save_raw_op(std::ofstream &file, const CKTcircuit &ckt, const OPResult &op,
                 const Circuitmap &map, int plot_index = 1,
                 const std::string &command_suffix = "")
{
    std::vector<std::string> nodeIndexToName = buildNodeIndexToNameMap(ckt, map);

    std::vector<RawVarDesc> vars;
    vars.reserve(ckt.external_nodes + map.map_branch_currents.size());
    for (int j = 1; j <= ckt.external_nodes; ++j) {
        vars.push_back({"v(" + nodeIndexToName[j] + ")", "voltage"});
    }
    // map_branch_currents is unordered; iteration order is implementation-defined,
    // but consistent within a single run — spicelib does not rely on any order.
    std::vector<int> branch_indices;
    branch_indices.reserve(map.map_branch_currents.size());
    for (const auto &[name, idx] : map.map_branch_currents) {
        vars.push_back({"i(" + name + ")", "current"});
        branch_indices.push_back(idx);
    }

    const int nvars = static_cast<int>(vars.size());
    const int npoints = 1;
    const std::string plotname = "op" + std::to_string(plot_index);

    write_raw_header_block(file, "EESpice output", plotname,
                           /*is_complex=*/false, nvars, npoints, command_suffix);
    write_raw_variables_block(file, vars);

    write_raw_binary_real(file, npoints, nvars,
        [&](int /*point*/, int v) -> double {
            // External nodes occupy solution(0 .. external_nodes-1); solution
            // is zero-indexed (arma::vec), node ids are 1-indexed, so v == 0
            // maps to solution(0), i.e., node id 1.
            if (v < ckt.external_nodes) {
                return op.solution(v);
            }
            return op.solution(branch_indices[v - ckt.external_nodes]);
        });
}

void save_raw_op(const std::string &filename, const CKTcircuit &ckt,
                 const OPResult &op, const Circuitmap &map)
{
    std::ofstream file = open_raw_file(filename);
    save_raw_op(file, ckt, op, map, /*plot_index=*/1);
    file.close();
    std::cout << "Data saved to " << filename << std::endl;
}

// -----------------------------------------------------------------------------
// TRAN writer
//
// Plotname: tran<N>, Flags: real. Scale variable "time" (type "time") at
// index 0. Then v(<node>) per external node, then i(<dev>) per branch
// current. Row-major binary: (time, v(1), v(2), ..., i(dev0), i(dev1), ...).
// Note: we deliberately do NOT emit the CSV-only "Time Step" column — that
// is not an ngspice concept and spicelib does not expect it.
// -----------------------------------------------------------------------------
void save_raw_tran(std::ofstream &file, const CKTcircuit &ckt,
                   const std::vector<Transient> &vec_trans,
                   const Circuitmap &map, int plot_index = 1,
                   const std::string &command_suffix = "")
{
    std::vector<std::string> nodeIndexToName = buildNodeIndexToNameMap(ckt, map);

    std::vector<RawVarDesc> vars;
    vars.reserve(1 + ckt.external_nodes + map.map_branch_currents.size());
    vars.push_back({"time", "time"});
    for (int j = 1; j <= ckt.external_nodes; ++j) {
        vars.push_back({"v(" + nodeIndexToName[j] + ")", "voltage"});
    }
    // Freeze the branch-current iteration order into a parallel vector so
    // the Variables block and the binary payload stay in lock-step
    // (map.map_branch_currents is unordered_map — iteration order is stable
    // within one run but not across runs; that is fine since the Variables
    // block is written in the same iteration).
    std::vector<int> branch_indices;
    branch_indices.reserve(map.map_branch_currents.size());
    for (const auto &[name, idx] : map.map_branch_currents) {
        vars.push_back({"i(" + name + ")", "current"});
        branch_indices.push_back(idx);
    }

    const int nvars = static_cast<int>(vars.size());
    const int npoints = static_cast<int>(vec_trans.size());
    const std::string plotname = "tran" + std::to_string(plot_index);

    write_raw_header_block(file, "EESpice output", plotname,
                           /*is_complex=*/false, nvars, npoints, command_suffix);
    write_raw_variables_block(file, vars);

    write_raw_binary_real(file, npoints, nvars,
        [&](int p, int v) -> double {
            if (v == 0) return vec_trans[p].time_trans;
            const int node_or_branch = v - 1;
            if (node_or_branch < ckt.external_nodes) {
                // External nodes: solution slots 0..external_nodes-1 align
                // with node ids 1..external_nodes (zero-indexed arma::vec,
                // one-indexed node map — see buildNodeIndexToNameMap).
                return vec_trans[p].solution(node_or_branch);
            }
            return vec_trans[p].solution(
                branch_indices[node_or_branch - ckt.external_nodes]);
        });
}

void save_raw_tran(const std::string &filename, const CKTcircuit &ckt,
                   const std::vector<Transient> &vec_trans,
                   const Circuitmap &map)
{
    std::ofstream file = open_raw_file(filename);
    save_raw_tran(file, ckt, vec_trans, map, /*plot_index=*/1);
    file.close();
    std::cout << "Data saved to " << filename << std::endl;
}

// -----------------------------------------------------------------------------
// DC writer
//
// Plotname: dc<N>, Flags: real. Scale variable = the sweep source id
// (DCResult::sweepName), type fixed to "voltage" — DCResult has no
// source-type tag (src/DC_calcs.hpp:24-27) and both ngspice and spicelib
// readers ignore the type string on the scale variable. Non-scale
// variables follow: v(<node>) per external node, then i(<dev>).
// -----------------------------------------------------------------------------
void save_raw_dc(std::ofstream &file, const CKTcircuit &ckt,
                 const std::vector<dc::DCResult> &vec_dc,
                 const Circuitmap &map, int plot_index = 1,
                 const std::string &command_suffix = "")
{
    if (vec_dc.empty()) {
        // Nothing to write — mirror the CSV writer's behavior at
        // saveCSV.hpp:83-86 and skip the body. File still ends up with
        // just the header which ngspice/spicelib will treat as empty.
        std::cerr << "Warning: vec_dc is empty, no DC solutions to write.\n";
        return;
    }

    std::vector<std::string> nodeIndexToName = buildNodeIndexToNameMap(ckt, map);

    std::vector<RawVarDesc> vars;
    vars.reserve(1 + ckt.external_nodes + map.map_branch_currents.size());
    vars.push_back({vec_dc.front().sweepName, "voltage"});
    for (int j = 1; j <= ckt.external_nodes; ++j) {
        vars.push_back({"v(" + nodeIndexToName[j] + ")", "voltage"});
    }
    std::vector<int> branch_indices;
    branch_indices.reserve(map.map_branch_currents.size());
    for (const auto &[name, idx] : map.map_branch_currents) {
        vars.push_back({"i(" + name + ")", "current"});
        branch_indices.push_back(idx);
    }

    const int nvars = static_cast<int>(vars.size());
    const int npoints = static_cast<int>(vec_dc.size());
    const std::string plotname = "dc" + std::to_string(plot_index);

    write_raw_header_block(file, "EESpice output", plotname,
                           /*is_complex=*/false, nvars, npoints, command_suffix);
    write_raw_variables_block(file, vars);

    write_raw_binary_real(file, npoints, nvars,
        [&](int p, int v) -> double {
            if (v == 0) return vec_dc[p].sweepValue;
            const int node_or_branch = v - 1;
            if (node_or_branch < ckt.external_nodes) {
                return vec_dc[p].solution(node_or_branch);
            }
            return vec_dc[p].solution(
                branch_indices[node_or_branch - ckt.external_nodes]);
        });
}

void save_raw_dc(const std::string &filename, const CKTcircuit &ckt,
                 const std::vector<dc::DCResult> &vec_dc,
                 const Circuitmap &map)
{
    std::ofstream file = open_raw_file(filename);
    save_raw_dc(file, ckt, vec_dc, map, /*plot_index=*/1);
    file.close();
    std::cout << "Data saved to " << filename << std::endl;
}

// -----------------------------------------------------------------------------
// AC writer
//
// Plotname: ac<N>, Flags: complex. Scale variable "frequency" (type
// "frequency") at index 0, coerced to (freq, 0.0) because the file is
// complex (rawfile.c:234-238: when ANY variable is complex, ALL variables
// are written as 2 doubles per point; reals get 0.0 imag).
//
// Non-scale variables are v(<node>) / i(<dev>), each stored as the native
// (real, imag) pair from the arma::cx_dvec solution. We deliberately emit
// real+imag here rather than mag/phase — CSV covers mag/phase already.
// -----------------------------------------------------------------------------
void save_raw_ac(std::ofstream &file, const CKTcircuit &ckt,
                 const std::vector<ac::ACResult> &vec_ac,
                 const Circuitmap &map, int plot_index = 1,
                 const std::string &command_suffix = "")
{
    std::vector<std::string> nodeIndexToName = buildNodeIndexToNameMap(ckt, map);

    std::vector<RawVarDesc> vars;
    std::vector<bool>       is_complex_var;
    const int n_branches = static_cast<int>(map.map_branch_currents.size());
    vars.reserve(1 + ckt.external_nodes + n_branches);
    is_complex_var.reserve(1 + ckt.external_nodes + n_branches);

    vars.push_back({"frequency", "frequency"});
    is_complex_var.push_back(false);  // frequency coerced to (f, 0.0)
    for (int j = 1; j <= ckt.external_nodes; ++j) {
        vars.push_back({"v(" + nodeIndexToName[j] + ")", "voltage"});
        is_complex_var.push_back(true);
    }
    std::vector<int> branch_indices;
    branch_indices.reserve(n_branches);
    for (const auto &[name, idx] : map.map_branch_currents) {
        vars.push_back({"i(" + name + ")", "current"});
        is_complex_var.push_back(true);
        branch_indices.push_back(idx);
    }

    const int nvars = static_cast<int>(vars.size());
    const int npoints = static_cast<int>(vec_ac.size());
    const std::string plotname = "ac" + std::to_string(plot_index);

    write_raw_header_block(file, "EESpice output", plotname,
                           /*is_complex=*/true, nvars, npoints, command_suffix);
    write_raw_variables_block(file, vars);

    write_raw_binary_complex(file, npoints, nvars,
        [&](int p, int v) -> std::complex<double> {
            if (v == 0) {
                // frequency scale — real-valued; real part carries value,
                // imag forced to 0.0 by write_raw_binary_complex via
                // is_complex_var[0]==false.
                return std::complex<double>(vec_ac[p].freq, 0.0);
            }
            const int node_or_branch = v - 1;
            std::complex<double> z;
            if (node_or_branch < ckt.external_nodes) {
                z = vec_ac[p].solution(node_or_branch);
            } else {
                z = vec_ac[p].solution(
                    branch_indices[node_or_branch - ckt.external_nodes]);
            }
            return z;
        },
        is_complex_var);
}

void save_raw_ac(const std::string &filename, const CKTcircuit &ckt,
                 const std::vector<ac::ACResult> &vec_ac,
                 const Circuitmap &map)
{
    std::ofstream file = open_raw_file(filename);
    save_raw_ac(file, ckt, vec_ac, map, /*plot_index=*/1);
    file.close();
    std::cout << "Data saved to " << filename << std::endl;
}

// -----------------------------------------------------------------------------
// Batch writer — multi-plot single file
//
// Opens ONE binary file and writes each successful run as a sequential plot
// (Plotname: dc1, dc2, tran1, ...), matching ngspice's native multi-plot
// layout (rawfile.c:76 uses "ab" append mode for the same purpose).
// Failed runs are reported to stderr and skipped — the binary format has no
// comment mechanism after the header, so representing a failure as a plot
// would lie to the reader. Per-type counters mirror save_csv_batch's scheme
// (saveCSV.hpp:322-325). Batch ordering is deterministic because batch.hpp
// collects results on the main thread after futures.wait().
// -----------------------------------------------------------------------------
void save_raw_batch(const std::vector<batch::BatchRunResult> &batch_results,
                    const std::string &output_path_override = "")
{
    const std::string path = output_path_override.empty()
                                 ? std::string("batch_results.raw")
                                 : output_path_override;

    // Ensure parent dir exists (matches eespice.cpp's resolve_output behavior).
    std::filesystem::path parent = std::filesystem::path(path).parent_path();
    if (!parent.empty()) {
        std::filesystem::create_directories(parent);
    }

    std::ofstream file = open_raw_file(path);

    int dc_counter = 1, tran_counter = 1, ac_counter = 1, op_counter = 1;
    int successful_count = 0, failed_count = 0;

    for (const auto &run_result : batch_results) {
        if (!run_result.success) {
            std::cerr << "Batch run failed (" << run_result.simulation_type
                      << "): " << run_result.error_type << " — "
                      << run_result.error_message << "; skipping in .raw output."
                      << std::endl;
            ++failed_count;
            continue;
        }
        const Circuitmap &map = run_result.ckt.map;
        const std::string suffix = make_batch_command_suffix(run_result.config);
        if (run_result.simulation_type == "op") {
            save_raw_op(file, run_result.ckt,
                        std::get<OPResult>(run_result.results),
                        map, op_counter++, suffix);
        } else if (run_result.simulation_type == "dc") {
            save_raw_dc(file, run_result.ckt,
                        std::get<std::vector<dc::DCResult>>(run_result.results),
                        map, dc_counter++, suffix);
        } else if (run_result.simulation_type == "tran") {
            save_raw_tran(file, run_result.ckt,
                          std::get<std::vector<Transient>>(run_result.results),
                          map, tran_counter++, suffix);
        } else if (run_result.simulation_type == "ac") {
            save_raw_ac(file, run_result.ckt,
                        std::get<std::vector<ac::ACResult>>(run_result.results),
                        map, ac_counter++, suffix);
        }
        ++successful_count;
    }
    file.close();

    std::cout << std::endl << "Batch Simulation Summary:" << std::endl;
    std::cout << "Total simulations: " << batch_results.size() << std::endl;
    std::cout << "Successful: " << successful_count << std::endl;
    std::cout << "Failed: " << failed_count << std::endl;
    std::cout << "Data saved to " << path << std::endl;
}

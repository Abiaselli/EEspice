#include <iostream>
#include <armadillo>
#include <omp.h>
#include <vector>
#include <string>
#include <chrono>
#include <iomanip>

struct Diode
{   
    std::string id_str;
    int id{};
    std::string nodePos_str, nodeNeg_str;
    int nodePos{}, nodeNeg{};
    double Is{};
    double VT{};
};

void DiodeCal(int node_x, int node_y, double Is, double VT, const arma::vec &solution)
{

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
        return;
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
}

void runSingleThreaded(const std::vector<Diode>& diodes, const arma::vec& solution, int iterations) {
    for (int iter = 0; iter < iterations; iter++) {
        for (size_t i = 0; i < diodes.size(); i++) {
            DiodeCal(diodes[i].nodePos, diodes[i].nodeNeg, diodes[i].Is, diodes[i].VT, solution);
        }
    }
}

void runParallel(const std::vector<Diode>& diodes, const arma::vec& solution, int iterations) {
    #pragma omp parallel for collapse(2)
    for (int iter = 0; iter < iterations; iter++) {
        for (size_t i = 0; i < diodes.size(); i++) {
            DiodeCal(diodes[i].nodePos, diodes[i].nodeNeg, diodes[i].Is, diodes[i].VT, solution);
        }
    }
}

int main(){
    // Example solution vector
    arma::vec solution = arma::vec(10, arma::fill::randu);

    // Vector of 8 diodes
    std::vector<Diode> diodes = {
        {"D1", 1, "n1", "n2", 1, 2, 1e-14, 0.026},
        {"D2", 2, "n2", "n3", 2, 3, 1e-14, 0.026},
        {"D3", 3, "n3", "n4", 3, 4, 1e-14, 0.026},
        {"D4", 4, "n4", "n5", 4, 5, 1e-14, 0.026},
        {"D5", 5, "n5", "n6", 5, 6, 1e-14, 0.026},
        {"D6", 6, "n6", "n7", 6, 7, 1e-14, 0.026},
        {"D7", 7, "n7", "n8", 7, 8, 1e-14, 0.026},
        {"D8", 8, "n8", "n9", 8, 9, 1e-14, 0.026}
    };

    const int iterations = 100000000;  // Run many iterations for accurate timing (>1s per test)

    std::cout << "OpenMP Benchmark: Diode Calculations\n";
    std::cout << "=====================================\n";
    std::cout << "Number of diodes: " << diodes.size() << "\n";
    std::cout << "Iterations: " << iterations << "\n";
    std::cout << "Total calculations: " << diodes.size() * iterations << "\n";
    std::cout << "CPU cores: " << omp_get_num_procs() << "\n";
    std::cout << "OpenMP max threads: " << omp_get_max_threads() << "\n\n";

    // Warm-up run to initialize caches
    runSingleThreaded(diodes, solution, 100);

    // Benchmark single-threaded execution
    auto start = std::chrono::high_resolution_clock::now();
    runSingleThreaded(diodes, solution, iterations);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> single_time = end - start;

    std::cout << "Single-threaded time: " << std::fixed << std::setprecision(6)
              << single_time.count() << " seconds\n";

    // Test different thread counts
    std::vector<int> thread_counts = {2, 4, 8, 16, 24};
    double best_speedup = 0;
    int best_threads = 1;

    std::cout << "\nTesting different thread counts:\n";
    std::cout << "================================\n";

    for (int num_threads : thread_counts) {
        if (num_threads > omp_get_num_procs()) {
            std::cout << "Skipping " << num_threads << " threads (exceeds available cores)\n";
            continue;
        }

        omp_set_num_threads(num_threads);

        start = std::chrono::high_resolution_clock::now();
        runParallel(diodes, solution, iterations);
        end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> parallel_time = end - start;

        double speedup = single_time.count() / parallel_time.count();
        double efficiency = (speedup / num_threads) * 100.0;

        std::cout << "Parallel (" << std::setw(2) << num_threads << " threads): "
                  << std::fixed << std::setprecision(6) << parallel_time.count()
                  << " seconds  |  Speedup: " << std::setprecision(2) << speedup
                  << "x  |  Efficiency: " << std::setprecision(1) << efficiency << "%\n";
        std::cout.flush();

        if (speedup > best_speedup) {
            best_speedup = speedup;
            best_threads = num_threads;
        }
    }

    std::cout << "\nBest configuration: " << best_threads << " threads with "
              << std::fixed << std::setprecision(2) << best_speedup << "x speedup\n";

    return 0;
}
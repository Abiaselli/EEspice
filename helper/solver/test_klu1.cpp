#define ARMA_DONT_USE_WRAPPER
#define ARMA_USE_MKL_ALLOC
#define ARMA_USE_SUPERLU
#include <iostream>
#include <armadillo>
#include <klu.h>
#include <chrono>

// Helper concept to ensure we are actually on a 64-bit system
static_assert(sizeof(arma::uword) == 8, "This optimization requires 64-bit Armadillo indices");
static_assert(sizeof(SuiteSparse_long) == 8, "This optimization requires 64-bit SuiteSparse indices");

int KLUsolver(arma::sp_mat& A, arma::vec& b){
    // --- ZERO-COPY KLU INTEGRATION ---
    // 1. Cast Pointers
    auto* Ap = reinterpret_cast<SuiteSparse_long*>(const_cast<arma::uword*>(A.col_ptrs));
    auto* Ai = reinterpret_cast<SuiteSparse_long*>(const_cast<arma::uword*>(A.row_indices));
    double* Ax = const_cast<double*>(A.values); // Values are double*, no cast needed

    // 2. Initialize KLU (Use klu_l_defaults for long)
    klu_l_common Common;
    klu_l_defaults(&Common);

    // 3. Analyze (Use klu_l_analyze)
    klu_l_symbolic* Symbolic = klu_l_analyze(A.n_rows, Ap, Ai, &Common);
    if (!Symbolic) {
        std::cerr << "KLU Analyze failed." << std::endl;
        return 1;
    }

    // 4. Factor (Use klu_l_factor)
    klu_l_numeric* Numeric = klu_l_factor(Ap, Ai, Ax, Symbolic, &Common);
    if (!Numeric) {
        std::cerr << "KLU Factor failed." << std::endl;
        klu_l_free_symbolic(&Symbolic, &Common);
        return 1;
    }

    // 5. Solve (Use klu_l_solve)
    klu_l_solve(Symbolic, Numeric, A.n_rows, 1, b.memptr(), &Common);

    // 6. Cleanup (Use klu_l_free...)
    klu_l_free_symbolic(&Symbolic, &Common);
    klu_l_free_numeric(&Numeric, &Common);

    return 0;
}

int main() {
    constexpr int iterations = 2021;
    arma::sp_mat A;
    A.load("lhs_arma.txt", arma::coord_ascii); // Load a sparse matrix from a file
    std::cout << "A nnz: " << A.n_nonzero << ", size: " << A.n_rows << "x" << A.n_cols << std::endl;

    arma::vec b;
    b.load("rhs.txt", arma::raw_ascii); // Load the right-hand side vector

    double solve_time = 0.0;

    // warm-up run
    for (int i = 0; i < 10; ++i) {
        arma::sp_mat A_copy = A; // Create a copy to ensure zero-copy integration
        arma::vec b_copy = b;     // Copy the RHS vector
        KLUsolver(A_copy, b_copy);
    }

    for(int i = 0; i < iterations; ++i) {
        arma::sp_mat A_copy = A; // Create a copy to ensure zero-copy integration
        arma::vec b_copy = b;     // Copy the RHS vector
        auto start = std::chrono::high_resolution_clock::now();
        int status = KLUsolver(A_copy, b_copy);
        auto end = std::chrono::high_resolution_clock::now();
        solve_time += std::chrono::duration<double>(end - start).count();
    }
    std::cout << " KLU solve time over " << iterations << " iterations: " << solve_time << " seconds." << std::endl;
    return 0;
}
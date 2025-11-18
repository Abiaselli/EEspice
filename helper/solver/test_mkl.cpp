#define ARMA_DONT_USE_WRAPPER
#define ARMA_USE_MKL_ALLOC
#define ARMA_USE_SUPERLU
#include <armadillo>
#include <mkl_lapacke.h>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <vector>
#include <cmath>

// Generate a sparse matrix from dense matrix with specified sparsity
arma::sp_mat generateSparseMatrix(const arma::mat &dense, double sparsity_ratio) {
    int size = dense.n_rows;
    arma::sp_mat sparse(size, size);

    // Keep diagonal elements (always important for circuit matrices)
    for (int i = 0; i < size; i++) {
        sparse(i, i) = dense(i, i);
    }

    // Add random off-diagonal elements to achieve desired sparsity
    int target_nonzeros = static_cast<int>(size * size * sparsity_ratio);
    int current_nonzeros = size; // diagonal already added

    arma::arma_rng::set_seed(12345);
    while (current_nonzeros < target_nonzeros) {
        int i = arma::randi<int>(arma::distr_param(0, size - 1));
        int j = arma::randi<int>(arma::distr_param(0, size - 1));

        if (i != j && sparse(i, j) == 0.0) {
            sparse(i, j) = dense(i, j);
            sparse(j, i) = dense(j, i); // Keep symmetric
            current_nonzeros += 2;
        }
    }

    return sparse;
}

void solveDN(arma::mat &LHS, arma::vec &RHS, arma::vec &solution) {
    solution = arma::solve(LHS, RHS, arma::solve_opts::fast);
}
void solveSP(arma::sp_mat &LHS, arma::vec &RHS, arma::vec &solution) {
    solution = arma::spsolve(LHS, RHS, "superlu");
}
void solveMKL(arma::mat &LHS, arma::vec &RHS, std::vector<lapack_int> &ipiv){
    // Get matrix dimensions
    lapack_int n = static_cast<lapack_int>(LHS.n_rows);
    lapack_int nrhs = 1; // Single RHS vector

    // Get raw pointers to Armadillo data
    double* A_ptr = LHS.memptr();
    double* b_ptr = RHS.memptr();

    // Call MKL LAPACK solver
    // LAPACK_COL_MAJOR = 102 (matches Armadillo's column-major storage)
    lapack_int info = LAPACKE_dgesv(
        LAPACK_COL_MAJOR,  // Matrix layout
        n,                  // Order of matrix
        nrhs,              // Number of RHS vectors
        A_ptr,             // Matrix A (will be overwritten)
        n,                 // Leading dimension of A
        ipiv.data(),       // Pivot indices
        b_ptr,             // RHS vector (will be overwritten with solution)
        n                  // Leading dimension of b
    );
    // Check for errors
    if (info < 0) {
        throw std::runtime_error("LAPACKE_dgesv: Invalid parameter at position " +
                                 std::to_string(-info));
    } else if (info > 0) {
        throw std::runtime_error("LAPACKE_dgesv: Matrix is singular, U(" +
                                 std::to_string(info) + "," + std::to_string(info) +
                                 ") is exactly zero");
    }
}

int main() {
    constexpr int size = 1004;
    constexpr int iter = 2033;
    constexpr double sparsity_ratio = 0.004; // 0.4% nonzeros, typical for circuit matrices
    arma::mat LHS = arma::randn<arma::mat>(size, size);
    arma::vec RHS = arma::randn<arma::vec>(size);
    arma::vec solution_dense(size);

    // warm-up
    for (int i = 0; i < 5; i++) {
        solveDN(LHS, RHS, solution_dense);
    }

    auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < iter; i++) {
        solveDN(LHS, RHS, solution_dense);
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    double dense_time = elapsed.count();
    double dense_time_per_solve = dense_time / iter;

    std::cout << "Dense solver time for size " << size << "x" << size << " over " << iter << " iterations: "
              << dense_time << " seconds (" << dense_time_per_solve * 1e6 << " μs per solve)" << std::endl;

    // MKL solver
    // warm-up
    for (int i = 0; i < 5; i++){
        arma::mat LHS_copy = LHS;
        arma::vec RHS_copy = RHS;
        std::vector<lapack_int> ipiv(size);
        solveMKL(LHS_copy, RHS_copy, ipiv);
    }

    double mkl_time = 0.0;
    double mkl_time_per_solve = 0.0;

    for (int i = 0; i < iter; i++) {
        arma::mat LHS_copy = LHS;
        arma::vec RHS_copy = RHS;
        std::vector<lapack_int> ipiv(size);

        auto start_mkl = std::chrono::high_resolution_clock::now();
        solveMKL(LHS_copy, RHS_copy, ipiv);
        auto end_mkl = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed_mkl = end_mkl - start_mkl;
        mkl_time += elapsed_mkl.count();
    }
    mkl_time_per_solve = mkl_time / iter;

    std::cout << "MKL solver time for size " << size << "x" << size << " over " << iter << " iterations: "
              << mkl_time << " seconds (" << mkl_time_per_solve * 1e6 << " μs per solve)" << std::endl;

    // Sparse solver
    arma::sp_mat LHS_sparse = generateSparseMatrix(LHS, sparsity_ratio);
    arma::vec solution_sparse(size);
    // warm-up
    for (int i = 0; i < 5; i++) {
        solveSP(LHS_sparse, RHS, solution_sparse);
    }

    auto start_sparse = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < iter; i++) {
        solveSP(LHS_sparse, RHS, solution_sparse);
    }
    auto end_sparse = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_sparse = end_sparse - start_sparse;
    double sparse_time = elapsed_sparse.count();
    double sparse_time_per_solve = sparse_time / iter;

    std::cout << "Sparse solver time for size " << size << "x" << size << " over " << iter << " iterations: "
              << sparse_time << " seconds (" << sparse_time_per_solve * 1e6 << " μs per solve)" << std::endl;

    return 0;
}
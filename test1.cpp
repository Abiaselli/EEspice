#define ARMA_DONT_USE_WRAPPER
#define ARMA_USE_MKL_ALLOC
#define ARMA_USE_SUPERLU
#include <iostream>
#include <armadillo>
#include <klu.h>
#include <chrono>

// --- CONTEXT STRUCTURE ---
// Holds the KLU state that needs to persist across iterations
struct KLUContext {
    klu_l_common Common;
    klu_l_symbolic* Symbolic = nullptr;
    klu_l_numeric* Numeric = nullptr;

    bool is_initialized = false;
};

// --- 1. SETUP (Run Once) ---
// Analyzes the sparsity pattern. 
// NOTE: 'A' must have the same structure as the matrices used later.
int KLUSetup(arma::sp_mat& A, KLUContext& ctx) {
    // Cast Pointers
    auto* Ap = reinterpret_cast<SuiteSparse_long*>(const_cast<arma::uword*>(A.col_ptrs));
    auto* Ai = reinterpret_cast<SuiteSparse_long*>(const_cast<arma::uword*>(A.row_indices));

    // Defaults
    klu_l_defaults(&ctx.Common);

    // Symbolic Analysis (The expensive part we want to cache)
    ctx.Symbolic = klu_l_analyze(A.n_rows, Ap, Ai, &ctx.Common);
    
    if (!ctx.Symbolic) {
        std::cerr << "KLU Analyze failed." << std::endl;
        return 1;
    }
    
    ctx.is_initialized = true;
    return 0;
}

// --- 2. SOLVE (Run Many Times) ---
// Factors the new values and solves the system
int KLUSolve(arma::sp_mat& A, arma::vec& b, KLUContext& ctx) {
    if (!ctx.is_initialized) return 1;

    // Pointer casts (Zero-Copy)
    auto* Ap = reinterpret_cast<SuiteSparse_long*>(const_cast<arma::uword*>(A.col_ptrs));
    auto* Ai = reinterpret_cast<SuiteSparse_long*>(const_cast<arma::uword*>(A.row_indices));
    double* Ax = const_cast<double*>(A.values);

    bool factored = false;

    // STEP A: Try Refactoring (Fastest)
    // We can only refactor if we already have a Numeric object from a previous run
    if (ctx.Numeric) {
        // klu_l_refactor returns TRUE (1) on success, FALSE (0) on failure
        int success = klu_l_refactor(Ap, Ai, Ax, ctx.Symbolic, ctx.Numeric, &ctx.Common);
        
        if (success) {
            factored = true;
        } else {
            // Refactor failed (pivots unstable due to value changes).
            // We must free the old Numeric object before creating a new one.
            // Note: Common.status will likely contain KLU_SINGULAR or similar.
            klu_l_free_numeric(&ctx.Numeric, &ctx.Common);
            ctx.Numeric = nullptr; 
        }
    }

    // STEP B: Fallback to Full Factorization (Slower)
    // We do this if it's the first run OR if refactor failed
    if (!factored) {
        ctx.Numeric = klu_l_factor(Ap, Ai, Ax, ctx.Symbolic, &ctx.Common);
        
        if (!ctx.Numeric) {
            std::cerr << "KLU Factor failed (Singular Matrix?)" << std::endl;
            return 1;
        }
    }

    // STEP C: Solve
    // klu_l_solve uses the Numeric object (whether refactored or freshly factored)
    int solve_status = klu_l_solve(ctx.Symbolic, ctx.Numeric, A.n_rows, 1, b.memptr(), &ctx.Common);
    
    if (solve_status == 0) {
         std::cerr << "KLU Solve failed." << std::endl;
         return 1;
    }

    return 0;
}

// --- 3. CLEANUP (Run Once) ---
void KLUCleanup(KLUContext& ctx) {
    if (ctx.Symbolic) klu_l_free_symbolic(&ctx.Symbolic, &ctx.Common);
    if (ctx.Numeric) klu_l_free_numeric(&ctx.Numeric, &ctx.Common);
}

int main() {
    constexpr int iterations = 2021;
    arma::sp_mat A;
    // Load a sparse matrix from a file
    // Ensure this file exists, or create a random one for testing
    // A.sprandu(1000, 1000, 0.1); // Fallback if file missing
    A.load("lhs_arma.txt", arma::coord_ascii); 
    
    std::cout << "A nnz: " << A.n_nonzero << ", size: " << A.n_rows << "x" << A.n_cols << std::endl;

    arma::vec b;
    // b.randu(1000); // Fallback if file missing
    b.load("rhs.txt", arma::raw_ascii);

    double solve_time = 0.0;
    KLUContext ctx;

    // --- PRE-LOOP SETUP ---
    // We analyze the structure of A once here.
    // Crucial Assumption: The Newton iterations change the VALUES of A, but not the PATTERN.
    if (KLUSetup(A, ctx) != 0) return -1;

    // warm-up run
    for (int i = 0; i < 10; ++i) {
        arma::sp_mat A_copy = A; 
        arma::vec b_copy = b;    
        KLUSolve(A_copy, b_copy, ctx);
    }

    // Measured Loop
    for(int i = 0; i < iterations; ++i) {
        // Simulation of Newton Step:
        // In a real app, you would update A_copy's values here.
        // Since we just copy A, the structure is guaranteed identical.
        arma::sp_mat A_copy = A; 
        arma::vec b_copy = b;     

        auto start = std::chrono::high_resolution_clock::now();
        
        // Run solver using cached Symbolic context
        int status = KLUSolve(A_copy, b_copy, ctx);
        
        auto end = std::chrono::high_resolution_clock::now();
        solve_time += std::chrono::duration<double>(end - start).count();
    }

    // --- FINAL CLEANUP ---
    KLUCleanup(ctx);

    std::cout << " KLU solve time over " << iterations << " iterations: " << solve_time << " seconds." << std::endl;
    return 0;
}
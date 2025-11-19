#include <iostream>
#include <armadillo>
#include <klu.h>

// Helper concept to ensure we are actually on a 64-bit system
static_assert(sizeof(arma::uword) == 8, "This optimization requires 64-bit Armadillo indices");
static_assert(sizeof(SuiteSparse_long) == 8, "This optimization requires 64-bit SuiteSparse indices");

int main() {
    // 1. Setup Data (Same as before)
    arma::umat locations = {
        {0, 0}, {1, 0}, {0, 1}, {2, 1}, {4, 1},
        {1, 2}, {2, 2}, {3, 2}, {4, 2}, {2, 3}, {1, 4}, {4, 4}
    };
    arma::vec values = {2., 3., 3., -1., 4., 4., -3., 1., 2., 2., 6., 1.};
    arma::sp_mat A(locations.t(), values);
    arma::vec b = {8., 45., -3., 3., 19.};

    // --- ZERO-COPY KLU INTEGRATION ---

    // 2. Cast Pointers
    // Armadillo exposes 'const uword*'. KLU 'l' functions expect 'SuiteSparse_long*'.
    // We use const_cast to strip const (KLU doesn't modify Ap/Ai, but signature requires it)
    // We use reinterpret_cast to treat unsigned long as signed long.
    
    auto* Ap = reinterpret_cast<SuiteSparse_long*>(const_cast<arma::uword*>(A.col_ptrs));
    auto* Ai = reinterpret_cast<SuiteSparse_long*>(const_cast<arma::uword*>(A.row_indices));
    double* Ax = const_cast<double*>(A.values); // Values are double*, no cast needed

    // 3. Initialize KLU (Use klu_l_defaults for long)
    klu_l_common Common;
    klu_l_defaults(&Common);

    // 4. Analyze (Use klu_l_analyze)
    klu_l_symbolic* Symbolic = klu_l_analyze(A.n_rows, Ap, Ai, &Common);
    
    if (!Symbolic) {
        std::cerr << "KLU Analyze failed." << std::endl;
        return 1;
    }

    // 5. Factor (Use klu_l_factor)
    klu_l_numeric* Numeric = klu_l_factor(Ap, Ai, Ax, Symbolic, &Common);

    if (!Numeric) {
        std::cerr << "KLU Factor failed." << std::endl;
        klu_l_free_symbolic(&Symbolic, &Common);
        return 1;
    }

    // 6. Solve (Use klu_l_solve)
    // In-place solve on b
    klu_l_solve(Symbolic, Numeric, A.n_rows, 1, b.memptr(), &Common);

    // 7. Cleanup (Use klu_l_free...)
    klu_l_free_symbolic(&Symbolic, &Common);
    klu_l_free_numeric(&Numeric, &Common);

    // Result is in b
    b.print("Solution:");

    return 0;
}
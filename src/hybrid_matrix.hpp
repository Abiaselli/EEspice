#ifndef HYBRID_MATRIX_HPP
#define HYBRID_MATRIX_HPP

#include <armadillo>
#include <variant>
#include <concepts>
#include <type_traits>

/**
 * HybridMatrix: A runtime-selectable matrix wrapper that can hold either
 * dense (arma::mat) or sparse (arma::sp_mat) matrices, providing a unified
 * interface for matrix operations in circuit simulation.
 *
 * This implementation uses C++20 features (std::variant, concepts, if constexpr)
 * to provide zero-overhead abstraction for matrix operations.
 */
class HybridMatrix {
private:
    std::variant<arma::mat, arma::sp_mat> data;
    bool is_sparse_matrix;
    size_t n_rows;
    size_t n_cols;

public:
    /**
     * Constructor: Creates either a dense or sparse matrix based on use_sparse flag
     * @param rows Number of rows
     * @param cols Number of columns
     * @param use_sparse If true, creates sparse matrix; otherwise dense
     */
    HybridMatrix(size_t rows, size_t cols, bool use_sparse = false)
        : n_rows(rows), n_cols(cols), is_sparse_matrix(use_sparse) {
        if (use_sparse) {
            data = arma::sp_mat(rows, cols);
        } else {
            data = arma::mat(rows, cols, arma::fill::zeros);
        }
    }

    /**
     * Default constructor: Creates empty dense matrix
     */
    HybridMatrix() : HybridMatrix(0, 0, false) {}

    /**
     * Copy constructor
     */
    HybridMatrix(const HybridMatrix& other)
        : data(other.data), is_sparse_matrix(other.is_sparse_matrix),
          n_rows(other.n_rows), n_cols(other.n_cols) {}

    /**
     * Assignment operator
     */
    HybridMatrix& operator=(const HybridMatrix& other) {
        if (this != &other) {
            data = other.data;
            is_sparse_matrix = other.is_sparse_matrix;
            n_rows = other.n_rows;
            n_cols = other.n_cols;
        }
        return *this;
    }

    /**
     * Check if matrix is sparse
     */
    bool is_sparse() const { return is_sparse_matrix; }

    /**
     * Get number of rows
     */
    size_t rows() const { return n_rows; }

    /**
     * Get number of columns
     */
    size_t cols() const { return n_cols; }

    /**
     * Add a value to matrix element (stamping operation)
     * This is the primary interface for circuit element stamping
     */
    void add_stamp(size_t row, size_t col, double value) {
        if (row >= n_rows || col >= n_cols) {
            throw std::out_of_range("Matrix index out of bounds in add_stamp");
        }

        std::visit([row, col, value](auto& mat) {
            mat(row, col) += value;
        }, data);
    }

    /**
     * Set a matrix element to a specific value
     */
    void set_element(size_t row, size_t col, double value) {
        if (row >= n_rows || col >= n_cols) {
            throw std::out_of_range("Matrix index out of bounds in set_element");
        }

        std::visit([row, col, value](auto& mat) {
            mat(row, col) = value;
        }, data);
    }

    /**
     * Get a matrix element value
     */
    double get_element(size_t row, size_t col) const {
        if (row >= n_rows || col >= n_cols) {
            throw std::out_of_range("Matrix index out of bounds in get_element");
        }

        return std::visit([row, col](const auto& mat) -> double {
            return mat(row, col);
        }, data);
    }

    /**
     * Solve linear system Ax = b
     * Automatically selects appropriate solver based on matrix type
     */
    arma::vec solve(const arma::vec& rhs) const {
        if (rhs.n_elem != n_rows) {
            throw std::invalid_argument("RHS vector size mismatch in solve");
        }

        return std::visit([&rhs](const auto& mat) -> arma::vec {
            using MatType = std::decay_t<decltype(mat)>;

            if constexpr (std::is_same_v<MatType, arma::sp_mat>) {
                // Use sparse solver for sparse matrices
                // "superlu" is typically fastest for circuit matrices
                return arma::spsolve(mat, rhs, "superlu");
            } else {
                // Use dense solver for dense matrices
                return arma::solve(mat, rhs, arma::solve_opts::fast);
            }
        }, data);
    }

    /**
     * Solve linear system with options for different solver methods
     * @param rhs Right-hand side vector
     * @param method Solver method: "auto", "lu", "qr", "cholesky" for dense;
     *               "superlu", "lapack" for sparse
     */
    arma::vec solve_with_method(const arma::vec& rhs, const std::string& method = "auto") const {
        if (rhs.n_elem != n_rows) {
            throw std::invalid_argument("RHS vector size mismatch in solve_with_method");
        }

        return std::visit([&rhs, &method](const auto& mat) -> arma::vec {
            using MatType = std::decay_t<decltype(mat)>;

            if constexpr (std::is_same_v<MatType, arma::sp_mat>) {
                // Sparse solver options
                if (method == "auto" || method == "superlu") {
                    return arma::spsolve(mat, rhs, "superlu");
                } else if (method == "lapack") {
                    return arma::spsolve(mat, rhs, "lapack");
                } else {
                    // Default to superlu for unknown methods
                    return arma::spsolve(mat, rhs, "superlu");
                }
            } else {
                // Dense solver options
                if (method == "auto") {
                    return arma::solve(mat, rhs, arma::solve_opts::fast);
                } else if (method == "lu") {
                    return arma::solve(mat, rhs);
                } else if (method == "qr") {
                    arma::mat Q, R;
                    arma::qr(Q, R, mat);
                    return arma::solve(R, Q.t() * rhs);
                } else {
                    // Default to standard solve
                    return arma::solve(mat, rhs, arma::solve_opts::fast);
                }
            }
        }, data);
    }

    /**
     * Create a deep copy of the matrix
     */
    HybridMatrix copy() const {
        return HybridMatrix(*this);
    }

    /**
     * Reset this matrix to match another matrix
     * Used for resetting to initial state during iterations
     */
    void reset_to(const HybridMatrix& source) {
        *this = source;
    }

    /**
     * Fill matrix with zeros
     */
    void zeros() {
        std::visit([](auto& mat) {
            mat.zeros();
        }, data);
    }

    /**
     * Fill matrix with specific value (only for dense matrices)
     */
    void fill(double value) {
        std::visit([value](auto& mat) {
            if constexpr (std::is_same_v<std::decay_t<decltype(mat)>, arma::mat>) {
                mat.fill(value);
            } else {
                // For sparse matrices, filling with non-zero values is inefficient
                // Only allow filling with zero
                if (std::abs(value) < 1e-15) {
                    mat.zeros();
                } else {
                    throw std::runtime_error("Cannot fill sparse matrix with non-zero value");
                }
            }
        }, data);
    }

    /**
     * Convert to dense matrix (for compatibility with existing code)
     */
    arma::mat to_dense() const {
        return std::visit([](const auto& mat) -> arma::mat {
            using MatType = std::decay_t<decltype(mat)>;

            if constexpr (std::is_same_v<MatType, arma::sp_mat>) {
                return arma::mat(mat);
            } else {
                return mat;
            }
        }, data);
    }

    /**
     * Get reference to underlying dense matrix (throws if sparse)
     */
    arma::mat& get_dense() {
        if (is_sparse_matrix) {
            throw std::runtime_error("Cannot get dense reference from sparse matrix");
        }
        return std::get<arma::mat>(data);
    }

    /**
     * Get const reference to underlying dense matrix (throws if sparse)
     */
    const arma::mat& get_dense() const {
        if (is_sparse_matrix) {
            throw std::runtime_error("Cannot get dense reference from sparse matrix");
        }
        return std::get<arma::mat>(data);
    }

    /**
     * Get reference to underlying sparse matrix (throws if dense)
     */
    arma::sp_mat& get_sparse() {
        if (!is_sparse_matrix) {
            throw std::runtime_error("Cannot get sparse reference from dense matrix");
        }
        return std::get<arma::sp_mat>(data);
    }

    /**
     * Get const reference to underlying sparse matrix (throws if dense)
     */
    const arma::sp_mat& get_sparse() const {
        if (!is_sparse_matrix) {
            throw std::runtime_error("Cannot get sparse reference from dense matrix");
        }
        return std::get<arma::sp_mat>(data);
    }

    /**
     * Resize matrix (preserves existing values where possible)
     */
    void resize(size_t new_rows, size_t new_cols) {
        n_rows = new_rows;
        n_cols = new_cols;

        std::visit([new_rows, new_cols](auto& mat) {
            mat.resize(new_rows, new_cols);
        }, data);
    }

    /**
     * Apply a visitor to the underlying matrix
     * Useful for operations not covered by the standard interface
     */
    template<typename Visitor>
    auto visit(Visitor&& vis) {
        return std::visit(std::forward<Visitor>(vis), data);
    }

    /**
     * Apply a const visitor to the underlying matrix
     */
    template<typename Visitor>
    auto visit(Visitor&& vis) const {
        return std::visit(std::forward<Visitor>(vis), data);
    }

    /**
     * Get an estimate of memory usage
     */
    size_t memory_usage() const {
        return std::visit([](const auto& mat) -> size_t {
            using MatType = std::decay_t<decltype(mat)>;

            if constexpr (std::is_same_v<MatType, arma::sp_mat>) {
                // Sparse matrix memory: values + row indices + column pointers
                return mat.n_nonzero * sizeof(double) +
                       mat.n_nonzero * sizeof(arma::uword) +
                       (mat.n_cols + 1) * sizeof(arma::uword);
            } else {
                // Dense matrix memory
                return mat.n_elem * sizeof(double);
            }
        }, data);
    }

    /**
     * Get number of non-zero elements
     */
    size_t n_nonzero() const {
        return std::visit([](const auto& mat) -> size_t {
            using MatType = std::decay_t<decltype(mat)>;

            if constexpr (std::is_same_v<MatType, arma::sp_mat>) {
                return mat.n_nonzero;
            } else {
                // Count non-zeros in dense matrix
                return arma::accu(mat != 0.0);
            }
        }, data);
    }

    /**
     * Get sparsity ratio (percentage of zero elements)
     */
    double sparsity() const {
        if (n_rows == 0 || n_cols == 0) return 0.0;
        size_t total_elements = n_rows * n_cols;
        size_t nonzeros = n_nonzero();
        return 100.0 * (1.0 - static_cast<double>(nonzeros) / total_elements);
    }

    /**
     * Print matrix information (for debugging)
     */
    void print_info() const {
        std::cout << "HybridMatrix: "
                  << (is_sparse_matrix ? "Sparse" : "Dense")
                  << " [" << n_rows << " x " << n_cols << "]"
                  << " Non-zeros: " << n_nonzero()
                  << " Sparsity: " << sparsity() << "%"
                  << " Memory: " << memory_usage() / 1024 << " KB"
                  << std::endl;
    }
};

#endif // HYBRID_MATRIX_HPP
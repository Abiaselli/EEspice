#ifndef HYBRID_MATRIX_HPP
#define HYBRID_MATRIX_HPP

#include <armadillo>
#include <variant>
#include <concepts>
#include <type_traits>
#include <unordered_map>
#include <vector>
#include <cstring>
#include <memory>

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

    // Pre-allocated sparse structure support
    // Maps (row<<32|col) to index in CSC values array for O(1) access
    // Shared pointer to avoid expensive copies during NR iterations
    std::shared_ptr<std::unordered_map<uint64_t, size_t>> position_map;
    bool pattern_locked = false;

    // Two-level baseline support for efficient NR iteration reset
    // For sparse matrices: store values arrays for memcpy
    std::vector<double> baseline_linear;   // State after stamping linear elements
    std::vector<double> baseline_step;     // State after stamping dynamic elements (capacitors)
    // For dense matrices: store full matrix copies
    arma::mat baseline_linear_dense;       // Dense matrix linear baseline
    arma::mat baseline_step_dense;         // Dense matrix step baseline

    // Explicit validity flags - don't infer from buffer size
    bool linear_baseline_valid = false;
    bool step_baseline_valid = false;

    static uint64_t make_key(size_t row, size_t col) {
        return (static_cast<uint64_t>(row) << 32) | static_cast<uint64_t>(col);
    }

    // Initialize position_map if needed (lazy initialization)
    void ensure_position_map() {
        if (!position_map) {
            position_map = std::make_shared<std::unordered_map<uint64_t, size_t>>();
        }
    }

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
     * Note: position_map is shared (not deep copied) since pattern never changes after locking
     */
    HybridMatrix(const HybridMatrix& other)
        : data(other.data), is_sparse_matrix(other.is_sparse_matrix),
          n_rows(other.n_rows), n_cols(other.n_cols),
          position_map(other.position_map),  // Shared pointer - just copies the pointer
          pattern_locked(other.pattern_locked),
          baseline_linear(other.baseline_linear),
          baseline_step(other.baseline_step),
          baseline_linear_dense(other.baseline_linear_dense),
          baseline_step_dense(other.baseline_step_dense),
          linear_baseline_valid(other.linear_baseline_valid),
          step_baseline_valid(other.step_baseline_valid) {}

    /**
     * Assignment operator
     */
    HybridMatrix& operator=(const HybridMatrix& other) {
        if (this != &other) {
            data = other.data;
            is_sparse_matrix = other.is_sparse_matrix;
            n_rows = other.n_rows;
            n_cols = other.n_cols;
            position_map = other.position_map;
            pattern_locked = other.pattern_locked;
            baseline_linear = other.baseline_linear;
            baseline_step = other.baseline_step;
            baseline_linear_dense = other.baseline_linear_dense;
            baseline_step_dense = other.baseline_step_dense;
            linear_baseline_valid = other.linear_baseline_valid;
            step_baseline_valid = other.step_baseline_valid;
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
     * For locked sparse matrices, preserves CSC structure
     */
    void zeros() {
        if (is_sparse_matrix && pattern_locked) {
            zero_values_keep_pattern();
            return;
        }
        std::visit([](auto& mat) {
            mat.zeros();
        }, data);
    }

    /**
     * Fill matrix with specific value (only for dense matrices)
     */
    void fill(double value) {
        std::visit([this, value](auto& mat) {
            if constexpr (std::is_same_v<std::decay_t<decltype(mat)>, arma::mat>) {
                mat.fill(value);
            } else {
                // For sparse matrices, filling with non-zero values is inefficient
                // Only allow filling with zero
                if (std::abs(value) < 1e-15) {
                    if (pattern_locked) {
                        zero_values_keep_pattern();
                    } else {
                        mat.zeros();
                    }
                } else {
                    throw std::runtime_error("Cannot fill sparse matrix with non-zero value");
                }
            }
        }, data);
    }

    /**
     * Zero all values while preserving sparse structure (CSC arrays unchanged)
     * For dense matrices, equivalent to zeros()
     * For locked sparse matrices, only clears values without destroying pattern
     */
    void zero_values_keep_pattern() {
        if (!is_sparse_matrix) {
            std::get<arma::mat>(data).zeros();
            return;
        }

        arma::sp_mat& sp = std::get<arma::sp_mat>(data);
        sp.sync();  // Ensure CSC is materialized
        double* vals = arma::access::rwp(sp.values);
        std::fill(vals, vals + sp.n_nonzero, 0.0);
    }

    /**
     * Convert internal representation to dense matrix (in-place conversion)
     * If already dense, this is a no-op. If sparse, converts to dense and updates flag.
     */
    void to_dense() {
        if (is_sparse_matrix) {
            // Convert sparse to dense in-place
            arma::sp_mat& sparse_mat = std::get<arma::sp_mat>(data);
            arma::mat dense_mat(sparse_mat);
            data = std::move(dense_mat);
            is_sparse_matrix = false;
        }
        // If already dense, do nothing
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

    // =========================================================================
    // Pre-allocated Sparse Structure API
    // =========================================================================
    // These methods enable O(1) element access for sparse matrices by:
    // 1. Recording all positions that will be stamped (pattern discovery)
    // 2. Building CSC structure once with those positions
    // 3. Using direct pointer access to values array during stamping
    // =========================================================================

    /**
     * Check if pattern has been locked
     */
    bool is_pattern_locked() const { return pattern_locked; }

    /**
     * Phase 1: Record a position that will receive stamps
     * Call this during setup phase to discover sparsity pattern.
     * No-op if pattern is already locked or matrix is dense.
     */
    void record_position(size_t row, size_t col) {
        if (pattern_locked || !is_sparse_matrix) return;
        if (row >= n_rows || col >= n_cols) {
            return;
        }

        ensure_position_map();
        uint64_t key = make_key(row, col);
        if (position_map->find(key) == position_map->end()) {
            (*position_map)[key] = position_map->size();  // Temporary index
        }
    }

    /**
     * Phase 2: Lock the pattern and build CSC structure
     * After calling this, no new positions can be recorded.
     * The sparse matrix is rebuilt with ALL positions (both existing and recorded),
     * preserving existing values. We use a tiny placeholder value to ensure Armadillo
     * doesn't filter out positions with zero values.
     */
    void lock_pattern() {
        if (pattern_locked || !is_sparse_matrix) return;

        ensure_position_map();

        // Get reference to current sparse matrix
        arma::sp_mat& old_sp = std::get<arma::sp_mat>(data);

        // Merge existing non-zero positions with recorded positions
        // First, add all existing non-zeros to position_map (if not already there)
        for (arma::sp_mat::const_iterator it = old_sp.begin(); it != old_sp.end(); ++it) {
            uint64_t key = make_key(it.row(), it.col());
            if (position_map->find(key) == position_map->end()) {
                (*position_map)[key] = position_map->size();
            }
        }

        if (position_map->empty()) {
            pattern_locked = true;
            return;
        }

        // Build triplets from all positions, preserving existing values
        // Use a tiny placeholder (1e-300) for positions with zero values
        // to prevent Armadillo from filtering them out during construction
        size_t nnz = position_map->size();
        arma::umat locations(2, nnz);
        arma::vec values(nnz);

        size_t idx = 0;
        for (const auto& [key, _] : *position_map) {
            arma::uword row = static_cast<arma::uword>(key >> 32);
            arma::uword col = static_cast<arma::uword>(key & 0xFFFFFFFF);
            locations(0, idx) = row;
            locations(1, idx) = col;
            // Preserve existing value, or use tiny placeholder if zero
            double val = old_sp(row, col);
            values(idx) = (val != 0.0) ? val : 1e-300;
            idx++;
        }

        // Create sparse matrix with structure (sorted by Armadillo)
        data = arma::sp_mat(locations, values, n_rows, n_cols);

        // Now set the placeholder values back to zero using direct access
        arma::sp_mat& sp = std::get<arma::sp_mat>(data);
        double* vals = arma::access::rwp(sp.values);
        for (size_t i = 0; i < sp.n_nonzero; i++) {
            if (std::abs(vals[i]) < 1e-290) {
                vals[i] = 0.0;
            }
        }

        // Rebuild position_map with actual CSC indices
        position_map->clear();
        for (arma::uword col = 0; col < sp.n_cols; col++) {
            for (arma::uword i = sp.col_ptrs[col]; i < sp.col_ptrs[col + 1]; i++) {
                arma::uword row = sp.row_indices[i];
                (*position_map)[make_key(row, col)] = i;
            }
        }

        pattern_locked = true;
        baseline_linear.resize(sp.n_nonzero, 0.0);
    }

    /**
     * Phase 3: Fast indexed stamping using pre-computed positions
     * Falls back to regular add_stamp if pattern not locked or dense.
     */
    void add_stamp_indexed(size_t row, size_t col, double value) {
        if (!pattern_locked || !is_sparse_matrix || !position_map) {
            add_stamp(row, col, value);
            return;
        }

        uint64_t key = make_key(row, col);
        auto it = position_map->find(key);
        if (it != position_map->end()) {
            arma::sp_mat& sp = std::get<arma::sp_mat>(data);
            double* vals = arma::access::rwp(sp.values);
            vals[it->second] += value;
        }
        // Silently ignore positions not in the pattern (shouldn't happen if setup is correct)
    }

    // =========================================================================
    // Two-Level Baseline API for Newton-Raphson Reset Optimization
    // =========================================================================
    // Level 1 (Linear baseline): State after stamping static linear elements
    // Level 2 (Step baseline): State after stamping dynamic elements (capacitors)
    // =========================================================================

    /**
     * Save current values as linear baseline
     * Call after linear elements have been stamped and pattern is locked.
     */
    void save_linear_baseline() {
        if (is_sparse_matrix && pattern_locked) {
            arma::sp_mat& sp = std::get<arma::sp_mat>(data);
            baseline_linear.resize(sp.n_nonzero);
            std::memcpy(baseline_linear.data(), sp.values, sp.n_nonzero * sizeof(double));
        } else if (!is_sparse_matrix) {
            baseline_linear_dense = std::get<arma::mat>(data);
        }
        linear_baseline_valid = true;
    }

    /**
     * Reset matrix to linear baseline (state after linear elements stamped)
     */
    void reset_to_linear_baseline() {
        if (is_sparse_matrix && pattern_locked && linear_baseline_valid) {
            arma::sp_mat& sp = std::get<arma::sp_mat>(data);
            std::memcpy(arma::access::rwp(sp.values), baseline_linear.data(),
                        sp.n_nonzero * sizeof(double));
        } else if (!is_sparse_matrix && linear_baseline_valid) {
            std::get<arma::mat>(data) = baseline_linear_dense;
        }
    }

    /**
     * Save current values as step baseline
     * Call after stamping dynamic elements (capacitors) at start of each time step.
     */
    void save_step_baseline() {
        if (is_sparse_matrix && pattern_locked) {
            arma::sp_mat& sp = std::get<arma::sp_mat>(data);
            baseline_step.resize(sp.n_nonzero);
            std::memcpy(baseline_step.data(), sp.values, sp.n_nonzero * sizeof(double));
        } else if (!is_sparse_matrix) {
            baseline_step_dense = std::get<arma::mat>(data);
        }
        step_baseline_valid = true;
    }

    /**
     * Reset matrix to step baseline (state after dynamic elements stamped)
     * Used at the start of each NR iteration.
     */
    void reset_to_step_baseline() {
        if (is_sparse_matrix && pattern_locked && step_baseline_valid) {
            arma::sp_mat& sp = std::get<arma::sp_mat>(data);
            std::memcpy(arma::access::rwp(sp.values), baseline_step.data(),
                        sp.n_nonzero * sizeof(double));
        } else if (!is_sparse_matrix && step_baseline_valid) {
            std::get<arma::mat>(data) = baseline_step_dense;
        }
    }

    /**
     * Initialize step baseline from linear baseline
     * Used for DC analysis where step baseline == linear baseline
     */
    void copy_linear_to_step_baseline() {
        if (is_sparse_matrix) {
            baseline_step = baseline_linear;
        } else {
            baseline_step_dense = baseline_linear_dense;
        }
        step_baseline_valid = linear_baseline_valid;
    }

    /**
     * Check if linear baseline has been saved
     */
    bool has_linear_baseline() const {
        return linear_baseline_valid;
    }

    /**
     * Check if step baseline has been saved
     */
    bool has_step_baseline() const {
        return step_baseline_valid;
    }

    /**
     * Get number of positions in the pattern
     */
    size_t pattern_size() const {
        return position_map ? position_map->size() : 0;
    }
};

#endif // HYBRID_MATRIX_HPP
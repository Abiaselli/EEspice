#include<iostream>
#include<armadillo>
#include "src/models.hpp"
#include "src/device.hpp"

// Helper function to compare two matrices and report differences
bool compare_matrices(const arma::mat& hybrid_result, const arma::mat& arma_result,
                     const std::string& test_name, bool verbose = false) {
    if (hybrid_result.n_rows != arma_result.n_rows || hybrid_result.n_cols != arma_result.n_cols) {
        std::cout << "L FAIL: " << test_name << " - Matrix dimensions mismatch\n";
        std::cout << "  HybridMatrix: " << hybrid_result.n_rows << "x" << hybrid_result.n_cols << "\n";
        std::cout << "  arma::mat: " << arma_result.n_rows << "x" << arma_result.n_cols << "\n";
        return false;
    }

    bool all_equal = arma::approx_equal(hybrid_result, arma_result, "absdiff", 1e-15);

    if (all_equal) {
        std::cout << " PASS: " << test_name << "\n";
        return true;
    } else {
        std::cout << "L FAIL: " << test_name << "\n";
        if (verbose) {
            std::cout << "  HybridMatrix result:\n" << hybrid_result << "\n";
            std::cout << "  arma::mat result:\n" << arma_result << "\n";
            std::cout << "  Difference:\n" << (hybrid_result - arma_result) << "\n";
        }
        // Always show which elements differ
        for (size_t i = 0; i < hybrid_result.n_rows; i++) {
            for (size_t j = 0; j < hybrid_result.n_cols; j++) {
                double diff = std::abs(hybrid_result(i, j) - arma_result(i, j));
                if (diff > 1e-15) {
                    std::cout << "  Mismatch at (" << i << "," << j << "): "
                              << hybrid_result(i, j) << " vs " << arma_result(i, j)
                              << " (diff: " << diff << ")\n";
                }
            }
        }
        return false;
    }
}

// Helper function to compare both LHS matrix and RHS vector for Vs_assigner
bool compare_vs_results(const arma::mat& hybrid_lhs, const arma::mat& arma_lhs,
                       const arma::vec& hybrid_rhs, const arma::vec& arma_rhs,
                       const std::string& test_name, bool verbose = false) {
    bool lhs_match = true;
    bool rhs_match = true;

    // Check LHS matrix dimensions
    if (hybrid_lhs.n_rows != arma_lhs.n_rows || hybrid_lhs.n_cols != arma_lhs.n_cols) {
        std::cout << "L FAIL: " << test_name << " - LHS matrix dimensions mismatch\n";
        std::cout << "  HybridMatrix: " << hybrid_lhs.n_rows << "x" << hybrid_lhs.n_cols << "\n";
        std::cout << "  arma::mat: " << arma_lhs.n_rows << "x" << arma_lhs.n_cols << "\n";
        lhs_match = false;
    } else if (!arma::approx_equal(hybrid_lhs, arma_lhs, "absdiff", 1e-15)) {
        std::cout << "L FAIL: " << test_name << " - LHS matrix content mismatch\n";
        if (verbose) {
            std::cout << "  HybridMatrix LHS:\n" << hybrid_lhs << "\n";
            std::cout << "  arma::mat LHS:\n" << arma_lhs << "\n";
        }
        for (size_t i = 0; i < hybrid_lhs.n_rows; i++) {
            for (size_t j = 0; j < hybrid_lhs.n_cols; j++) {
                double diff = std::abs(hybrid_lhs(i, j) - arma_lhs(i, j));
                if (diff > 1e-15) {
                    std::cout << "  LHS mismatch at (" << i << "," << j << "): "
                              << hybrid_lhs(i, j) << " vs " << arma_lhs(i, j)
                              << " (diff: " << diff << ")\n";
                }
            }
        }
        lhs_match = false;
    }

    // Check RHS vector dimensions
    if (hybrid_rhs.n_elem != arma_rhs.n_elem) {
        std::cout << "L FAIL: " << test_name << " - RHS vector dimensions mismatch\n";
        std::cout << "  HybridMatrix RHS size: " << hybrid_rhs.n_elem << "\n";
        std::cout << "  arma::mat RHS size: " << arma_rhs.n_elem << "\n";
        rhs_match = false;
    } else if (!arma::approx_equal(hybrid_rhs, arma_rhs, "absdiff", 1e-15)) {
        std::cout << "L FAIL: " << test_name << " - RHS vector content mismatch\n";
        if (verbose) {
            std::cout << "  HybridMatrix RHS:\n" << hybrid_rhs.t() << "\n";
            std::cout << "  arma::mat RHS:\n" << arma_rhs.t() << "\n";
        }
        for (size_t i = 0; i < hybrid_rhs.n_elem; i++) {
            double diff = std::abs(hybrid_rhs(i) - arma_rhs(i));
            if (diff > 1e-15) {
                std::cout << "  RHS mismatch at index " << i << ": "
                          << hybrid_rhs(i) << " vs " << arma_rhs(i)
                          << " (diff: " << diff << ")\n";
            }
        }
        rhs_match = false;
    }

    if (lhs_match && rhs_match) {
        std::cout << " PASS: " << test_name << "\n";
        return true;
    }
    return false;
}

// Test case structure
struct TestCase {
    std::string name;
    int node_x;
    int node_y;
    double G;
};

void test_R_assigner() {
    std::cout << "=================================================\n";
    std::cout << "====== Testing R_assigner Function Equivalence =====\n";
    std::cout << "=================================================\n\n";

    int total_tests = 0;
    int passed_tests = 0;

    // Test Case 1: Both nodes are ground (should do nothing)
    {
        std::cout << "Test Case 1: Both nodes ground (node_x=0, node_y=0)\n";
        const size_t size = 5;
        HybridMatrix hybrid_lhs(size, size, false);  // Dense mode
        arma::mat arma_lhs = arma::zeros(size, size);

        R_assigner(0, 0, 1.0, hybrid_lhs);
        R_assigner(0, 0, 1.0, arma_lhs);

        total_tests++;
        hybrid_lhs.to_dense();
        if (compare_matrices(hybrid_lhs.get_dense(), arma_lhs, "Ground-to-ground resistor")) {
            passed_tests++;
        }
        std::cout << "\n";
    }

    // Test Case 2: node_x = 0, node_y != 0 (single diagonal stamp)
    {
        std::cout << "Test Case 2: node_x ground, node_y non-zero (node_x=0, node_y=2)\n";
        const size_t size = 5;
        HybridMatrix hybrid_lhs(size, size, false);
        arma::mat arma_lhs = arma::zeros(size, size);

        R_assigner(0, 2, 0.5, hybrid_lhs);
        R_assigner(0, 2, 0.5, arma_lhs);

        total_tests++;
        hybrid_lhs.to_dense();
        if (compare_matrices(hybrid_lhs.get_dense(), arma_lhs, "Resistor from ground to node 2")) {
            passed_tests++;
        }
        std::cout << "\n";
    }

    // Test Case 3: node_x != 0, node_y = 0 (single diagonal stamp)
    {
        std::cout << "Test Case 3: node_x non-zero, node_y ground (node_x=3, node_y=0)\n";
        const size_t size = 5;
        HybridMatrix hybrid_lhs(size, size, false);
        arma::mat arma_lhs = arma::zeros(size, size);

        R_assigner(3, 0, 2.0, hybrid_lhs);
        R_assigner(3, 0, 2.0, arma_lhs);

        total_tests++;
        hybrid_lhs.to_dense();
        if (compare_matrices(hybrid_lhs.get_dense(), arma_lhs, "Resistor from node 3 to ground")) {
            passed_tests++;
        }
        std::cout << "\n";
    }

    // Test Case 4: Both nodes non-zero (2x2 conductance pattern)
    {
        std::cout << "Test Case 4: Both nodes non-zero (node_x=2, node_y=4)\n";
        const size_t size = 5;
        HybridMatrix hybrid_lhs(size, size, false);
        arma::mat arma_lhs = arma::zeros(size, size);

        R_assigner(2, 4, 0.1, hybrid_lhs);
        R_assigner(2, 4, 0.1, arma_lhs);

        total_tests++;
        hybrid_lhs.to_dense();
        if (compare_matrices(hybrid_lhs.get_dense(), arma_lhs, "Resistor between node 2 and node 4")) {
            passed_tests++;
        }
        std::cout << "\n";
    }

    // Test Case 5: Multiple resistors (verify accumulation works correctly)
    {
        std::cout << "Test Case 5: Multiple resistors with accumulation\n";
        const size_t size = 5;
        HybridMatrix hybrid_lhs(size, size, false);
        arma::mat arma_lhs = arma::zeros(size, size);

        // Add several resistors
        R_assigner(1, 2, 1.0, hybrid_lhs);
        R_assigner(1, 2, 1.0, arma_lhs);

        R_assigner(2, 3, 0.5, hybrid_lhs);
        R_assigner(2, 3, 0.5, arma_lhs);

        R_assigner(1, 0, 0.25, hybrid_lhs);
        R_assigner(1, 0, 0.25, arma_lhs);

        R_assigner(3, 0, 2.0, hybrid_lhs);
        R_assigner(3, 0, 2.0, arma_lhs);

        total_tests++;
        hybrid_lhs.to_dense();
        if (compare_matrices(hybrid_lhs.get_dense(), arma_lhs, "Multiple resistors with accumulation")) {
            passed_tests++;
        }
        std::cout << "\n";
    }

    // Test Case 6: Large conductance value
    {
        std::cout << "Test Case 6: Large conductance value\n";
        const size_t size = 5;
        HybridMatrix hybrid_lhs(size, size, false);
        arma::mat arma_lhs = arma::zeros(size, size);

        R_assigner(1, 3, 1e6, hybrid_lhs);
        R_assigner(1, 3, 1e6, arma_lhs);

        total_tests++;
        hybrid_lhs.to_dense();
        if (compare_matrices(hybrid_lhs.get_dense(), arma_lhs, "Large conductance value (1e6)")) {
            passed_tests++;
        }
        std::cout << "\n";
    }

    // Test Case 7: Small conductance value
    {
        std::cout << "Test Case 7: Small conductance value\n";
        const size_t size = 5;
        HybridMatrix hybrid_lhs(size, size, false);
        arma::mat arma_lhs = arma::zeros(size, size);

        R_assigner(2, 5, 1e-12, hybrid_lhs);
        R_assigner(2, 5, 1e-12, arma_lhs);

        total_tests++;
        hybrid_lhs.to_dense();
        if (compare_matrices(hybrid_lhs.get_dense(), arma_lhs, "Small conductance value (1e-12)")) {
            passed_tests++;
        }
        std::cout << "\n";
    }

    // Test Case 8: Negative conductance (for testing edge cases)
    {
        std::cout << "Test Case 8: Negative conductance\n";
        const size_t size = 5;
        HybridMatrix hybrid_lhs(size, size, false);
        arma::mat arma_lhs = arma::zeros(size, size);

        R_assigner(1, 4, -0.5, hybrid_lhs);
        R_assigner(1, 4, -0.5, arma_lhs);

        total_tests++;
        hybrid_lhs.to_dense();
        if (compare_matrices(hybrid_lhs.get_dense(), arma_lhs, "Negative conductance")) {
            passed_tests++;
        }
        std::cout << "\n";
    }

    // Test Case 9: Adjacent nodes
    {
        std::cout << "Test Case 9: Adjacent nodes (node_x=1, node_y=2)\n";
        const size_t size = 5;
        HybridMatrix hybrid_lhs(size, size, false);
        arma::mat arma_lhs = arma::zeros(size, size);

        R_assigner(1, 2, 0.333, hybrid_lhs);
        R_assigner(1, 2, 0.333, arma_lhs);

        total_tests++;
        hybrid_lhs.to_dense();
        if (compare_matrices(hybrid_lhs.get_dense(), arma_lhs, "Adjacent nodes resistor")) {
            passed_tests++;
        }
        std::cout << "\n";
    }

    // Test Case 10: Complex circuit with many resistors
    {
        std::cout << "Test Case 10: Complex circuit with 10+ resistors\n";
        const size_t size = 10;
        HybridMatrix hybrid_lhs(size, size, false);
        arma::mat arma_lhs = arma::zeros(size, size);

        // Create a more complex resistor network
        std::vector<TestCase> resistors = {
            {"R1", 1, 2, 1.0},
            {"R2", 2, 3, 0.5},
            {"R3", 3, 4, 2.0},
            {"R4", 1, 0, 0.1},
            {"R5", 4, 0, 0.2},
            {"R6", 2, 5, 1.5},
            {"R7", 5, 6, 0.8},
            {"R8", 6, 0, 0.3},
            {"R9", 3, 7, 1.2},
            {"R10", 7, 8, 0.6},
            {"R11", 8, 9, 0.4},
            {"R12", 9, 0, 0.5},
            {"R13", 1, 5, 0.7},
            {"R14", 3, 6, 0.9}
        };

        for (const auto& r : resistors) {
            R_assigner(r.node_x, r.node_y, r.G, hybrid_lhs);
            R_assigner(r.node_x, r.node_y, r.G, arma_lhs);
        }

        total_tests++;
        hybrid_lhs.to_dense();
        if (compare_matrices(hybrid_lhs.get_dense(), arma_lhs, "Complex circuit with 14 resistors")) {
            passed_tests++;
        }
        std::cout << "\n";
    }

    // Summary
    std::cout << "=========================================\n";
    std::cout << "R_assigner RESULTS: " << passed_tests << "/" << total_tests << " tests passed\n";
    if (passed_tests == total_tests) {
        std::cout << " ALL R_assigner TESTS PASSED\n";
    } else {
        std::cout << "L SOME R_assigner TESTS FAILED\n";
    }
    std::cout << "=========================================\n\n";
}

void test_Vs_assigner() {
    std::cout << "=================================================\n";
    std::cout << "====== Testing Vs_assigner Function Equivalence =====\n";
    std::cout << "=================================================\n\n";

    int total_tests = 0;
    int passed_tests = 0;

    // Test Case 1: Voltage source between two non-ground nodes
    {
        std::cout << "Test Case 1: V-source between non-ground nodes (node_x=2, node_y=3, V=5.0)\n";
        const size_t size = 3;
        HybridMatrix hybrid_lhs(size, size, false);
        arma::mat arma_lhs = arma::zeros(size, size);
        arma::vec hybrid_rhs = arma::zeros(size);
        arma::vec arma_rhs = arma::zeros(size);

        Vs_assigner(2, 3, 5.0, hybrid_lhs, hybrid_rhs);
        Vs_assigner(2, 3, 5.0, arma_lhs, arma_rhs);

        total_tests++;
        hybrid_lhs.to_dense();
        if (compare_vs_results(hybrid_lhs.get_dense(), arma_lhs, hybrid_rhs, arma_rhs,
                              "V-source between nodes 2 and 3")) {
            passed_tests++;
        }
        std::cout << "\n";
    }

    // Test Case 2: Voltage source from node to ground
    {
        std::cout << "Test Case 2: V-source from node to ground (node_x=2, node_y=0, V=3.3)\n";
        const size_t size = 3;
        HybridMatrix hybrid_lhs(size, size, false);
        arma::mat arma_lhs = arma::zeros(size, size);
        arma::vec hybrid_rhs = arma::zeros(size);
        arma::vec arma_rhs = arma::zeros(size);

        Vs_assigner(2, 0, 3.3, hybrid_lhs, hybrid_rhs);
        Vs_assigner(2, 0, 3.3, arma_lhs, arma_rhs);

        total_tests++;
        hybrid_lhs.to_dense();
        if (compare_vs_results(hybrid_lhs.get_dense(), arma_lhs, hybrid_rhs, arma_rhs,
                              "V-source from node 2 to ground")) {
            passed_tests++;
        }
        std::cout << "\n";
    }

    // Test Case 3: Voltage source from ground to node
    {
        std::cout << "Test Case 3: V-source from ground to node (node_x=0, node_y=2, V=1.8)\n";
        const size_t size = 3;
        HybridMatrix hybrid_lhs(size, size, false);
        arma::mat arma_lhs = arma::zeros(size, size);
        arma::vec hybrid_rhs = arma::zeros(size);
        arma::vec arma_rhs = arma::zeros(size);

        Vs_assigner(0, 2, 1.8, hybrid_lhs, hybrid_rhs);
        Vs_assigner(0, 2, 1.8, arma_lhs, arma_rhs);

        total_tests++;
        hybrid_lhs.to_dense();
        if (compare_vs_results(hybrid_lhs.get_dense(), arma_lhs, hybrid_rhs, arma_rhs,
                              "V-source from ground to node 2")) {
            passed_tests++;
        }
        std::cout << "\n";
    }

    // Test Case 4: Multiple voltage sources (matrix grows: 3x3 -> 4x4 -> 5x5)
    {
        std::cout << "Test Case 4: Multiple voltage sources (3x3 -> 4x4 -> 5x5)\n";
        const size_t size = 3;
        HybridMatrix hybrid_lhs(size, size, false);
        arma::mat arma_lhs = arma::zeros(size, size);
        arma::vec hybrid_rhs = arma::zeros(size);
        arma::vec arma_rhs = arma::zeros(size);

        // First voltage source
        Vs_assigner(1, 0, 5.0, hybrid_lhs, hybrid_rhs);
        Vs_assigner(1, 0, 5.0, arma_lhs, arma_rhs);

        // Second voltage source
        Vs_assigner(2, 3, 3.3, hybrid_lhs, hybrid_rhs);
        Vs_assigner(2, 3, 3.3, arma_lhs, arma_rhs);

        total_tests++;
        hybrid_lhs.to_dense();
        if (compare_vs_results(hybrid_lhs.get_dense(), arma_lhs, hybrid_rhs, arma_rhs,
                              "Multiple voltage sources")) {
            passed_tests++;
        }
        std::cout << "\n";
    }

    // Test Case 5: Negative voltage value
    {
        std::cout << "Test Case 5: Negative voltage value (V=-5.0)\n";
        const size_t size = 3;
        HybridMatrix hybrid_lhs(size, size, false);
        arma::mat arma_lhs = arma::zeros(size, size);
        arma::vec hybrid_rhs = arma::zeros(size);
        arma::vec arma_rhs = arma::zeros(size);

        Vs_assigner(2, 1, -5.0, hybrid_lhs, hybrid_rhs);
        Vs_assigner(2, 1, -5.0, arma_lhs, arma_rhs);

        total_tests++;
        hybrid_lhs.to_dense();
        if (compare_vs_results(hybrid_lhs.get_dense(), arma_lhs, hybrid_rhs, arma_rhs,
                              "Negative voltage value")) {
            passed_tests++;
        }
        std::cout << "\n";
    }

    // Test Case 6: Zero voltage value
    {
        std::cout << "Test Case 6: Zero voltage value (V=0.0)\n";
        const size_t size = 3;
        HybridMatrix hybrid_lhs(size, size, false);
        arma::mat arma_lhs = arma::zeros(size, size);
        arma::vec hybrid_rhs = arma::zeros(size);
        arma::vec arma_rhs = arma::zeros(size);

        Vs_assigner(1, 2, 0.0, hybrid_lhs, hybrid_rhs);
        Vs_assigner(1, 2, 0.0, arma_lhs, arma_rhs);

        total_tests++;
        hybrid_lhs.to_dense();
        if (compare_vs_results(hybrid_lhs.get_dense(), arma_lhs, hybrid_rhs, arma_rhs,
                              "Zero voltage value")) {
            passed_tests++;
        }
        std::cout << "\n";
    }

    // Test Case 7: Large voltage value
    {
        std::cout << "Test Case 7: Large voltage value (V=1e6)\n";
        const size_t size = 3;
        HybridMatrix hybrid_lhs(size, size, false);
        arma::mat arma_lhs = arma::zeros(size, size);
        arma::vec hybrid_rhs = arma::zeros(size);
        arma::vec arma_rhs = arma::zeros(size);

        Vs_assigner(1, 3, 1e6, hybrid_lhs, hybrid_rhs);
        Vs_assigner(1, 3, 1e6, arma_lhs, arma_rhs);

        total_tests++;
        hybrid_lhs.to_dense();
        if (compare_vs_results(hybrid_lhs.get_dense(), arma_lhs, hybrid_rhs, arma_rhs,
                              "Large voltage value (1e6)")) {
            passed_tests++;
        }
        std::cout << "\n";
    }

    // Test Case 8: Adjacent nodes
    {
        std::cout << "Test Case 8: Adjacent nodes (node_x=1, node_y=2, V=2.5)\n";
        const size_t size = 3;
        HybridMatrix hybrid_lhs(size, size, false);
        arma::mat arma_lhs = arma::zeros(size, size);
        arma::vec hybrid_rhs = arma::zeros(size);
        arma::vec arma_rhs = arma::zeros(size);

        Vs_assigner(1, 2, 2.5, hybrid_lhs, hybrid_rhs);
        Vs_assigner(1, 2, 2.5, arma_lhs, arma_rhs);

        total_tests++;
        hybrid_lhs.to_dense();
        if (compare_vs_results(hybrid_lhs.get_dense(), arma_lhs, hybrid_rhs, arma_rhs,
                              "Adjacent nodes V-source")) {
            passed_tests++;
        }
        std::cout << "\n";
    }

    // Summary
    std::cout << "=========================================\n";
    std::cout << "Vs_assigner RESULTS: " << passed_tests << "/" << total_tests << " tests passed\n";
    if (passed_tests == total_tests) {
        std::cout << " ALL Vs_assigner TESTS PASSED\n";
    } else {
        std::cout << "L SOME Vs_assigner TESTS FAILED\n";
    }
    std::cout << "=========================================\n\n";
}

int main(){
    // Run R_assigner tests
    test_R_assigner();

    // Run Vs_assigner tests
    test_Vs_assigner();

    std::cout << "\n";
    std::cout << "=================================================\n";
    std::cout << "================== FINAL SUMMARY ================\n";
    std::cout << "=================================================\n";
    std::cout << "All function equivalence tests completed!\n";
    std::cout << "Review results above for details.\n";
    std::cout << "=================================================\n";

    return 0;
}

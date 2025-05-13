#pragma once
#include <cstdint>
#include <stdexcept>
#include <vector>

namespace Math {

constexpr std::uint64_t factorial(unsigned n) {
    if (n > 20) throw std::invalid_argument("Input exceeds 20");
    return (n <= 1) ? 1 : n * factorial(n - 1);
}

// First-order divided difference
inline double DividedDiff(const double q0, const double q1, const double h){
    // DD1
    return (q1 - q0) / h;
}

// Higher-order divided differences calculation using recursion
// Calculates the nth order divided difference
// qValues should be provided in chronological order from oldest to newest
// hValues should be the corresponding time step sizes
inline double DividedDiffN(const std::vector<double>& qValues, const std::vector<double>& hValues, unsigned int order) {
    if (order == 0) return qValues.back();
    
    if (order == 1) {
        if (qValues.size() < 2) throw std::invalid_argument("Not enough q values for DD1");
        if (hValues.size() < 1) throw std::invalid_argument("Not enough h values for DD1");
        return DividedDiff(qValues[qValues.size()-2], qValues.back(), hValues.back());
    }
    
    if (qValues.size() < order + 1) throw std::invalid_argument("Not enough q values for the requested order");
    if (hValues.size() < order) throw std::invalid_argument("Not enough h values for the requested order");
    
    // Calculate sum of h values for denominator
    double h_sum = 0.0;
    for (size_t i = hValues.size() - order; i < hValues.size(); ++i) {
        h_sum += hValues[i];
    }
    
    // Calculate lower order divided differences
    std::vector<double> left_q(qValues.begin(), qValues.end() - 1);
    std::vector<double> right_q(qValues.begin() + 1, qValues.end());
    
    std::vector<double> left_h(hValues.begin(), hValues.end() - 1);
    std::vector<double> right_h(hValues.begin() + 1, hValues.end());
    
    // Recursively calculate the divided differences
    double dd_left = DividedDiffN(left_q, left_h, order - 1);
    double dd_right = DividedDiffN(right_q, right_h, order - 1);
    
    return (dd_right - dd_left) / h_sum;
}

// Optimized version for computation of divided differences of any order
// This version uses a non-recursive algorithm with lower memory overhead
// qValues should be provided in chronological order from oldest to newest
// hValues should be corresponding time steps
inline double DividedDiffNOptimized(const std::vector<double>& qValues, const std::vector<double>& hValues, unsigned int order) {
    const size_t n = qValues.size();
    
    if (order == 0) return qValues.back();
    if (n <= order) throw std::invalid_argument("Not enough q values for the requested order");
    if (hValues.size() < n - 1) throw std::invalid_argument("Not enough h values for time steps");
    
    // Special case for first-order
    if (order == 1) return (qValues[n-1] - qValues[n-2]) / hValues[n-2];
    
    // For higher orders, use the table-based approach
    std::vector<double> dd(n);
    
    // Initialize zeroth order divided differences (just the values)
    for (size_t i = 0; i < n; ++i) {
        dd[i] = qValues[i];
    }
    
    // Compute each order of divided difference
    for (size_t k = 1; k <= order; ++k) {
        for (size_t i = 0; i < n - k; ++i) {
            // Calculate time span for this difference
            double dt = 0.0;
            for (size_t j = i; j < i + k; ++j) {
                dt += hValues[j];
            }
            // Update divided difference
            dd[i] = (dd[i+1] - dd[i]) / dt;
        }
    }
    
    return dd[0];
}

} // namespace Math
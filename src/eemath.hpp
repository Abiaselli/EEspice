#pragma once
#include <cstdint>
#include <stdexcept>
#include <vector>

namespace Math {

constexpr std::uint64_t factorial(unsigned n) {
    if (n > 20) throw std::invalid_argument("Input exceeds 20");
    return (n <= 1) ? 1 : n * factorial(n - 1);
}

/// \brief Compute the n-th order divided difference y[x0,...,x_n].
///
/// qValues must contain at least order+1 entries:
///   { y(x0), y(x1), ..., y(x_order) }
/// hValues must contain at least order entries:
///   { x1-x0, x2-x1, ..., x_order - x_{order-1} }
///
/// \param qValues   vector of function values y(x_i)
/// \param hValues   vector of successive spacings x_{i+1}-x_i
/// \param order     the order n of the divided difference
/// \return          the divided difference y[x0,x1,...,x_order]
double DividedDiff(const std::vector<double>& qValues,
                   const std::vector<double>& hValues,
                   unsigned int order)
{
    if (qValues.size() < order + 1)
        throw std::invalid_argument(
            "DividedDiff: need at least order+1 entries in qValues");
    if (hValues.size() < order)
        throw std::invalid_argument(
            "DividedDiff: need at least order entries in hValues");

    // Working copy of the first (order+1) y-values
    std::vector<double> d(qValues.begin(),
                          qValues.begin() + (order + 1));

    // Build prefix-sums so that
    //   prefix[k] == x_k - x_0 = sum_{i=0..k-1} hValues[i]
    std::vector<double> prefix(order + 1, 0.0);
    for (unsigned int i = 0; i < order; ++i) {
        prefix[i + 1] = prefix[i] + hValues[i];
    }

    // In-place construction of the divided-difference table:
    //   for k = 1..order:
    //     for i = 0..order-k:
    //       d[i] = (d[i+1] - d[i]) / (x_{i+k} - x_i)
    //            = (d[i+1] - d[i]) / (prefix[i+k] - prefix[i])
    for (unsigned int k = 1; k <= order; ++k) {
        for (unsigned int i = 0; i + k <= order; ++i) {
            double numerator   = d[i + 1] - d[i];
            double denominator = prefix[i + k] - prefix[i];
            d[i] = numerator / denominator;
        }
    }

    // After order passes, d[0] == y[x0, x1, ..., x_order]
    return d[0];
}

} // namespace Math
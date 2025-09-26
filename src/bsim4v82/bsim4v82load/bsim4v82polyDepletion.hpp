#pragma once
#include <cmath>
#include "sim_variables.hpp"
#include "bsim4v82/bsim4v82const.hpp"

namespace bsim4{
/* function to compute poly depletion effect */
int BSIM4polyDepletion(
    double  phi,
    double  ngate,
    double  epsgate,
    double  coxe,
    double  Vgs,
    double *Vgs_eff,
    double *dVgs_eff_dVg)
{
    double T1, T2, T3, T4, T5, T6, T7, T8;

    /* Poly Gate Si Depletion Effect */
    if ((ngate > 1.0e18) &&
        (ngate < 1.0e25) && (Vgs > phi) && (epsgate!=0)
       ){
        T1 = 1.0e6 * CHARGE * epsgate * ngate / (coxe * coxe);
        T8 = Vgs - phi;
        T4 = std::sqrt(1.0 + 2.0 * T8 / T1);
        T2 = 2.0 * T8 / (T4 + 1.0);
        T3 = 0.5 * T2 * T2 / T1; /* T3 = Vpoly */
        T7 = 1.12 - T3 - 0.05;
        T6 = std::sqrt(T7 * T7 + 0.224);
        T5 = 1.12 - 0.5 * (T7 + T6);
        *Vgs_eff = Vgs - T5;
        *dVgs_eff_dVg = 1.0 - (0.5 - 0.5 / T4) * (1.0 + T7 / T6);
    }
    else {
        *Vgs_eff = Vgs;
        *dVgs_eff_dVg = 1.0;
    }
    return(0);
}

inline void dexp(const double A, double& B, double& C) noexcept {
    if (A > EXP_THRESHOLD) {
        B = MAX_EXP * (1.0 + (A) - EXP_THRESHOLD);
        C = MAX_EXP;
    } else if (A < -EXP_THRESHOLD) {
        B = MIN_EXP;
        C = 0;
    } else {
        B = std::exp(A);
        C = B;
    }
}

} // namespace bsim4
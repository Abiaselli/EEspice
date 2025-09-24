#pragma once
#include <chrono>
#include <vector>

// Global variables
/* Upper transient iteration limit for iteration count in time-step control algorithm*/
/* some convergence issues that get resolved by increasing max iter */
constexpr int ITL4 = 100;
bool converged = false;

/* Scale factor for predicting*/
constexpr double TRTOL = 7.0;

/* Error Charge Factor*/
constexpr double CHGTOL = 1.0e-14;

/* Relative error*/
constexpr double RELTOL = 1.0e-3;

/* Absolute error for voltage in SPICE book*/
constexpr double VNTOL = 1.0e-6;
/* Absolute error for current in SPICE book*/
constexpr double ABSTOL = 1.0e-12;

// Relative charge tolerance for truncation error
constexpr double relq = 0.01;
// Relative current tolerance for truncation error in SPICE OPUS
constexpr double lteretol = 0.01;
// Absolute tolerance for truncation error in SPICE OPUS
constexpr double lteabstol = 1.0e-6;

constexpr double h_op = 1.0e-25; // Initial timestep in OP

constexpr double LargeEpsilon = 1.0e-6;
constexpr double SmallEpsilon = 1.0e-12;


// To control debug mode
bool debugMode = false;

// Batch mode
bool batchMode = false;

// Statistics of the transient simulation
int NR_ITE = 0;
int total_timepoint = 0;

// Constants in eespice
constexpr double CONSTsqrt2 = 1.4142135623730950488016887242097;
constexpr double CONSTpi = 3.1415926535897932384626433832795;       //pi
constexpr double CONSTnap = 2.7182818284590452353602874713527;
constexpr double CONSTlog10e = 0.43429448190325182765112891891661;
constexpr double CONSTlog2e = 1.4426950408889634073599246810019;

/* https://physics.nist.gov/cgi-bin/cuu/Value?c
 * value = 299 792 458 m s-1 (exact) */
constexpr double CONSTc = 299792458;

// https://www.nist.gov/pml/weights-and-measures/si-units-temperature
constexpr double CONSTCtoK = 273.15;
constexpr double CONSTKtoC = -273.15;

/* https://physics.nist.gov/cgi-bin/cuu/Value?e
 *                value = 1.602 176 6208 x 10-19 C
 * standard uncertainty = 0.000 000 0098 x 10-19 C */
constexpr double CHARGE = 1.6021766208e-19;

/* https://physics.nist.gov/cgi-bin/cuu/Value?k
 *                value = 1.380 648 52 x 10-23 J K-1
 * standard uncertainty = 0.000 000 79 x 10-23 J K-1 */
constexpr double CONSTboltz = 1.38064852e-23;

/* https://physics.nist.gov/cgi-bin/cuu/Value?h
 *                value = 6.626 070 040 x 10-34 J s
 * standard uncertainty = 0.000 000 081 x 10-34 J s */
constexpr double CONSTplanck = 6.626070040e-34;

constexpr double nomTemp = 300.15; // Nominal temperature in Kelvin

constexpr double CONSTmuZero = (4.0 * CONSTpi * 1E-7); /* MuZero H/m */

/* epsilon zero  e0*u0*c*c=1 */
constexpr double CONSTepsZero = (1.0 / (CONSTmuZero * CONSTc * CONSTc)); /* F/m */

/* This value is not really constant over temperature and frequency, but
 * 3.9 is the most common "all-purpose" value */
constexpr double CONSTepsrSiO2 = 3.9;

constexpr double CONSTepsSiO2 = (CONSTepsrSiO2 * CONSTepsZero);  /* epsilon SiO2 F/m */
constexpr double REFTEMP = (27.0 + CONSTCtoK); /* 27 degrees C in K */


constexpr double CONSTroot2 = CONSTsqrt2;
constexpr double CONSTvt0 = CONSTboltz * REFTEMP / CHARGE;
constexpr double CONSTKoverQ = CONSTboltz / CHARGE;
constexpr double CONSTe = CONSTnap;

// Integration method
constexpr int BACKWARD_EULER = 0;
constexpr int TRAPEZOIDAL = 1;
constexpr int GEAR = 2;






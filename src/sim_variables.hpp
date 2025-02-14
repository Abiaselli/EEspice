#pragma once
#include <chrono>
#include <vector>
#include "BS_thread_pool/BS_thread_pool.hpp"

// Global variables
/* Upper transient iteration limit for iteration count in time-step control algorithm*/
const constexpr int ITL4 = 10;

/* Scale factor for predicting*/
const constexpr double TRTOL = 7.0;

/* Error Charge Factor*/
const constexpr double CHGTOL = 1.0e-14;

/* Relative error*/
const constexpr double RELTOL = 1.0e-3;

/* Absolute error for voltage in SPICE book*/
const constexpr double VNTOL = 1.0e-6;
/* Absolute error for current in SPICE book*/
const constexpr double ABSTOL = 1.0e-12;

// Relative charge tolerance for truncation error
const constexpr double relq = 0.01;
// Relative current tolerance for truncation error in SPICE OPUS
const constexpr double lteretol = 0.01;
// Absolute tolerance for truncation error in SPICE OPUS
const constexpr double lteabstol = 1.0e-6;

const constexpr double h_op = 1.0e-25; // Initial timestep in OP

const constexpr double LargeEpsilon = 1.0e-6;
const constexpr double SmallEpsilon = 1.0e-12;

BS::thread_pool pool(3);

// To control debug mode
bool debugMode = false;

// Statistics of the transient simulation
int NR_ITE = 0;
int total_timepoint = 0;
int total_NR_iteration = 0;







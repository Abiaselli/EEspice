#pragma once
#include <chrono>
#include <vector>

// Global variables
/* Upper transient iteration limit for iteration count in time-step control algorithm*/
const constexpr int ITL4 = 10;

/* Scale factor for predicting*/
const constexpr double TRTOL = 7.0;

/* Error Charge Factor*/
const constexpr double CHGTOL = 1e-14;

/* Relative error*/
const constexpr double RELTOL = 1e-3;

/* Absolute error for voltage in SPICE book*/
const constexpr double VNTOL = 1e-6;
/* Absolute error for current in SPICE book*/
const constexpr double ABSTOL = 1e-12;

// Relative charge tolerance for truncation error
const constexpr double relq = 0.01;
// Relative current tolerance for truncation error in SPICE OPUS
const constexpr double lteretol = 0.01;
// Absolute tolerance for truncation error in SPICE OPUS
const constexpr double lteabstol = 1e-6;

const constexpr double h_op = 1e-25; // Initial timestep in OP

/*If timestep reach the limit*/
// bool TMAX_reach = false;
// bool TMIN_reach = false;
// auto totalNR = std::chrono::duration<double, std::milli>::zero();
// auto totalEV = std::chrono::duration<double, std::milli>::zero();
// auto totalSOLVE = std::chrono::duration<double, std::milli>::zero();


int NR_ITE{};
// bool ITE_mid;
// bool ITE_up;
// bool ITE_down;

double RD = 1.0;
double RG = 1.0;
double RS = 1.0;


// To control debug mode
bool debugMode = false;

// auto timemid = std::chrono::duration<double, std::milli>::zero();
// auto timeup = std::chrono::duration<double, std::milli>::zero();
// auto timedown = std::chrono::duration<double, std::milli>::zero();

// std::chrono::time_point<std::chrono::high_resolution_clock> mid1, mid2, up1, up2, down1, down2;


int total_timepoint = 0;


// std::vector<double> time_solve;
// std::vector<double> time_ev;
// std::vector<double> time_NR;
// std::vector<double> time_co;
// std::vector<double> time_part1;
// std::vector<double> time_part2;
// std::vector<double> time_part3;
// std::vector<double> time_part4;




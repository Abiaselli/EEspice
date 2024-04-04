#pragma once

// Global variables
/* Scale factor for predicting*/
int TRTOL = 7;

/* Error Charge Factor*/
double CHGTOL = 1e-14;

/* Relative error*/
double RELTOL = 1e-3;

/* Absolute error for voltage in SPICE book*/
double VNTOL = 1e-6;
/* Absolute error for current in SPICE book*/
double ABSTOL = 1e-12;

// Relative charge tolerance for truncation error
double relq = 0.01;
// Relative current tolerance for truncation error in SPICE OPUS
double lteretol = 0.01;
// Absolute tolerance for truncation error in SPICE OPUS
double lteabstol = 1e-6;

/*If timestep reach the limit*/
// bool TMAX_reach = false;
// bool TMIN_reach = false;




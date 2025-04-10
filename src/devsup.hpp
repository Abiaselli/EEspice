#pragma once
#include <cmath>

/* 
 * Limit the per-iteration change of VDS 
 */
double
DEVlimvds(double vnew, double vold)
{

    if(vold >= 3.5) {
        if(vnew > vold) {
            vnew = std::min(vnew,(3 * vold) +2);
        } else {
            if (vnew < 3.5) {
                vnew = std::max(vnew,2.0);
            }
        }
    } else {
        if(vnew > vold) {
            vnew = std::min(vnew,4.0);
        } else {
            vnew = std::max(vnew,-.5);
        }
    }
    return(vnew);
}

/* 
 * Limit the per-iteration change of FET voltages 
 *
 * This code has been fixed by Alan Gillespie: a new
 * definition for vtstlo. 
 */
double
DEVfetlim(double vnew, double vold, double vto)
{
    double vtsthi;
    double vtstlo;
    double vtox;
    double delv;
    double vtemp;

    vtsthi = std::abs(2*(vold-vto))+2;
    vtstlo = std::abs(vold-vto)+1;
    vtox = vto + 3.5;
    delv = vnew-vold;

    if (vold >= vto) {
        if(vold >= vtox) {
            if(delv <= 0) {
                /* going off */
                if(vnew >= vtox) {
                    if(-delv >vtstlo) {
                        vnew =  vold - vtstlo;
                    }
                } else {
                    vnew = std::max(vnew,vto+2);
                }
            } else {
                /* staying on */
                if(delv >= vtsthi) {
                    vnew = vold + vtsthi;
                }
            }
        } else {
            /* middle region */
            if(delv <= 0) {
                /* decreasing */
                vnew = std::max(vnew,vto-.5);
            } else {
                /* increasing */
                vnew = std::min(vnew,vto+4);
            }
        }
    } else {
        /* off */
        if(delv <= 0) {
            if(-delv >vtsthi) {
                vnew = vold - vtsthi;
            } 
        } else {
            vtemp = vto + .5;
            if(vnew <= vtemp) {
                if(delv >vtstlo) {
                    vnew = vold + vtstlo;
                }
            } else {
                vnew = vtemp;
            }
        }
    }
    return(vnew);
}

/*  
 * Limit the per-iteration change of PN junction voltages 
 *
 * This code has been fixed by Alan Gillespie adding limiting
 * for negative voltages.
 */
double
DEVpnjlim(double vnew, double vold, double vt, double vcrit, int *icheck)
{
    double arg;

    if((vnew > vcrit) && (std::abs(vnew - vold) > (vt + vt))) {
        if(vold > 0.0) {
            arg = (vnew - vold) / vt;
            if(arg > 0.0) {
                vnew = vold + vt * (2+std::log(arg-2));
            } else {
                vnew = vold - vt * (2+std::log(2-arg));
            }
        } else {
            vnew = vt * std::log(vnew/vt);
        }
        *icheck = 1;
    } else {
       if (vnew < 0.0) {
           if (vold > 0.0) {
               arg = -1*vold-1;
           } else {
               arg = 2.0*vold-1;
           }
           if (vnew < arg) {
              vnew = arg;
              *icheck = 1;
           } else {
              *icheck = 0;
           }
        } else {
           *icheck = 0;
        }
    }
    return(vnew);
}

#ifndef PHASEFIELD_H_
#define PHASEFIELD_H_

#include <iostream>
#include <stdio.h>
#include <math.h>

void Init_RND(int n, double x, double rnd_init[]);

void d_dxx(double P_gp, double gp2[], double* deriv);

#endif

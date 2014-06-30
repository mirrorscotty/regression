#ifndef STRESS_STRAIN_H
#define STRESS_STRAIN_H

#include "pasta.h"
#include "matrix.h"

double strain(double, double);
double dstrain(double, double);
double stress_model(double t, matrix*);
matrix* maxwell_stress(double, double, maxwell*, matrix*, matrix*, double, double);
matrix* fit_stress(double, double, maxwell*, double, double);
double storage_mod(double, double, double);
double loss_mod(double, double, double);

#endif


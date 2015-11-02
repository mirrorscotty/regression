#ifndef STRESS_STRAIN_H
#define STRESS_STRAIN_H

#include "material-data.h"
#include "matrix.h"

double strain(double, double);
double dstrain(double, double);
double stress_model(double t, matrix*);
matrix* maxwell_stress(maxwell*, matrix*, matrix*, double, double);
matrix* fit_stress(double, double, maxwell*, double, double);
double storage_mod(double, double, double);
double loss_mod(double, double, double);

matrix* maxwell_stress_rozzi(matrix*, matrix*, double, double);
matrix* fit_stress_rozzi(double, double, double, double);

#endif


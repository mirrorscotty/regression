#ifndef REGRESS_H
#define REGRESS_H

#include "matrix.h"

matrix* regress(matrix*, matrix*);
matrix* polyfit(matrix*, matrix*, int);
matrix* fitnlm(double (*)(double, matrix*), matrix*, matrix*, matrix*);

#endif

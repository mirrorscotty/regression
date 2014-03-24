#ifndef FITNLM_H
#define FITNLM_H

#include "matrix.h"

matrix* fitnlm(double (*)(double, matrix*), matrix*, matrix*, matrix*);

#endif


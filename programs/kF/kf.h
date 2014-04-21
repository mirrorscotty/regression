#ifndef KF_H
#define KF_H

#include "matrix.h"

#define CONSTX0 0
#define CONSTXe 18.261700
#define CONSTnterms 50
#define BETA0 1e-4

double CrankEquation(double, double, double, double, int);
double CrankkF(double, double, double, double);
double CrankModel(double, matrix*);

matrix* LoadIGASorpTime(char*);
matrix* LoadIGASorpXdb(char*, double);

double CalcXe(matrix*, matrix*, double);

double fitsubset(matrix*, matrix*, int, int);
void calckf(matrix*, matrix*, double, char*);
void calckfstep(matrix*, matrix*, double, char*);
void fitkf(matrix*, matrix*, char*);

#endif


#ifndef KF_H
#define KF_H

#include "matrix.h"
#include "mechanical.h"

#define CONSTX0 0
#define CONSTXe 18.261700
#define CONSTnterms 50 
#define BETA0 1e-4

double CrankEquation(double, double, double, double, int);
double CrankkF(double, double, double, double, double);
double CrankModel(double, matrix*);

matrix* LoadIGASorpTime(char*);
matrix* LoadIGASorpXdb(char*, double);
matrix* LoadIGASorpRH(char*);

double CalcXe(int, matrix*, matrix*, double);
double NCalcXe(int, matrix*, matrix*, double);
double CalcXeIt(int, matrix*, matrix*, double);

double fitsubset(matrix*, matrix*, int, int);
matrix* calckf(matrix*, matrix*, double);
matrix* calckfstep(matrix*, matrix*, double);
matrix* fitkf(matrix*, matrix*);

int FindInitialPointkF(matrix*);
int FindInitialPointRH(matrix*);

double DeborahNumber(int, int, matrix*, matrix*, double, double, maxwell*);
double NewLength(int, int, matrix*, matrix*, double, double);
matrix* DeborahMatrix(int, matrix*, matrix*, double, double, maxwell*);
matrix* LengthMatrix(int, matrix*, matrix*, double, double);


#endif


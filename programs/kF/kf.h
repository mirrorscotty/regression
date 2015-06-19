#ifndef KF_H
#define KF_H

#include "matrix.h"
#include "material-data.h"

#define CONSTX0 0
#define CONSTXe 18.261700
#define CONSTnterms 50
#define BETA0 1e-4

#define SLABWIDTH 6e-3
#define SLABLENGTH 8e-3

double CrankEquation(double, double, double, double, int);
double CrankkF(double, double, double, double, double);
double CrankModel(double, matrix*);

vector* LoadIGASorpTime(char*);
vector* LoadIGASorpXdb(char*, double);
vector* LoadIGASorpRH(char*);

double CalcXe(int, matrix*, matrix*, double);
double NCalcXe(int, vector*, vector*, double);
double CalcXeIt(int, vector*, vector*, double);

double fitsubset(matrix*, matrix*, int, int);
vector* calckf(vector*, vector*, double);
matrix* calckfstep(matrix*, matrix*, double);
matrix* fitkf(matrix*, matrix*);

int FindInitialPointkF(vector*);
int FindInitialPointRH(vector*);

double DeborahNumber(int, int, vector*, vector*, double, double, maxwell*);
double NewLength(int, int, vector*, vector*, double, double);
vector* DeborahMatrix(int, vector*, vector*, double, double, maxwell*);
vector* LengthMatrix(int, vector*, vector*, double, double);
vector* LengthConstD(int, vector*, double, double);
vector* LengthWaterLoss(int, vector*, double, double, double);
vector* LengthDensityChange(int, vector*, double, double, double);
vector* DOswinVector(int, vector*, double);

vector* MassFlux(int, vector*, vector*, double);
vector* MomentumFlux(int, vector*, vector*, vector*, double, maxwell*);
vector* PastaMassFlux(int, vector*, vector*, double, double);

#endif


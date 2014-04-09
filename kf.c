/**
 * @file kf.c
 * Solve the Crank equation for kF using a nonlinear regression scheme.
 */

#include "matrix.h"
#include "fitnlm.h"
#include <math.h>
#include <stdio.h>

double CrankEquation(double kf, double t, double X0, double Xe, int nterms)
{
    double value = 0;
    int n;
    for(n=0; n<nterms; n++) {
        value += 8/(pow(2*n+1, 2) * M_PI*M_PI) *
            exp(-kf * t * pow(2*n+1, 2));
    }

    value = value * (X0-Xe) + Xe;

    return value;
}

double CrankModel(double t, matrix *beta)
{
    double X0 = 0.3168046007,
           Xe = 0.1123, 
           kf = val(beta, 0, 0);
    int nterms = 50;

    return CrankEquation(kf, t, X0, Xe, nterms);
}

double fitsubset(matrix *x, matrix *y, int rowstart, int rowend)
{
    matrix *xx, *yy, *beta, *beta0;
    int nrows = rowend-rowstart,
        i;

    xx = CreateMatrix(nrows, 1);
    yy = CreateMatrix(nrows, 1);
    beta0 = CreateMatrix(1, 1);

    setval(beta0, 1e-4, 0, 0);

    for(i=0; i<nrows; i++) {
        setval(xx, val(x, i+rowstart, 0), i, 0);
        setval(yy, val(y, i+rowstart, 0), i, 0);
    }

    beta = fitnlm(&CrankModel, xx, yy, beta0);
    return val(beta, 0, 0);
}

int main(int argc, char *argv[])
{
    matrix *data, *Xdb, *t, *kf;
    int tcol = 1,
        xdbcol = 4,
        chunksize = 500,
        rowstart = 0,
        nchunks, i;

    data = mtxloadcsv(argv[1], 0);

    Xdb = ExtractColumn(data, xdbcol);
    t = ExtractColumn(data, tcol);

    /* Truncated because of integer division */
    nchunks = (nRows(t)-rowstart)/chunksize;
    kf = CreateMatrix(nchunks, 1);

    for(i=0; i<nchunks; i++)
        setval(kf, fitsubset(t, Xdb, rowstart+i*chunksize, rowstart+(i+1)*chunksize), i, 0);

    mtxprntfile(kf, "kF.csv");
    return 0;
}


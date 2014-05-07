/**
 * @file flux.c
 * Calculate the mass and momentum flux for the sample at each data point
 */

#include "kf.h"
#include "matrix.h"
#include "pasta.h"

/* TODO: Calculate average flux over several data points. */
matrix* MassFlux(int initial, matrix *t, matrix *Xdb, double Mdry)
{
    matrix *J; /* Moisture flux matrix */
    double area = SLABLENGTH*SLABWIDTH,
           DMwatDt,
           Ji;
    int i;

    J = CreateMatrix(nRows(Xdb), 1);

    for(i=initial; i<nRows(J); i++) {
        DMwatDt = (val(Xdb, i, 0) - val(Xdb, i-1, 0))
            / (val(t, i, 0) - val(t, i-1, 0)) * Mdry;
        Ji = 0.5*DMwatDt/area;
        setval(J, Ji, i, 0);
    }

    return J;
}

/* TODO: Double check this function to make sure it's giving good results */
matrix* MomentumFlux(int initial, matrix *t,
                     matrix *Xdb, matrix *L,
                     double T, maxwell *m)
{
    matrix *M;
    double Mi = 0, /* Set momentum flux to zero initially */
           L0,
           ti;
    int i, j;

    M = CreateMatrix(nRows(t), 1);
    L0 = val(L, initial, 0);

    for(i=initial; i<nRows(t); i++) {
        for(j=initial; j<i; j++) {
            ti = val(t, j, 0);
            Mi += MaxwellModulus(m, ti, T, val(Xdb, j, 0))
                * (val(L, j, 0) - val(L, j-1, 0))/(val(t, j, 0)-val(t, j-1, 0))
                / L0 * (val(t, j, 0)-val(t, j-1, 0));
        }
        setval(M, Mi, i, 0);
        Mi = 0;
    }

    return M;
}


#include <math.h>
#include "pasta.h"
#include "matrix.h"
#include "regress.h"

double w; /* Frequency global variable */

double strain(double e0, /* Strain magnitude */
              double t) /* time */
{
    return e0*sin(t*w);
}

double dstrain(double e0, double t)
{
    return e0*w*cos(t*w);
}

double stress_model(double t, matrix *beta)
{
    double s0 = val(beta, 0, 0),
           shift = val(beta, 1, 0);

    return s0*sin(t*w+shift);
}

matrix* maxwell_stress(double e0, double freq, maxwell *m, matrix *t, matrix *de, double T, double M)
{
    matrix *s;
    int i, j;
    double dt = val(t, 1, 0) - val(t, 0, 0),
           stress;

    s = CreateMatrix(nRows(t), 1);

    for(i=0; i<nRows(t); i++) {
        for(j=0; j<i; j++) {
            stress = 0;
            stress += MaxwellModulus(m, val(t, j, 0), T, M)
                * val(de, j, 0) * dt;
        }
        setval(s, stress, i, 0);
    }

    return s;
}

matrix* fit_stress(double e0, double freq, maxwell *m, double T, double X)
{
    int i, 
        npts = 1000;
    double dt = .1,
           s0guess = e0,
           shiftguess = 0;
    matrix *t, *e, *beta0, *beta, *de, *s;

    w = freq;
    
    t = CreateMatrix(npts, 1);
    e = CreateMatrix(npts, 1);
    de = CreateMatrix(npts, 1);
    for(i=0; i<npts; i++) {
        setval(t, i*dt, i, 0);
        setval(e, strain(e0, i*dt), i, 0);
        setval(de, dstrain(e0, i*dt), i, 0);
    }

    beta0 = CreateMatrix(2, 1);
    setval(beta0, s0guess, 0, 0);
    setval(beta0, shiftguess, 1, 0);

    s = maxwell_stress(e0, freq, m, t, de, T, X);

    beta = fitnlm(&stress_model, t, s, beta0);

    return beta;
}

double storage_mod(double e0, double s0, double shift)
{
    return s0/e0 * cos(shift);
}

double loss_mod(double e0, double s0, double shift)
{
    return s0/e0 * sin(shift);
}


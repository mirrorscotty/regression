/**
 * @file flux.c
 * Calculate the mass and momentum flux for the sample at each data point
 */

#include "kf.h"
#include "matrix.h"
#include "pasta.h"
#include "choi-okos.h"

#define NPTS 50 /* Number of points to average the flux over */

/**
 * Calculate the mass flux of water leaving the surface of the pasta slab.
 * @param initial Row number of the first data point to consider
 * @param t Column matrix of time values [s]
 * @param Xdb Column matrix of moisture contents [kg/kg db]
 * @param Mdry Mass of the bone-dry sample [kg]
 * @returns Column matrix of mass flux values [kg/(m^2 s)]
 */
matrix* MassFlux(int initial, matrix *t, matrix *Xdb, double Mdry)
{
    matrix *J; /* Moisture flux matrix */
    double area = SLABLENGTH*SLABWIDTH, /* Sample area */
           DMwatDt, /* Derivative of water mass w.r.t time */
           Ji; /* Value of J at data point i */
    int i, j; /* Loop indicies */

    J = CreateMatrix(nRows(Xdb), 1);

    /* Only calculate flux every few data points */
    for(i=initial; i<nRows(J); i+=NPTS) {
        /* Calculate range of change of mass of water in the pasta. */
        DMwatDt = (val(Xdb, i, 0) - val(Xdb, i-NPTS, 0))
            / (val(t, i, 0) - val(t, i-NPTS, 0)) * Mdry;
        /* Calculate total flux over the time period */
        Ji = 0.5*DMwatDt/area;
        /* Calculate the average flux at each data point and save it to the
         * output matrix */
        for(j=i-NPTS; j<i; j++) 
            setval(J, Ji/NPTS, j, 0);
    }

    return J;
}

/* TODO: Double check this function to make sure it's giving good results */
matrix* MomentumFlux(int initial, matrix *t,
                     matrix *Xdb, matrix *L,
                     double T, maxwell *m)
{
    matrix *M, *tadj;
    double Mi = 0, /* Set momentum flux to zero initially */
           L0,
           ti;
    int i, j;

    M = CreateMatrix(nRows(t), 1);
    L0 = val(L, initial, 0);

    tadj = CreateMatrix(nRows(t), 1);
    for(i=0; i<nRows(tadj); i++)
        setval(tadj, val(t, i, 0)-val(t, initial, 0), i, 0);

    for(i=initial; i<nRows(t); i+=1) {
        ti = val(tadj, i, 0);
        for(j=initial; j<i; j++) {
            Mi += MaxwellModulus(m, ti-val(tadj, j, 0), T, val(Xdb, j, 0))
                * (val(L, j, 0) - val(L, j-1, 0)) / L0;
            //printf("G(t) = %g\t", MaxwellModulus(m, ti, T, val(Xdb, j, 0)));
        }

        for(j=i-1; j<i; j++)
            setval(M, Mi/1, j, 0);
        Mi = 0;
    }

    return M;
}

/**
 * Calculate the mass flux of pasta at the surface of the slab.
 * @param initial Row number of the first data point to consider
 * @param t Column matrix of time values [s]
 * @param L Column matrix of slab thicknesses [m]
 * @param L0 Initial thickness [m]
 * @param T Drying temperature [K]
 * @returns Column matrix of mass flux values [kg/(m^2 s)]
 */
matrix* PastaMassFlux(int initial, matrix *t, matrix *L, double L0, double T)
{
    matrix *M; /* Moisture flux matrix */
    double DLDt, /* Derivative of thickness with respect to time */
           v, /* Velocity */
           rhop; /* Pasta density (excluding water content) */
    int i, j; /* Loop indicies */
    choi_okos *co;

    /* Calculate pasta density */
    co = CreateChoiOkos(PASTACOMP);
    rhop = rho(co, T);
    DestroyChoiOkos(co);

    M = CreateMatrix(nRows(t), 1);

    /* Calculate the flux every few points only. then average the values out to
     * get a flux value for each time. */
    for(i=initial; i<nRows(M); i+=NPTS) {
        /* Calculate change in slab thickness with respect to time */
        DLDt = (val(L, i, 0) - val(L, i-NPTS, 0))
            / (val(t, i, 0) - val(t, i-NPTS, 0));
        /* Calculate velocity, assuming symmetry */
        v = 0.5*DLDt;
        /* Average out the value over NPTS and set each value equal to the
         * average */
        for(j=i-NPTS; j<i; j++) 
            setval(M, v*rhop/NPTS, j, 0);
    }

    return M;
}


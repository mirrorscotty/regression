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
vector* MassFlux(int initial, vector *t, vector *Xdb, double Mdry)
{
    vector *J; /* Moisture flux matrix */
    double area = SLABLENGTH*SLABWIDTH, /* Sample area */
           DMwatDt, /* Derivative of water mass w.r.t time */
           Ji; /* Value of J at data point i */
    int i, j; /* Loop indicies */

    J = CreateVector(len(Xdb));

    /* Only calculate flux every few data points */
    for(i=initial; i<len(J); i+=NPTS) {
        /* Calculate range of change of mass of water in the pasta. */
        DMwatDt = (valV(Xdb, i) - valV(Xdb, i-NPTS))
            / (valV(t, i) - valV(t, i-NPTS)) * Mdry;
        /* Calculate total flux over the time period */
        Ji = 0.5*DMwatDt/area;
        /* Calculate the average flux at each data point and save it to the
         * output matrix */
        for(j=i-NPTS; j<i; j++)
            setvalV(J, j, Ji/NPTS);
    }

    return J;
}

/* TODO: Double check this function to make sure it's giving good results */
vector* MomentumFlux(int initial, vector *t,
                     vector *Xdb, vector *L,
                     double T, maxwell *m)
{
    vector *M, *tadj;
    double Mi = 0, /* Set momentum flux to zero initially */
           L0,
           ti;
    int i, j;

    M = CreateVector(len(t));
    L0 = valV(L, initial);

    tadj = CreateVector(len(t));
    for(i=0; i<len(tadj); i++)
        setvalV(tadj, i, valV(t, i)-valV(t, initial));

    for(i=initial; i<len(t); i+=1) {
        ti = valV(tadj, i);
        for(j=initial; j<i; j++) {
            Mi += MaxwellModulus(m, ti-valV(tadj, j), T, valV(Xdb, j))
                * (valV(L, j) - valV(L, j-1)) / L0;
            //printf("G(t) = %g\t", MaxwellModulus(m, ti, T, val(Xdb, j, 0)));
        }

        for(j=i-1; j<i; j++)
            setvalV(M, j, Mi/1);
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
vector* PastaMassFlux(int initial, vector *t, vector *L, double L0, double T)
{
    vector *M; /* Moisture flux vector */
    double DLDt, /* Derivative of thickness with respect to time */
           v, /* Velocity */
           rhop; /* Pasta density (excluding water content) */
    int i, j; /* Loop indicies */
    choi_okos *co;

    /* Calculate pasta density */
    co = CreateChoiOkos(PASTACOMP);
    rhop = rho(co, T);
    DestroyChoiOkos(co);

    M = CreateVector(len(t));

    /* Calculate the flux every few points only. then average the values out to
     * get a flux value for each time. */
    for(i=initial; i<len(M); i+=NPTS) {
        /* Calculate change in slab thickness with respect to time */
        DLDt = (valV(L, i) - valV(L, i-NPTS))
            / (valV(t, i) - valV(t, i-NPTS));
        /* Calculate velocity, assuming symmetry */
        v = 0.5*DLDt;
        /* Average out the value over NPTS and set each value equal to the
         * average */
        for(j=i-NPTS; j<i; j++)
            setvalV(M, j, v*rhop/NPTS);
    }

    return M;
}


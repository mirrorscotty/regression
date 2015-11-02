/**
 * @file stress-strain.c
 * Set of functions to determine the storage and loss moduli of viscoelastic
 * materials.
 */

#include <math.h>
#include "material-data.h"
#include "matrix.h"
#include "regress.h"
#include "stress-strain.h"

/** Frequency global variable */
double w;

/**
 * Calculate the stress on a viscoelastic material using the Maxwell model
 * relaxation function with temperature and moisture effects.
 *
 * \f[
 * \sigma(t) = \int_{\tau_0}^t\! G(t-\tau)\dot{\epsilon}(\tau)\,\mathrm{d}\tau
 * \f]
 *
 * @param m Maxwell material parameters
 * @param t Column matrix of time values [s]
 * @param de Matrix containing the time derivative of strain [1/s]
 * @param T Temperature [K]
 * @param M Moisture content [kg/kg db]
 * @returns Column matrix of stress values [-]
 */
matrix* maxwell_stress_rozzi(matrix *t, matrix *de,
                       double T, double M)
{
    matrix *s;
    int i, j; /* Loop indicies */
    double dt = val(t, 1, 0) - val(t, 0, 0), /* Delta t */
           stress; /* Intermediate value used for integration */

    /* Make a matrix to store the output in */
    s = CreateMatrix(nRows(t), 1);

    /* Numerically integrate the equation for each time step */
    for(i=0; i<nRows(t); i++) {
        for(j=0; j<i; j++) {
            stress = 0;
            stress += MaxwellRelaxLaura(val(t, j, 0), T, M)
                * val(de, j, 0) * dt;
        }
        setval(s, stress, i, 0);
    }

    return s;
}

/**
 * Fit the measured stress to calculate stress magnitude and phase lag.
 * Stress-strain data is generated based on the supplied strain magnitude,
 * oscillation frequency, Maxwell parameters, temperature, and moisture content.
 * The stress magnitude and phase lag are then calculated using nonlinear
 * regression.
 * @param e0 Strain magnitude [-]
 * @param freq Oscillation frequency [1/s]
 * @param m Maxwell parameters
 * @param T Temperature [K]
 * @param X Moisture content [kg/kg db]
 * @returns A 2x1 matrix. Element 1,1 is strain magnitude, and element 2,1 is
 *      phase lag.
 */
matrix* fit_stress_rozzi(double e0, double freq, double T, double X)
{
    int i, /* Loop index */
        npts = 1000; /* Number of points to use for fitting the data */
    double dt = .1, /* Time step size to use when generating data */
           s0guess = e0*MaxwellRelaxLaura(.01, T, X), /* Initial guess for the stress magnitude */
           shiftguess = .3; /* Initial guess for phase lag */
    matrix *t, /* Time matrix */
           *e, /* Strain */
           *beta0, /* Initial guess for fitting coefficients */
           *beta, /* Fitting parameter matrix */
           *de, /* Time deriviative of strain */
           *s; /* Stress */

    /* Set the global variable to the supplied frequency */
    w = freq;

    /* Create matricies for time, strain, and strain rate */
    t = CreateMatrix(npts, 1);
    e = CreateMatrix(npts, 1);
    de = CreateMatrix(npts, 1);

    /* Calculate values for time, strain, and strain rate */
    for(i=0; i<npts; i++) {
        setval(t, i*dt, i, 0);
        setval(e, strain(e0, i*dt), i, 0);
        setval(de, dstrain(e0, i*dt), i, 0);
    }

    /* Create the initial guess matrix for the regression parameters */
    beta0 = CreateMatrix(2, 1);
    setval(beta0, s0guess, 0, 0);
    setval(beta0, shiftguess, 1, 0);

    /* Calculate the values for stress at each point in time based on the
     * Maxwell material model */
    s = maxwell_stress_rozzi(t, de, T, X);

    /* Fit the stress to the appropriate equation to find stress magnitude and
     * phase lag */
    beta = fitnlm(&stress_model, t, s, beta0);

    /* Return the results */
    return beta;
}


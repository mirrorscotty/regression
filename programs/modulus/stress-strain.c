/**
 * @file stress-strain.c
 * Set of functions to determine the storage and loss moduli of viscoelastic
 * materials.
 */

#include <math.h>
#include "pasta.h"
#include "matrix.h"
#include "regress.h"

/** Frequency global variable */
double w;

/**
 * Calculate the imposed strain based on the strain magnitude, oscillation
 * frequency, and current time. Frequency is supplied via a global variable.
 * @param e0 Strain magnitude [-]
 * @param t Time [s]
 * @returns Strain [-]
 */
double strain(double e0, double t)
{
    return e0*sin(t*w);
}

/**
 * Time derivative of imposed strain.
 * @param e0 Strain magnitude [-]
 * @param t Time [s]
 * @returns Time derivative of strain [1/s]
 */
double dstrain(double e0, double t)
{
    return e0*w*cos(t*w);
}

/**
 * Function used to fit the calculate stress to the function:
 * \f[ \sigma = \sigma_0 \sin(t w + \delta) \f]
 * @param t Time [s]
 * @param beta Coefficient matrix. The first one is strain magnitude, and the
 *      second is phase lag.
 * @returns Stress [-]
 */
double stress_model(double t, matrix *beta)
{
    double s0 = val(beta, 0, 0),
           shift = val(beta, 1, 0);

    return s0*sin(t*w+shift);
}

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
matrix* maxwell_stress(maxwell *m, matrix *t, matrix *de,
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
            stress += MaxwellModulus(m, val(t, j, 0), T, M)
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
matrix* fit_stress(double e0, double freq, maxwell *m, double T, double X)
{
    int i, /* Loop index */
        npts = 1000; /* Number of points to use for fitting the data */
    double dt = .1, /* Time step size to use when generating data */
           s0guess = e0, /* Initial guess for the stress magnitude */
           shiftguess = 0; /* Initial guess for phase lag */
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
    s = maxwell_stress(m, t, de, T, X);

    /* Fit the stress to the appropriate equation to find stress magnitude and
     * phase lag */
    beta = fitnlm(&stress_model, t, s, beta0);

    /* Return the results */
    return beta;
}

/**
 * Calculate the storage modulus of a viscoelastic material given the strain
 * magnitude, stress magnitude, and phase lag between the two.
 * \f[ E' = \frac{\sigma_0}{\epsilon_0} \cos\delta \f]
 * @param e0 Strain magnitude
 * @param s0 Stress magnitude
 * @param shift Phase lag
 * @returns Storage modulus
 */
double storage_mod(double e0, double s0, double shift)
{
    return s0/e0 * cos(shift);
}

/**
 * Calculate the loss modulus of a viscoelastic material given the strain
 * magnitude, stress magnitude, and phase lag between the two.
 * \f[ E'' = \frac{\sigma_0}{\epsilon_0} \sin\delta \f]
 * @param e0 Strain magnitude
 * @param s0 Stress magnitude
 * @param shift Phase lag
 * @returns Loss modulus
 */
double loss_mod(double e0, double s0, double shift)
{
    return s0/e0 * sin(shift);
}


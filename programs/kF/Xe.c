#include "kf.h"
#include "matrix.h"
#include "regress.h"
#include <stdio.h>
#include <math.h>

/**
 * Calculate the equilibrium moisture content. This function determines the best
 * value of Xe to make a plot of \f$\ln\frac{X-X_e}{X_0-X_e}\f$ vs time linear.
 * In order to do this, it fits the data to the equation \f$y = a t + b\f$, where
 * \f$y = \ln(X-X_e)\f$ and \f$b = \ln(X_0-X_e)\f$, and then solves the equation
 * \f$F(X_e) = b - \ln(X_0-X_e) = 0\f$ using Newton's method.
 *
 * @param t Column matrix containing time during drying [s]
 * @param Xdb Column matrix of moisture content [kg/kg db]
 * @param Xe0 Initial guess for equilibrium moisture content.
 * @returns Equilibrium moisture content [kg/kg db]
 */
double CalcXe(int initial, matrix *t, matrix *Xdb, double Xe0)
{
    double f, df, /* Function values and derivatives */
           b, /* Fitting parameters. Only the constant matters */
           tol = 1e-10, /* Tolerance for Newton's method */
           Xe = Xe0, /* Set Xe to the initial guess */
           Xep, /* Previous guess */
           X0; /* Initial moisture content */
    matrix *beta, /* Matrix of fitting values */
           *y, /* Set equal to ln(X - Xe) */
           *Xadj, 
           *tadj;
    int i, /* Loop index */
        iter = 0; /* Current iteration */

    /* Set the initial moisture content */
    X0 = val(Xdb, initial, 0);

    /* Make smaller matricies that contain only the "good" data. */
    tadj = CreateMatrix(nRows(Xdb) - initial, 1);
    Xadj = CreateMatrix(nRows(Xdb) - initial, 1);

    for(i=initial; i<nRows(t); i++) {
        setval(tadj, val(t, i, 0), i-initial, 0);
        setval(Xadj, val(Xdb, i, 0), i-initial, 0);
    }

    /* Actually find Xe */
    do {
        /* Make a y matrix containing ln(Xdb - Xe) */
        y = CreateMatrix(nRows(Xadj), 1);
        for(i=0; i<nRows(Xadj); i++)
            setval(y, log(val(Xadj, i, 0) - Xe), i, 0);

        /* Calculate b */
        beta = polyfit(tadj, y, 1);
        b = val(beta, 0, 0);

        /* Calculate f and df */
        f = b - log(X0 - Xe);
        df = 1/(X0 - Xe);

        /* Calculate the new value of Xe */
        Xep = Xe;
        Xe = Xe - f/df;

        /* Clean up */
        DestroyMatrix(y);
        DestroyMatrix(beta);

        /* Keep track of how many iterations we've gone through */
        iter++;

        /* Print out the current value */
        printf("Xe = %g\r", Xe);
    } while( fabs(Xe - Xep) > tol ); /* Check our value */

    /* Print out how many iterations it took to find Xe */
    printf("Solution converged after %d iterations.\n", iter);

    return Xe;
}


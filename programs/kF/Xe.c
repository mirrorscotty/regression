#include "kf.h"
#include "matrix.h"
#include "regress.h"
#include <stdio.h>
#include <math.h>

/**
 * Global variable used to keep track of initial moisture content between
 * NCalcXe and XeModel.
 *
 * @see NCalcXe XeModel
 */
double Xinit = 0;

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
 *
 * @see polyfit
 */
double CalcXe(int initial, matrix *t, matrix *Xdb, double Xe0)
{
    double f, df, /* Function values and derivatives */
           b, /* Fitting parameters. Only the constant matters */
           tol = 1e-10, /* Tolerance for Newton's method */
           Xe = Xe0, /* Set Xe to the initial guess */
           Xep, /* Previous guess */
           X0, /* Initial moisture content */
           r2;
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
        r2 = rsquared(tadj, y, beta);
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
        printf("Xe = %g, R^2 = %g\n", Xe, r2);
        if(Xe < 0) {
            printf("Failure to converge after %d iterations.\n", iter);
            return 0;
        }
    } while( fabs(Xe - Xep) > tol ); /* Check our value */

    /* Print out how many iterations it took to find Xe */
    printf("Solution converged after %d iterations.\n", iter);

    return Xe;
}

/**
 * Model for fitting \f[\ln\frac{X-Xe}{X0-Xe} = at \f] in order to find Xe. This
 * is used by the fitnlm function to calculate the fitting parameters.
 * @param t Time [s]
 * @param beta Set of fitting parameters. Row 0, col 0 is Xe and row 1, col 0 is
 *      the kf parameter.
 * @returns Moisture content [kg/kg db]
 *
 * @see fitnlm
 */
double XeModel(double t, matrix *beta)
{
    double Xe = val(beta, 0, 0),
           kf = val(beta, 1, 0);

    return (Xinit - Xe) * exp(kf*t) + Xe;
}

/**
 * Calculate equilibrium moisture content (Xe) using nonlinear regression.
 * @param initial Row number of the first data point to use in the calculation
 * @param t Column matrix of time values [s]
 * @param Xdb Column matrix of moisture contents [kg/kg db]
 * @param Xe0 Initial guess for Xe. The initial value for kf is hard coded below
 * @returns Equilibrium moisture content [kg/kg db]
 *
 * @see fitnlm XeModel
 */
double NCalcXe(int initial, matrix *t, matrix *Xdb, double Xe0)
{
    double f, df, /* Function values and derivatives */
           b, /* Fitting parameters. Only the constant matters */
           tol = 1e-10, /* Tolerance for Newton's method */
           Xe = Xe0, /* Set Xe to the initial guess */
           Xep; /* Previous guess */
    matrix *beta, /* Matrix of fitting values */
           *beta0,
           *Xadj, 
           *tadj;
    int i, /* Loop index */
        iter = 0; /* Current iteration */

    /* Set the initial moisture content */
    Xinit = val(Xdb, initial, 0);
    beta0 = CreateMatrix(2,1);
    setval(beta0, Xe0, 0, 0);
    setval(beta0, .007, 1, 0);

    /* Make smaller matricies that contain only the "good" data. */
    tadj = CreateMatrix(nRows(Xdb) - initial, 1);
    Xadj = CreateMatrix(nRows(Xdb) - initial, 1);

    for(i=initial; i<nRows(t); i++) {
        setval(tadj, val(t, i, 0), i-initial, 0);
        setval(Xadj, val(Xdb, i, 0), i-initial, 0);
    }

    /* Actually find Xe */
    beta = fitnlm(&XeModel, tadj, Xadj, beta0);

    Xe = val(beta, 0, 0);
    mtxprnt(beta);
    mtxprnt(beta0);

    DestroyMatrix(beta);
    DestroyMatrix(beta0);
    DestroyMatrix(tadj);
    DestroyMatrix(Xadj);

    return Xe;
}

/**
 * Calculate the equilibrium moisture content. This algorithm takes an initial
 * guess for Xe and then fits a linear equation to
 * \f$\ln\frac{X-X_e}{X_0-X_e}\f$ vs. \f$ t\f$. It then determines the \f$R^2\f$
 * value for that set of coefficients and then iteratively improves the fit
 * using Newton's method.
 * @param initial Row number of the first data point to use
 * @param t Column matrix of time values [s]
 * @param Xdb Column matrix of moisture contents [kg/kg db]
 * @param Xe0 Initial guess for equilibrium moisture content [kg/kg db]
 * @returns Equilibrium moisture content [kg/kg db]
 *
 * @see regress
 */
double CalcXeIt(int initial, matrix *t, matrix *Xdb, double Xe0)
{
    int iter = 0, /* Keep track of the number of iterations */
        i; /* Loop index */
    double Xinit, /* Initial moisture content */
           Xe, /* Current guess for Xe */
           Xep, /* Previous guess for Xe */
           dR, /* First derivative of R^2 with respect to Xe */
           d2R, /* Second derivative of R^2 */
           tol = 1e-7, /* How close Xe and Xep need to be before we stop */
           h = 1e-10, /* Used for numerical differentiation */
           m = 2; /* Used to increase the rate of convergence */
    matrix *beta, /* Beta value from regress */
           *beta_ph, /* Same, but calculated at Xe + h */
           *beta_mh, /* Beta at Xe - h */
           *tmp_ph, /* Same as beta_ph, but with an extra zero to make rsquared happy */
           *tmp,
           *tmp_mh,
           *y, /* y values */
           *yph, /* Same as y, but at Xe + h */
           *ymh, /* y at Xe - h */
           *tadj, /* New matrix of t values that starts at the initial row */
           *Xadj; /* Same as tadj, but for Xdb */

    /* Set the initial moisture content */
    Xinit = val(Xdb, initial, 0);
    /* Set the first value of Xe to Xe0 */
    Xe = Xe0;

    /* Make smaller matricies that contain only the "good" data. */
    tadj = CreateMatrix(nRows(Xdb) - initial, 1);
    Xadj = CreateMatrix(nRows(Xdb) - initial, 1);
    for(i=initial; i<nRows(t); i++) {
        /* In addition to just copying the data, subtract the initial time from
         * each value to make sure that the intercept for the model goes through
         * the origin */
        setval(tadj, val(t, i, 0)-val(t, initial, 0), i-initial, 0);
        setval(Xadj, val(Xdb, i, 0), i-initial, 0);
    }

    /* Actually find Xe */
    do {
        /* Make a y matrix containing ln((Xdb - Xe)/(X0-Xe)) */
        y = CreateMatrix(nRows(Xadj), 1);
        for(i=0; i<nRows(Xadj); i++)
            setval(y, log((val(Xadj, i, 0) - Xe)/(Xinit - Xe)), i, 0);
        /* Calculate the kf parameter */
        beta = regress(y, tadj);

        /* Do the same, but at Xe - h */
        ymh = CreateMatrix(nRows(Xadj), 1);
        for(i=0; i<nRows(Xadj); i++)
            setval(ymh, log((val(Xadj, i, 0) - Xe-h)/(Xinit - Xe-h)), i, 0);
        beta_ph = regress(ymh, tadj);

        /* At Xe + h */
        yph= CreateMatrix(nRows(Xadj), 1);
        for(i=0; i<nRows(Xadj); i++)
            setval(yph, log((val(Xadj, i, 0) - Xe+h)/(Xinit - Xe+h)), i, 0);
        beta_mh = regress(yph, tadj);

        /* Add in a constant parameter of zero to the beta matrix. Do this for
         * each of the beta matricies we've got. */
        tmp_ph = CreateMatrix(2,1);
        setval(tmp_ph, 0, 0, 0);
        setval(tmp_ph, val(beta_ph, 0, 0), 1, 0);
        tmp = CreateMatrix(2,1);
        setval(tmp, 0, 0, 0);
        setval(tmp, val(beta, 0, 0), 1, 0);
        tmp_mh = CreateMatrix(2,1);
        setval(tmp_mh, 0, 0, 0);
        setval(tmp_mh, val(beta_mh, 0, 0), 1, 0);

        /* Calculate f and df */
        dR = (rsquared(tadj, yph, tmp_ph) - rsquared(tadj, ymh, tmp_mh))/(2*h);
        d2R = (rsquared(tadj, yph, tmp_ph) - 2*rsquared(tadj, y, tmp) + rsquared(tadj, ymh, tmp_mh))/(h*h);

        /* Calculate the new value of Xe */
        Xep = Xe;
        Xe = Xep + m*dR/d2R;

        /* Clean up */
        DestroyMatrix(y);
        DestroyMatrix(yph);
        DestroyMatrix(ymh);
        DestroyMatrix(beta);
        DestroyMatrix(tmp_ph);
        DestroyMatrix(tmp);
        DestroyMatrix(tmp_mh);

        /* Keep track of how many iterations we've gone through */
        iter++;

        /* Print out the current value */
        printf("Xe = %g\r", Xe);

        /* If Xe ever goes negative, admit defeat. */
        if(Xe < 0) {
            printf("Failure to converge after %d iterations.\n", iter);
            return Xe;
        }
    } while( fabs(Xe - Xep) > tol ); /* Check our value */

    /* Print out how many iterations it took to find Xe */
    printf("Solution converged after %d iterations.\n", iter);

    return Xe;
}


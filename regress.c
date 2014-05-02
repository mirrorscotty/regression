/**
 * @file regress.c
 * Set of functions for linear regression.
 */

#include <stdlib.h>
#include <math.h>

#include "regress.h"
#include "matrix.h"

/**
 * Equivalent of the Matlab "regress" function. Solves for the fitting
 * parameters using matrix algebra. Each column of the X matrix is a set of data
 * to used to fit a single parameter. To fit a constant, X should contain a
 * column of ones.
 * \f[
 * \underline{b} = (\underline{\underline{X}}^T\underline{\underline{X}})^{-1}
 *     \underline{\underline{X}}^T\underline{y}
 * \f]
 * @param y Column vector of dependent variable values
 * @param X Matrix of independent variable values, one variable per column
 * @returns Column matrix of fitted parameters. Each row corresponds to a 
 *      column in the supplies X matrix
 */
matrix* regress(matrix *y, matrix *X)
{
    matrix *Xt, *XtX, *XtXinv, *XtXinvXt, *beta;

    Xt = mtxtrn(X);
    XtX = mtxmul(Xt, X);
    XtXinv = CalcInv(XtX);
    XtXinvXt = mtxmul(XtXinv, Xt);

    beta = mtxmul(XtXinvXt, y);

    DestroyMatrix(Xt);
    DestroyMatrix(XtX);
    DestroyMatrix(XtXinv);
    DestroyMatrix(XtXinvXt);

    return beta;
}

/**
 * Matlab "polyfit" function. This fits the x-y data to a polynomial of
 * arbitrary order using the regress function.
 * @param x Column vector of dependent variable values
 * @param y Column vector of independent variable values
 * @param order Degree of the polynomial to fit to
 * @returns Column matrix of fitted parameters. Element n corresponds to the
 *      coefficient in front of x^n.
 */
matrix* polyfit(matrix* x, matrix* y, int order)
{
    matrix *X, *beta;
    int i, j, nelem;

    X = NULL;
    nelem = nRows(x);

    X = CreateMatrix(nelem, order+1);

    for(i=0; i<=nelem; i++) {
        for(j=0; j<=order; j++) {
            setval(X, pow(val(x, i, 0), j), i, j);
        }
    }

    beta = regress(y, X);
    DestroyMatrix(X);

    return beta;
}

/**
 * Calculate the coefficient of determination. This works only for output from
 * polyfit, or if the beta matrix supplied is of the same form. The \f$ R^2\f$
 * value is calculated using the following formula:
 * \f[
 * R^2 = 1-\frac{SSres}{SStot}
 * \f]
 * where \f$ SSres = \sum_0^N (y_i - f(x_i))^2 \f$ and
 * \f$ SStot = \sum_0^N (y_i - \overline{y})^2 \f$. The ybar value is the
 * average of the supplied y values.
 * @param x Column matrix of x values
 * @param y Column matrix of y values
 * @param beta Column matrix of fitting parameters
 * @returns R^2
 *
 * @see polyfit
 */
double rsquared(matrix* x, matrix* y, matrix *beta)
{
    double ybar = 0, /* Average y value */
           SStot = 0, /* Total sum of squares */
           SSres = 0, /* Residual sum of squares */
           f; /* Model value at xi */
    int i, j; /* Loop indicies */

    /* Calculate the average y value */
    for(i=0; i<nRows(y); i++)
        ybar += val(y, i, 0);
    ybar = ybar/nRows(y);

    /* Calculate the total sum of squares */
    for(i=0; i<nRows(y); i++)
        SStot += pow(val(y, i, 0) - ybar, 2);

    /* Calculate the residual sum of squares */
    for(i=0; i<nRows(y); i++) {
        /* Find the function value at xi. This assumes the function is of the
         * form f(x) = b0 + b1*x + b2*x^2 + ... + bn*X^n */
        f = 0;
        for(j=0; j<nRows(beta); j++)
            f += val(beta, j, 0) * pow(val(x, i, 0), j);

        SSres += pow(val(y, i, 0) - f, 2);
    }

    /* R^2 = 1-SSres/SStot */
    return 1-SSres/SStot;
}


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
    matrix *Y, *X, *beta;
    int i, j, nelem;

    X = NULL;
    nelem = nRows(x);

    X = CreateMatrix(nelem, order+1);

    for(i=0; i<=nelem; i++) {
        for(j=0; j<=order; j++) {
            setval(X, pow(val(x, i, 0), j), i, j);
        }
    }
    //Y = mtxtrn(y);

    beta = regress(y, X);
    DestroyMatrix(X);
    //DestroyMatrix(Y);

    return beta;
}
 

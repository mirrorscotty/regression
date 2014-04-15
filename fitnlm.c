/**
 * @file fitnlm.c
 * Non-linear least squares analysis using the Gauss-Newton algorithm.
 */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "regress.h"
#include "matrix.h"

/**
 * Calculate the J matrix used in the equations
 * \f[
 * J_{ij} = -\frac{\partial f_i}{\partial \beta_j}
 * \f]
 * @param model Model equation to fit to
 * @param x Column vector of independent variable values
 * @param beta Column matrix of fitting parameters
 * @returns Jacobian matrix
 */
matrix *CalcJacobian(double (*model)(double x, matrix *beta), matrix *x, matrix *beta)
{
    double h = 1e-10, /* Values used for numeric differentiation */
           Jij, /* i-jth element of the J matrix */
           xi; /* ith element of the supplied x values */
    matrix *J; /* Calculating this */
    int i, j, /* i is the row index for x, and j for beta */
        nx = nRows(x), /* Number of x values */
        nbeta = nRows(beta); /* Number of fitting parameters */

    /* Calculate all of the beta(i) + h values and store them for later */
    matrix** betah;
    betah = (matrix**) calloc(sizeof(matrix*), nbeta);
    for(j=0; j<nbeta; j++) {
        betah[j] = CopyMatrix(beta);
        addval(betah[j], h, j, 0);
    }

    J = CreateMatrix(nx, nbeta);
    for(i=0; i<nx; i++) {
        xi = val(x, i, 0);
        for(j=0; j<nbeta; j++) {
            /* Calculate the derivative of the model with respect to each
             * parameter. Each row is one x value. */
            Jij = (model(xi, betah[j]) - model(xi, beta)) / h;
            /* Save each Jij value into the J matrix */
            setval(J, Jij, i, j);
        }
    }
    
    /* Get rid of the betah crap we made earlier */
    for(j=0; j<nbeta; j++)
        DestroyMatrix(betah[j]);
    free(betah);

    return J;
}

/**
 * Calculate the difference between the predicted value and the actual value
 * \f[
 * \Delta y_i = y_i - f(x_i, \underline{\beta})
 * \f]
 * @param model Model equation to fit to
 * @param x Column vector of indep values
 * @param y Column vector of dependent values
 * @param beta Column vector of fitting parameters
 */
matrix* CalcDy(double (*model)(double x, matrix *beta), matrix *x, matrix *y, matrix* beta)
{
    matrix* dy;
    int i;
    double yi, xi; /* Individual x and y values */

    dy = CreateMatrix(nRows(x), 1);

    for(i=0; i<nRows(x); i++) {
        yi = val(y, i, 0); /* Current y value */
        xi = val(x, i, 0); /* Current x value */
        setval(dy, yi - model(xi, beta), i, 0); /* y - model(x, beta) */
    }

    return dy;
}

/**
 * Fit the given model to the x-y data provided
 * \f[
 * J_{ij}J_{is} \Delta\beta_{s} = J_{ij} \Delta y_i
 * \f]
 * @param model Equation to fit
 * @param x Column matrix of x values
 * @param y Column matrix of y values
 * @param beta0: Matrix of coefficients for the model
 * @returns Column vector of fitted coefficients
 */
matrix* fitnlm(double (*model)(double x, matrix *beta), matrix *x, matrix *y, matrix *beta0)
{
    matrix *J, /* Jacobian */
           *Jt, /* Transpose of J */
           *dy, /* dy value calculated in the function above */
           *A, *b, /* Coefficient and right hand side of the equation */
           *beta, /* Current values for the fitting parameters */
           *dbeta; /* Amount beta needs to change by after each iteration */
    /* Maximum amount of changed allowed for a single element of dbeta */
    double tol = .001;
    int i,
        maxiter = 500, /* Maximum number of iterations allowed */
        iter = 0; /* Current iteration */

    /* Make a copy of beta so we don't overwrite the supplied values */
    beta = CopyMatrix(beta0);
    /* Set dbeta to NULL so that the program doesn't segfault */
    dbeta = NULL;

    /* Loop until the change between iterations is less than the tolerance */
    do {
        /* Only delete dbeta if there's something to delete */
        if(dbeta)
            DestroyMatrix(dbeta);

        /* Calculate the values of each matrix. */
        dy = CalcDy(model, x, y, beta);
        J = CalcJacobian(model, x, beta);
        Jt = mtxtrn(J);
        A = mtxmul(Jt, J);
        b = mtxmul(Jt, dy);

        /* Solve the system of equations for how far off the fitting parameters
         * are. */
        dbeta = SolveMatrixEquation(A, b);

        /* beta = beta + dbeta */
        for(i=0; i<nRows(beta); i++)
            addval(beta, val(dbeta, i, 0), i, 0);

        /* Delete all the clutter we created */
        DestroyMatrix(J);
        DestroyMatrix(Jt);
        DestroyMatrix(A);
        DestroyMatrix(dy);
        DestroyMatrix(b);

        /* Check to see how many iterations we've gone through and quit if it
         * doesn't look like we're going to come up with an answer */
        if(iter++>maxiter) {
            printf("Maximum number of iterations reached, exiting.\n");
            break;
        }
    } while(fabs(mtxextrm(dbeta)) > tol); /* Check error */

    DestroyMatrix(dbeta); /* Get rid of the last unneeded matrix */

    return beta;
}


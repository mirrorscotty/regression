/**
 * @file fitburgers.c
 */

#include <stdio.h>
#include <math.h>

#include "regress.h"

#include "matrix.h"

#include "diffusivity.h"
#include "constants.h"
#include "isotherms.h"

/**
 * Modified version of the diffusion model from the Handbook of Food
 * Engineering.
 * This function is set up for use with fitnlm().
 * @param Xdb Moisture content [kg/kg db]
 * @param beta A 1x1 matrix containing the value of tau (tortuosity) [-]
 * @returns effective diffusivity
 */
double CreepModel(matrix *x, matrix *beta)
{
    double aM, aP, xi, J;
    double J0 = val(beta, 0, 0),
           J1 = val(beta, 1, 0),
           J2 = val(beta, 2, 0),
           l1 = val(beta, 3, 0),
           l2 = val(beta, 4, 0),
           mu0 = val(beta, 5, 0),
           aM0 = val(beta, 6, 0),
           M0 = val(beta, 7, 0),
           aP0 = val(beta, 8, 0),
           P0 = val(beta, 9, 0),

           t = val(x, 0, 0),
           M = val(x, 0, 1),
           P = val(x, 0, 2);

    aM = exp(aM0*(M-M0));
    aP = exp(aP0*(P-P0));
    xi = t*aM*aP;

    J = J0 + J1*(1-exp(-xi/l1)) + J2*(1-exp(-xi/l2)) + xi/mu0;
    return J;

}

/**
 * Load in a data file and calculate tortuosity.
 */
int main(int argc, char *argv[])
{
    matrix *data, *X, *y, *beta0, *beta;
    vector *t, *Xdb, *P;
    int tcol = 0, /* Column to get time from */
        xdbcol = 1, /* Column for Xdb */
        pcol = 2, /* Pressure column */
        jcol = 3; /* Creep compliance column */

    /* If a filename isn't supplied, spit out usage info and exit */
    if(argc != 2) {
        puts("Usage:");
        puts("fitburgers <datafile.csv>");
        return 0;
    }

    /* Load the csv file into a matrix */
    data = mtxloadcsv(argv[1], 0);

    /* Pull out the relevant data */
    t = ExtractColumnAsVector(data, tcol);
    Xdb = ExtractColumnAsVector(data, xdbcol);
    P = ExtractColumnAsVector(data, pcol);
    y = ExtractColumn(data, jcol);
    beta0 = ParseMatrix("[1.63e-6;1.45e-7;1.56e-7;2.282;25.78;1.42e9;-73;.14;1;2e5]");

    X = CatColVector(3, t, Xdb, P);

    beta = fitnlmM(&CreepModel, X, y, beta0);
    mtxprnt(beta);

    return 0;
}


/**
 * @file gab.c
 * Fit a set of water activity and moisture content data to the GAB equation
 * using nonlinear regression.
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "matrix.h"
#include "regress.h"

/**
 * GAB Equation suitable for the fitnlm function.
 * \f[
 * X_{db} = X_m \frac{C k a_w}{(1-k a_w)(1-k a_w + C k a_w)}
 * \f]
 * @param aw Water activity [-]
 * @param beta Column matrix of fitting parameters. Element 1 is C, element 2
 *      is k, and element 3 is Xm.
 * @returns Moisture content [kg/kg db]
 */
double oswin(matrix *X, matrix *beta)
{
    double k0 = val(beta, 0, 0),
           k1 = val(beta, 1, 0),
           n0 = val(beta, 2, 0),
           n1 = val(beta, 3, 0),

           aw = val(X, 0, 0),
           T = val(X, 0, 1),
    
           Xdb;
    printf("k0 = %g, k1 = %g, n0 = %g, n1 = %g, aw = %g, T = %g\n",
            k0, k1, n0, n1, aw, T);

    Xdb = (k0 + k1*T) * pow(aw/(1-aw), n0 + n1*T);

    return Xdb;
}

/**
 * Fit the Oswin parameters given water activity
 */
int main(int argc, char *argv[])
{
    double Tcol = 1,
           awcol = 3,
           Xdbcol = 2;
    matrix *data, *aw, *Xdb, *T, *beta0, *beta, *X;

    if(argc != 2) {
        puts("Usage:");
        printf("%s <aw.csv>\n", argv[0]);
        exit(0);
    }
    //data = mtxloadcsv("Andrieu.csv", 0);
    data = mtxloadcsv(argv[1], 1);

    /* Get the water activity from column 1 and the moisture content from
     * column 6. */
    aw = ExtractColumn(data, awcol);
    Xdb = ExtractColumn(data, Xdbcol);
    T = ExtractColumn(data, Tcol);

    X = AugmentMatrix(aw, T);
    mtxprnt(X);

    /* Use the parameters from Gina's thesis as a starting point */
    beta0 = CreateOnesMatrix(4, 1);
    setval(beta0, 0.1571, 0, 0);
    setval(beta0, -0.0012, 1, 0);
    setval(beta0, 0.2076, 2, 0);
    setval(beta0, 0.0043, 3, 0);

    /* Attempt to fit the gab parameters to the supplied data */
    beta = fitnlmM(&oswin, X, Xdb, beta0);

    /* Print out the fitted values */
    printf("k0 = %g\nk1 = %g\nn0 = %g\nn1 = %g\n",
            val(beta, 0, 0),
            val(beta, 1, 0),
            val(beta, 2, 0),
            val(beta, 3, 0));

    return 0;
}


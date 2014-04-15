/**
 * @file gab.c
 * Fit a set of water activity and moisture content data to the GAB equation
 * using nonlinear regression.
 */

#include <stdio.h>
#include "matrix.h"
#include "fitnlm.h"

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
double gab(double aw, matrix* beta)
{
    double C, k, Xm, Xdb;

    C = val(beta, 0, 0);
    k = val(beta, 1, 0);
    Xm = val(beta, 2, 0);
    
    Xdb = C*k*Xm*aw/((1-k*aw)*(1-k*aw+C*k*aw));

    return Xdb;
}

/**
 * Fit the GAB parameters given water activity
 */
int main(int argc, char *argv[])
{
    matrix *data, *aw, *Xdb, *tmp0, *tmp1, *beta0, *beta;

    if(argc != 2) {
        puts("Usage:");
        puts("gab <aw.csv>");
    }
    //data = mtxloadcsv("Andrieu.csv", 0);
    data = mtxloadcsv(argv[1], 0);

    /* Get the water activity from column 1 and the moisture content from
     * column 6. */
    aw = ExtractColumn(data, 0);
    Xdb = ExtractColumn(data, 5);

    /* Stick the two matricies together and delete any rows that contain
     * empty values */
    tmp0 = AugmentMatrix(aw, Xdb);
    tmp1 = DeleteNaNRows(tmp0);
    DestroyMatrix(aw);
    DestroyMatrix(Xdb);

    /* Pull the two columns back apart */
    aw = ExtractColumn(tmp1, 0);
    Xdb = ExtractColumn(tmp1, 1);
    DestroyMatrix(tmp0);
    DestroyMatrix(tmp1);

    /* Set up the beta matrix with some initial guesses at the GAB constants.
     * The solver needs these to be pretty close to the actual values, or it
     * will fail to converge */
    beta0 = CreateOnesMatrix(3, 1);
    setval(beta0, 6, 0, 0);
    setval(beta0, .5, 1, 0);
    setval(beta0, .04, 2, 0);

    /* Attempt to fit the gab parameters to the supplied data */
    beta = fitnlm(&gab, aw, Xdb, beta0);

    /* Print out the fitted values */
    printf("C = %g\nk = %g\nXm = %g\n",
            val(beta, 0, 0),
            val(beta, 1, 0),
            val(beta, 2, 0));

    return 0;
}


/**
 * @file fitburgers.c
 */

#include <stdio.h>
#include <math.h>

#include "regress.h"

#include "matrix.h"

#include "material-data.h"

/**
 * Load in a data file and calculate tortuosity.
 */
int main(int argc, char *argv[])
{
    int i;
    matrix *data, *X, *y, *beta0, *beta;
    vector *T, *Xdb;
    int xdbcol = 0, /* Column for Xdb */
        tempcol = 1, /* Temperature column */
        deffcol = 2; /* Effective diffusivity column */

    /* If a filename isn't supplied, spit out usage info and exit */
    if(argc != 2) {
        puts("Usage:");
        puts("fitachantadiff <datafile.csv>");
        return 0;
    }

    /* Load the csv file into a matrix */
    data = mtxloadcsv(argv[1], 1);

    /* Pull out the relevant data */
    T = ExtractColumnAsVector(data, tempcol);
    Xdb = ExtractColumnAsVector(data, xdbcol);
    y = ExtractColumn(data, deffcol);
    beta0 = ParseMatrix("[1.78e-5;36543.88;1.47e-10]");
    for(i=0; i<nRows(beta0); i++)
        setval(beta0, sqrt(val(beta0, i, 0)), i, 0);

    X = CatColVector(2, Xdb, T);

    beta = fitnlmM(&AchantaDiffModel, X, y, beta0);

    printf("D0: %g\nEa: %g\nDvap: %g\n",
            pow(val(beta, 0, 0), 2),
            pow(val(beta, 1, 0), 2),
            pow(val(beta, 2, 0), 2));

    DestroyMatrix(beta);
    DestroyMatrix(beta0);
    DestroyMatrix(data);
    DestroyMatrix(X);
    DestroyMatrix(y);

    DestroyVector(T);
    DestroyVector(Xdb);

    return 0;
}


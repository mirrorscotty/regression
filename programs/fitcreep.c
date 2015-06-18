#include "matrix.h"
#include "regress.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

int main(int argc, char *argv[])
{
    matrix *input, *X, *y, *b;
    vector *J, *tau;
    int i, j;
    double ti;

    if(argc < 3) {
        printf("Usage:\n"
               "fitcreep: <file> <t1> <t2> ... <tn>\n"
               "<file>: Filename containing the creep function data.\n"
               "<t1>: First retardation time\n"
               "<t2>: Second retardation time\n"
               "...\n"
               "<tn>: Nth retardation time.\n");
        exit(0);
    }

    /* The first row should probably be a header, and the second one might be
     * junk as well. */
    input = mtxloadcsv(argv[1], 2);

    tau = CreateVector(argc-2);
    J = CreateVector(argc-2);
    for(i=2; i<argc; i++)
        setvalV(tau, i-2, atof(argv[i]));

    y = ExtractColumn(input, 1);
    X = CreateMatrix(nRows(input), argc-1);

    for(i=0; i<nRows(input); i++) {
        ti = val(input, i, 0);
        setval(X, 1, i, 0);
        for(j=0; j<len(tau); j++)
            setval(X, 1-exp(-ti/valV(tau, j)), i, j+1);
    }

    b = regress(y, X);
    J = ExtractColumnAsVector(b, 0);

    PrintVector(tau);
    PrintVector(J);

    DestroyMatrix(input);
    DestroyMatrix(y);
    DestroyMatrix(X);
    DestroyMatrix(b);
    DestroyVector(J);
    DestroyVector(tau);

    return 0;
}


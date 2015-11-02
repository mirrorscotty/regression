#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "matrix.h"
#include "material-data.h"
#include "stress-strain.h"

int main(int argc, char *argv[])
{
    double e0, /* Imposed strain magnitude */
           s0, /* Measured stress magnitude */
           wmin, /* Strain oscillation frequency [1/s] */
           wmax, /* Strain oscillation frequency [1/s] */
           shift, /* Phase lag between stress and strain [-] */
           T, /* Temperature of material [K] */
           Xdb; /* Moisture content [kg/kg db] */
    maxwell *m; /* Set of Maxwell material parameters */
    matrix *beta; /* Coefficient matrix for regression */
    vector *frequency, *storage, *loss;
    matrix *output;
    int npts = 100, i;
    char *outfile;

    /* Print a usage statement if not enough arguments are supplied */
    if(argc != 6) {
        puts("Usage:");
        puts("modulus <e0> <w> <T> <Xdb>");
        puts("e0: Imposed strain magnitude");
        puts("w: Minimum frequency of strain oscillation");
        puts("w: Maximum frequency of strain oscillation");
        puts("T: Material temperature [K]");
        puts("Xdb: Material moisture content [kg/kg db]");

        exit(0);
    }

    /* Store command line parameters in variables */
    e0 = atof(argv[1]);
    wmin = atof(argv[2]);
    wmax = atof(argv[3]);
    T = atof(argv[4]);
    Xdb = atof(argv[5]);

    frequency = linspaceV(wmin, wmax, npts);
    storage = CreateVector(npts);
    loss = CreateVector(npts);

    for(i=0; i<npts; i++) {
        /* Fit the measured stress to the equation: s = s0 * sin(t*w+shift) */
        beta = fit_stress_rozzi(e0, valV(frequency, i), T, Xdb);

        /* Grab stress magnitude and phase lag from the coefficient matrix */
        s0 = val(beta, 0, 0);
        shift = val(beta, 1, 0);

        setvalV(storage, i, storage_mod(e0, s0, shift));
        setvalV(loss, i, loss_mod(e0, s0, shift));
    }

    outfile = (char*) calloc(sizeof(char), 20);
    sprintf(outfile, "output-%g-%g.csv", T, Xdb);

    output = CatColVector(3, frequency, storage, loss);
    mtxprntfilehdr(output, outfile, "freq(hz),storage,loss\n");

    return 0;
}


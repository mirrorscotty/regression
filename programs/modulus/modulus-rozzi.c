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
           w, /* Strain oscillation frequency [1/s] */
           shift, /* Phase lag between stress and strain [-] */
           T, /* Temperature of material [K] */
           Xdb; /* Moisture content [kg/kg db] */
    maxwell *m; /* Set of Maxwell material parameters */
    matrix *beta; /* Coefficient matrix for regression */

    /* Print a usage statement if not enough arguments are supplied */
    if(argc != 5) {
        puts("Usage:");
        puts("modulus <e0> <w> <T> <Xdb>");
        puts("e0: Imposed strain magnitude");
        puts("w: Frequency of strain oscillation");
        puts("T: Material temperature [K]");
        puts("Xdb: Material moisture content [kg/kg db]");

        exit(0);
    }

    /* Store command line parameters in variables */
    e0 = atof(argv[1]);
    w = atof(argv[2]);
    T = atof(argv[3]);
    Xdb = atof(argv[4]);

    /* Fit the measured stress to the equation: s = s0 * sin(t*w+shift) */
    beta = fit_stress_rozzi(e0, w, T, Xdb);

    /* Grab stress magnitude and phase lag from the coefficient matrix */
    s0 = val(beta, 0, 0);
    shift = val(beta, 1, 0);

    /* Calculate the storage and loss moduli and print them */
    printf("Storage Modulus: %g\nLoss Modulus: %g\n",
            storage_mod(e0, s0, shift), loss_mod(e0, s0, shift));

    return 0;
}


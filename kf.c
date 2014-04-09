/**
 * @file kf.c
 * Solve the Crank equation for kF using a nonlinear regression scheme.
 */

#include "matrix.h"
#include "fitnlm.h"
#include <math.h>
#include <stdio.h>

/**
 * Crank equation for a diffusion in a sheet.
 * \f[
 * \frac{X_{db}-X_e}{X_0-X_e}
 *     = \frac{8}{\pi^2}\sum_{n=0}^\infty\frac{1}{(2n+1)^2}
 *             \exp\left\{-k_F t (2n+1)^2\right\}
 * \f]
 * \f[
 * k_F = \frac{\pi^2 D}{l^2}
 * \f]
 * (Crank 1956)
 * @param kf Diffusivity constant (D*pi^2/l^2) where D is diffusivity and l is
 *      the slab thickness [1/s]
 * @param t Time [sec]
 * @param X0 Initial moisture content [kg/kg db]
 * @param Xe Equilibrium moisture content [kg/kg db]
 * @param nterms Number of terms of the equation to calculate
 * @returns Moisture content [kg/kg db]
 */
double CrankEquation(double kf, double t, double X0, double Xe, int nterms)
{
    double value = 0;
    int n;

    for(n=0; n<nterms; n++) {
        value += 8/(pow(2*n+1, 2) * M_PI*M_PI) *
            exp(-kf * t * pow(2*n+1, 2));
    }

    value = value * (X0-Xe) + Xe;

    return value;
}

/**
 * Function to allow the Crank equation to be used in the fitnlm function
 * @param t Time [s]
 * @param beta 1x1 matrix containing the value for kF
 * @returns Moisture content [kg/kg db]
 */
double CrankModel(double t, matrix *beta)
{
    double X0 = 0.3168046007, /* Initial moisture content */
           Xe = 0.1123,  /* Equilibrium moisture content */
           kf = val(beta, 0, 0); /* Get kF from the beta matrix */
    int nterms = 50; /* Number of terms to use */

    return CrankEquation(kf, t, X0, Xe, nterms);
}

/**
 * Fit a certain number of rows out of the supplied data
 * @param x Column matrix of x values
 * @param y Column matrix of y values
 * @param rowstart First row to use
 * @param rowend Last row to use
 * @returns Fitted value for kF [1/s]
 */
double fitsubset(matrix *x, matrix *y, int rowstart, int rowend)
{
    matrix *xx, /* Matrix to contain the x values of interest */
           *yy, /* Same for the y values */
           *beta, /* beta matrix for fitnlm */
           *beta0; /* Initial value for beta */
    int nrows = rowend-rowstart, /* Number of rows to fit */
        i; /* Loop index */

    /* Make all the matricies */
    xx = CreateMatrix(nrows, 1);
    yy = CreateMatrix(nrows, 1);
    beta0 = CreateMatrix(1, 1);

    /* Set the initial value of kF to 1e-4 */
    setval(beta0, 1e-4, 0, 0);

    /* Fill the xx and yy matricies with data from x and y */
    for(i=0; i<nrows; i++) {
        setval(xx, val(x, i+rowstart, 0), i, 0);
        setval(yy, val(y, i+rowstart, 0), i, 0);
    }

    /* Fit the data */
    beta = fitnlm(&CrankModel, xx, yy, beta0);

    /* Return the value for kF */
    return val(beta, 0, 0);
}

int main(int argc, char *argv[])
{
    matrix *data, *Xdb, *t, *kf;
    int tcol = 1, /* Column to get time from (in sec) */
        xdbcol = 4, /* Column for Xdb */
        chunksize = 500, /* Number of data points to fit per coefficient */
        rowstart = 0, /* First row to use */
        nchunks, i; /* Number of kF values to fit, loop index */

    /* If a filename isn't supplied, spit out usage info and exit */
    if(argc != 2) {
        puts("Usage:");
        puts("kF <datafile.csv>");
        return 0;
    }

    /* Load the csv file into a matrix */
    data = mtxloadcsv(argv[1], 0);

    /* Pull out the relevant data */
    Xdb = ExtractColumn(data, xdbcol);
    t = ExtractColumn(data, tcol);

    /* Truncated because of integer division */
    nchunks = (nRows(t)-rowstart)/chunksize;
    /* Make a matrix to store all of the kF data */
    kf = CreateMatrix(nchunks, 1);

    /* Fit kF using the Crank equation and nonlinear regresssion */
    for(i=0; i<nchunks; i++)
        setval(kf, fitsubset(t, Xdb, rowstart+i*chunksize, rowstart+(i+1)*chunksize), i, 0);

    /* Output the data to kF.csv */
    mtxprntfile(kf, "kF.csv");

    return 0;
}


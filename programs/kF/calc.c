/**
 * @file calc.c
 * Functions to calculate kF
 */

#include "kf.h"
#include "matrix.h"
#include "regress.h"

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
    setval(beta0, BETA0, 0, 0);

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

/**
 * Calculates kF using Newton's method
 *
 * TODO: Allow an initial point to be supplied.
 *
 * @param t Vector of times [s]
 * @param Xdb Vector of moisture contents [kg/kg db]
 * @param Xe Equilibrium moisture content [kg/kg db]
 * @returns Vector of values. Col 1: Time [s], Col 2: Moisture Content
 *      [kg/kg db], Col 3: kF [1/s]
 */
vector* calckf(vector *t, vector *Xdb, double Xe)
{
    vector *kF;
    int i; /* loop index */
    double kFi,
           X0 = valV(Xdb, 0),
           beta = BETA0;

    /* Make a matrix to store the resulting values in */
    kF = CreateVector(len(t));

    /* Calculate kF at each point directly using Newton's method */
    for(i=0; i<len(kF); i++) {
        beta = BETA0;
        kFi = CrankkF(valV(t, i), valV(Xdb, i), X0, Xe, beta);
        setvalV(kF, i, kFi);
        //beta = val(kF, i, 0);
    }

    return kF;
}

/**
 * Calculate kF using Newton's method, but resetting the time to zero and using
 * the previous data point as the initial moisture content. This should
 * (hopefully) make the kF data more accurate.
 * @param t Column matrix of times [s]
 * @param Xdb Column matrix of moisture contents [kg/kg db]
 * @param Xe Equilibrum moisture content [kg/kg db]
 * @param file Filename to save the data to
 * @returns Matrix of values. Col 1: Time [s], Col 2: Moisture Content
 *      [kg/kg db], Col 3: kF [1/s]
 */
matrix* calckfstep(matrix *t, matrix *Xdb, double Xe)
{
    matrix *kF, *data1, *data2;
    int i; /* loop index */
    double X0, /* Moisture content at the previous data point */
           kFi, /* Value of kF for that data point */
           dt; /* Difference in times between two points */

    /* Make a kF matrix */
    kF = CreateMatrix(nRows(t), 1);

    /* Set the values for X0 and dt. X0 will change with each loop iteration */
    X0 = val(Xdb, 0, 0);
    dt = val(t, 0, 0);

    /* Calculate all the kF values */
    for(i=0; i<nRows(kF); i++) {
        kFi = CrankkF(dt, val(Xdb, i, 0), X0, Xe, BETA0);
        setval(kF, kFi, i, 0);
        X0 = val(Xdb, i, 0);
    }

    /* Output data */
    data1 = AugmentMatrix(t, Xdb);
    data2 = AugmentMatrix(data1, kF);

    /* Clean up */
    DestroyMatrix(kF);
    DestroyMatrix(data1);

    return data2;
}

/**
 * Calculates kF using nonlinear regression. For this function, Xe must be set
 * in the macro at the top of the file. I'm far too lazy to rewrite my fitting
 * function to accept an additional parameter.
 * @param t Vector of time values [s]
 * @param Xdb Vector of moisture contents [kg/kg db]
 * @param file Filename to save the data to
 * @returns Matrix of values. Col 1: Time [s], Col 2: Moisture Content
 *      [kg/kg db], Col 3: kF [1/s]
 */
matrix* fitkf(matrix *t, matrix *Xdb)
{
    matrix *kf;
    int chunksize = 3, /* Number of data points to fit per coefficient */
        rowstart = 0, /* First row to use */
        nchunks, /* Number of kF values to fit */
        i, j; /* , loop indicies */
    double Xavg;

    /* Truncated because of integer division */
    nchunks = (nRows(t)-rowstart)/chunksize;
    /* Make a matrix to store all of the kF data */
    kf = CreateMatrix(nchunks, 2);

    /* Fit kF using the Crank equation and nonlinear regresssion */
    for(i=0; i<nchunks; i++) {
        Xavg = 0;
        for(j=rowstart+i*chunksize; j<rowstart+(i+1)*chunksize; j++)
            Xavg += val(Xdb, j, 0);
        Xavg = Xavg/chunksize;

        setval(kf, Xavg, i, 0);
        setval(kf, fitsubset(t, Xdb, rowstart+i*chunksize, rowstart+(i+1)*chunksize), i, 1);
    }

    /* Clean up */
    DestroyMatrix(Xdb);

    return kf;
}


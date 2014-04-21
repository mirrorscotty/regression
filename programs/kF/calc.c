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
 * @param t Column matrix of times [s]
 * @param Xdb column matrix of moisture contents [kg/kg db]
 * @param Xe Equilibrium moisture content [kg/kg db]
 * @param file Filename to save the data to
 */
void calckf(matrix *t, matrix *Xdb, double Xe, char *file)
{
    matrix *kF, *data1, *data2;
    int tcol = 1, /* Column to get time from (in sec) */
        xdbcol = 4, /* Column for Xdb */
        i; /* loop index */
    double X0 = val(Xdb, 0, 0);

    /* Make a matrix to store the resulting values in */
    kF = CreateMatrix(nRows(t), 1);

    /* Calculate kF at each point directly using Newton's method */
    for(i=0; i<nRows(kF); i++)
        setval(kF, CrankkF(val(t, i, 0), val(Xdb, i, 0), X0, Xe), i, 0);

    /* Print the results to a file */
    data1 = AugmentMatrix(t, Xdb);
    data2 = AugmentMatrix(data1, kF);
    mtxprntfile(data2, file);

    /* Clean up */
    DestroyMatrix(kF);
    DestroyMatrix(data1);
    DestroyMatrix(data2);
    
    return;
}

/**
 * Calculate kF using Newton's method, but resetting the time to zero and using
 * the previous data point as the initial moisture content. This should
 * (hopefully) make the kF data more accurate.
 * @param t Column matrix of times [s]
 * @param Xdb Column matrix of moisture contents [kg/kg db]
 * @param Xe Equilibrum moisture content [kg/kg db]
 * @param file Filename to save the data to
 */
void calckfstep(matrix *t, matrix *Xdb, double Xe, char *file)
{
    matrix *kF, *data1, *data2;
    int tcol = 1, /* Column to get time from (in sec) */
        xdbcol = 4, /* Column for Xdb */
        i; /* loop index */
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
        kFi = CrankkF(dt, val(Xdb, i, 0), X0, Xe);
        setval(kF, kFi, i, 0);
        X0 = val(Xdb, i, 0);
    }

    /* Output data */
    data1 = AugmentMatrix(t, Xdb);
    data2 = AugmentMatrix(data1, kF);
    mtxprntfile(data2, file);
    
    /* Clean up */
    DestroyMatrix(kF);
    DestroyMatrix(data1);
    DestroyMatrix(data2);

    return;
}

/**
 * Calculates kF using nonlinear regression. For this function, Xe must be set
 * in the macro at the top of the file. I'm far too lazy to rewrite my fitting
 * function to accept an additional parameter.
 * @param t Column matrix of time values [s]
 * @param Xdb Column matrix of moisture contents [kg/kg db]
 * @param file Filename to save the data to
 */
void fitkf(matrix *t, matrix *Xdb, char *file)
{
    matrix *kf;
    int tcol = 1, /* Column to get time from (in sec) */
        xdbcol = 4, /* Column for Xdb */
        chunksize = 3, /* Number of data points to fit per coefficient */
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

    /* Output the data to kF.csv */
    mtxprntfile(kf, file);

    /* Clean up */
    DestroyMatrix(kf);
    DestroyMatrix(Xdb);

    return;
}


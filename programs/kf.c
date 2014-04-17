/**
 * @file kf.c
 * Solve the Crank equation for kF using a nonlinear regression scheme.
 */

#include "matrix.h"
#include "regress.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define CONSTX0 0.3168046007
#define CONSTXe 0.1123
#define CONSTnterms 50
#define BETA0 1e-4

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

double CrankkF(double t, double X)
{
    double X0 = CONSTX0, /* Initial moisture content */
           Xe = CONSTXe,  /* Equilibrium moisture content */
           kf = BETA0,
           kfp = 0,
           f, df, h = 1e-10, tol = 1e-10;
    int nterms = CONSTnterms; /* Number of terms to use */

    do {
        f = CrankEquation(kf, t, X0, Xe, nterms) - X;
        df = (CrankEquation(kf+h, t, X0, Xe, nterms) 
                - CrankEquation(kf-h, t, X0, Xe, nterms))/(2*h);
        kfp = kf;
        kf = kf - f/df;
    } while(fabs(kfp - kf) > tol);

    return kf;
}

/**
 * Function to allow the Crank equation to be used in the fitnlm function
 * @param t Time [s]
 * @param beta 1x1 matrix containing the value for kF
 * @returns Moisture content [kg/kg db]
 */
double CrankModel(double t, matrix *beta)
{
    double X0 = CONSTX0, /* Initial moisture content */
           Xe = CONSTXe,  /* Equilibrium moisture content */
           kf = val(beta, 0, 0); /* Get kF from the beta matrix */
    int nterms = CONSTnterms; /* Number of terms to use */

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
 * @param filename File to load data from
 */
void calckf(char* filename)
{
    matrix *data, *Xdb, *t, *kF, *data1, *data2;
    int tcol = 1, /* Column to get time from (in sec) */
        xdbcol = 4, /* Column for Xdb */
        chunksize = 500, /* Number of data points to fit per coefficient */
        rowstart = 0, /* First row to use */
        nchunks, i; /* Number of kF values to fit, loop index */
    double Xavg;

    /* Load the csv file into a matrix */
    data = mtxloadcsv(filename, 0);

    /* Pull out the relevant data */
    Xdb = ExtractColumn(data, xdbcol);
    t = ExtractColumn(data, tcol);

    kF = CreateMatrix(nRows(t), 1);

    for(i=0; i<nRows(kF); i++)
        setval(kF, CrankkF(val(t, i, 0), val(Xdb, i, 0)), i, 0);

    data1 = AugmentMatrix(t, Xdb);
    data2 = AugmentMatrix(data1, kF);
    mtxprntfile(data2, "allkf.csv");
    
    return;
}

/**
 * Calculates kF using nonlinear regression
 * @param filename File to load data from
 */
void fitkf(char* filename)
{
    matrix *data, *Xdb, *t, *kf;
    int tcol = 1, /* Column to get time from (in sec) */
        xdbcol = 4, /* Column for Xdb */
        chunksize = 3, /* Number of data points to fit per coefficient */
        rowstart = 0, /* First row to use */
        nchunks, i, j; /* Number of kF values to fit, loop index */
    double Xavg;

    /* Load the csv file into a matrix */
    data = mtxloadcsv(filename, 0);

    /* Pull out the relevant data */
    Xdb = ExtractColumn(data, xdbcol);
    t = ExtractColumn(data, tcol);

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
    mtxprntfile(kf, "kF.csv");

    return;
}

/**
 * Calculate the equilibrium moisture content. This function determines the best
 * value of Xe to make a plot of \f$\ln\frac{X-X_e}{X_0-X_e}\f$ vs time linear.
 * In order to do this, it fits the data to the equation \f$y = a t + b\f$, where
 * \f$y = \ln(X-X_e)\f$ and \f$b = \ln(X_0-X_e)\f$, and then solves the equation
 * \f$F(X_e) = b - \ln(X_0-X_e) = 0\f$ using Newton's method.
 *
 * @param t Column matrix containing time during drying [s]
 * @param Xdb Column matrix of moisture content [kg/kg db]
 * @param Xe0 Initial guess for equilibrium moisture content.
 * @returns Equilibrium moisture content [kg/kg db]
 */
double CalcXe(matrix *t, matrix *Xdb, double Xe0)
{
    double f, df, /* Function values and derivatives */
           b, /* Fitting parameters. Only the constant matters */
           tol = 1e-10, /* Tolerance for Newton's method */
           Xe = Xe0, /* Set Xe to the initial guess */
           Xep, /* Previous guess */
           X0; /* Initial moisture content */
    matrix *beta, /* Matrix of fitting values */
           *y; /* Set equal to ln(X - Xe) */
    int i, /* Loop index */
        iter = 0;

    /* Set the initial moisture content */
    X0 = val(Xdb, 0, 0);

    do {
        /* Calculate b */
        y = CreateMatrix(nRows(Xdb), 1);
        for(i=0; i<nRows(Xdb); i++)
            setval(y, log(val(Xdb, i, 0) - Xe), i, 0);
        beta = polyfit(t, y, 1);
        b = val(beta, 0, 0);

        /* Calculate f and df */
        f = b - log(X0 - Xe);
        df = 1/(X0 - Xe);

        /* Calculate the new value of Xe */
        Xep = Xe;
        Xe = Xe - f/df;

        /* Clean up */
        DestroyMatrix(y);
        DestroyMatrix(beta);

        /* Keep track of how many iterations we've gone through */
        iter++;

        /* Print out the current value */
        printf("Xe = %g\r", Xe);
    } while( fabs(Xe - Xep) > tol ); /* Check our value */

    printf("Solution converged after %d iterations.\n", iter);

    return Xe;
}

matrix* LoadIGASorpTime(char *file)
{
    int row0 = 17, /* First row that contains numbers */
        col = 0, /* Get mass data from column 1 */
        i; /* Loop index */
    matrix *data, *min, *t;

    data = mtxloadcsv(file, row0);
    min = ExtractColumn(data, col);
    t = CreateMatrix(nRows(min), 1);

    for(i=0; i<nRows(min); i++)
        setval(t, val(min, i, 0)*60, i, 0);

    DestroyMatrix(min);
    DestroyMatrix(data);
    
    return t;
}

matrix* LoadIGASorpXdb(char *file, double Xdry)
{
    int row0 = 17, /* First row that contains numbers */
        col = 1, /* Get mass data from column 2 */
        i; /* Loop index */
    matrix *data, *M, *Xdb;

    data = mtxloadcsv(file, row0);
    M = ExtractColumn(data, col);
    Xdb = CreateMatrix(nRows(M), 1);

    for(i=0; i<nRows(M); i++)
        setval(Xdb, (val(M, i, 0)-Xdry)/Xdry, i, 0);

    DestroyMatrix(M);
    DestroyMatrix(data);
    
    return Xdb;
}

int main(int argc, char *argv[])
{
    matrix *t, *X;
    double Xe;

    /* If a filename isn't supplied, spit out usage info and exit */
    if(argc != 3) {
        puts("Usage:");
        puts("kF <datafile.csv> <Xdry>");
        puts("datafile.csv is the file to load data from.");
        puts("Xdry is the moisture content of the dry sample.");
        return 0;
    }

    t = LoadIGASorpTime(argv[1]);
    X = LoadIGASorpXdb(argv[1], atof(argv[2]));
    
    Xe = CalcXe(t, X, 0);
    printf("Xe = %g\n", Xe);

    //calckf(argv[1]);
    //fitkf(argv[1]);
}


#include "diffusivity.h"
#include <math.h>
#include "mechanical.h"
#include "matrix.h"
#include "choi-okos.h"
#include "constants.h"
#include "kf.h"

/**
 * Determine the lowest number row in kF such that all rows after it contain
 * positive values.
 * @param kF Column matrix of kF values.
 * @returns Row number
 */
int FindInitialPointkF(matrix *kF)
{
    int i;
    for(i=nRows(kF); i>0; i--)
        if(val(kF, i, 0) < 0)
            return i+1;
    return 0;
}

/**
 * Determine the lowest number row in RH such that all rows after it fall within
 * the tolerance.
 * @param RH Column matrix of relative humidity values
 * @returns Row number
 */
int FindInitialPointRH(matrix *RH)
{
    int i, /* Loop index */
        nrows = nRows(RH); /* Number of data points */
    double avg = 0, /* Set the average to zero initially */
           tol = 0.05; /* How close to the average we need to be */

    /* Calculate the average relative humidity */
    for(i=0; i<nrows; i++)
        avg += val(RH, i, 0);
    avg = avg/nrows;

    /* Find the lowest number row such that all rows after it fall within the
     * tolerance */
    for(i=nrows/2; i>0; i--)
        if(fabs(val(RH, i, 0) - avg) > tol)
            return i+1;

    return 0;
}

/**
 * Calculate the Deborah number. The characteristic length used to calculate
 * diffusion time scale is taken to be half the thickness of the slab.
 * @param initial Row number for the first data point to use.
 * @param point Row number for the current data point.
 * @param Xdb Column matrix of moisture contents [kg/kg db]
 * @param kF Column matrix of kF values [1/s]
 * @param L0 Initial slab thickness (full thickness) [m]
 * @param T Temperature (assumed constant) [K]
 * @param m Set of Maxwell parameters used for determining mean relaxation time.
 * @returns Deborah number [-]
 */
double DeborahNumber(int initial,
                     int point,
                     matrix *Xdb,
                     matrix *kF,
                     double L0,
                     double T,
                     maxwell* m)
{
    double Dkf0, /* Diffusivity at t=0 (calculated from kF values) */
           D0, /* Diffusivity at t=0 (calculated from model) */
           Di, /* Diffusivity at current time (from model) */
           D, /* Diffusivity (adjusted) */
           Li, /* Length at current time (from kF value) */
           Xdb0, /* Initial moisture content */
           kf0, /* Initial kF */
           Xdbi, /* Current moisture content */
           kFi, /* Current kF value */
           tD, /* Characteristic diffusion time */
           tr; /* Mean relaxation time */

    /* Grab the initial moisture content and kF values */
    Xdb0 = val(Xdb, initial, 0);
    kf0 = val(kF, initial, 0);

    /* Get the ones at the time we're interested in. */
    Xdbi = val(Xdb, point, 0);
    kFi = val(kF, point, 0);

    /* Calculate diffusivities */
    Dkf0 = kf0*L0*L0/(M_PI*M_PI);
    D0 = DiffCh10(Xdb0, T);
    Di = DiffCh10(Xdbi, T);

    /* Normalize the model diffusivity based on the initial diffusivity from the
     * kF value */
    D = Di/D0*Dkf0;

    /* Calculate the current length */
    Li = sqrt(M_PI*M_PI*D/kFi);
    Li = Li/2;

    /* Calculate characteristic diffusion time and mean relaxation time. */
    tD = Li*Li/D;
    tr = MeanRelaxTime(m);

    /* De = tr/tD */
    return tr/tD;
}

/**
 * Determine the current thickness of the sample based on kF value and
 * calculated diffusivity. This function returns the full thickness of the slab.
 * @param initial Row number for the first data point to use.
 * @param point Row number for the current data point.
 * @param Xdb Column matrix of moisture contents [kg/kg db]
 * @param kF Column matrix of kF values [1/s]
 * @param L0 Initial length [m]
 * @param T Temperature (assumed constant) [K]
 * @returns Thickness [-]
 */
double NewLength(int initial,
                 int point,
                 matrix *Xdb,
                 matrix *kF,
                 double L0,
                 double T)
{
    double Dkf0, /* Diffusivity at t=0 (calculated from kF values) */
           D0, /* Diffusivity at t=0 (calculated from model) */
           D, /* Diffusivity (adjusted) */
           Di, /* Diffusivity at current time (from model) */
           Li, /* Length at current time (from kF value) */
           Xdb0, /* Initial moisture content */
           kf0, /* Initial kF */
           Xdbi, /* Current moisture content */
           kFi; /* Current kF value */

    /* Grab the initial moisture content and kF values */
    Xdb0 = val(Xdb, initial, 0);
    kf0 = val(kF, initial, 0);

    /* Get the ones at the time we're interested in. */
    Xdbi = val(Xdb, point, 0);
    kFi = val(kF, point, 0);

    /* Calculate diffusivities */
    Dkf0 = kf0*L0*L0/(M_PI*M_PI);
    D0 = DiffCh10(Xdb0, T);
    Di = DiffCh10(Xdbi, T);

    /* Normalize the model diffusivity based on the initial diffusivity from the
     * kF value */
    D = Di/D0*Dkf0;

    /* Calculate the current length */
    Li = sqrt(M_PI*M_PI*D/kFi);

    return Li;
}

/**
 * Create a matrix of Deborah numbers, one of each data point.
 * @param Xdb Row matrix of moisture contents [kg/kg db]
 * @param kF Row matrix of kF values [1/s]
 * @param L0 Initial length [m]
 * @param T Drying temperature [K]
 * @param m Set of Maxwell material parameters.
 * @returns Column matrix of Deborah numbers.
 */
matrix* DeborahMatrix(int initial,
                      matrix *Xdb,
                      matrix *kF,
                      double L0,
                      double T,
                      maxwell* m)
{
    int i;

    matrix* De;

    De = CreateMatrix(nRows(kF), 1);
    for(i=0; i<initial; i++)
        setval(De, 0, i, 0);
    for(i=initial; i<nRows(De); i++)
        setval(De, DeborahNumber(initial, i, Xdb, kF, L0, T, m), i, 0);

    return De;
}

/**
 * Create a matrix of lengths, one for each data point. Each length is
 * calculated from the kF and diffusivity values at that point in time.
 * @param Xdb Row matrix of moisture contents [kg/kg db]
 * @param kF Row matrix of kF values [1/s]
 * @param L0 Initial length [m]
 * @param T Drying temperature [K]
 * @returns Column matrix of lengths
 */
matrix* LengthMatrix(int initial,
                     matrix *Xdb,
                     matrix *kF,
                     double L0,
                     double T)
{
    int i;

    matrix *L;

    L = CreateMatrix(nRows(kF), 1);
    for(i=0; i<initial; i++)
        setval(L, L0, i, 0);
    for(i=initial; i<nRows(L); i++)
        setval(L, NewLength(initial, i, Xdb, kF, L0, T), i, 0);

    return L;
}

/**
 * Calculate the length change due to water loss, assuming maximum shrinkage.
 * @param initial First row to look at
 * @param Xdb Column matrix of moisture contents [kg/kg db]
 * @param L0 Initial length [m]
 * @param Mdry Mass of the bone-dry sample
 * @param T Temperature at which the sample was dried.
 * @returns Column matrix of sample thicknesses
 */
matrix* LengthWaterLoss(int initial,
                        matrix *Xdb,
                        double L0,
                        double Mdry,
                        double T)
{
    double length = SLABLENGTH, /* Sample length [m] */
           width = SLABWIDTH, /* Sample width [m] */
           rhow, /* Density of water [kg/m^3] */
           Xdbi, /* Individual moisture content values */
           X0 = val(Xdb, initial, 0); /* Initial moisture content */
    int i; /* Loop index */
    choi_okos *co; /* Choi-okos values for water density */
    matrix *L; /* Calculated matrix of thicknesses */

    /* Calculate the density of water at this temperature */
    co = CreateChoiOkos(WATERCOMP);
    rhow = rho(co, T);
    DestroyChoiOkos(co);

    L = CreateMatrix(nRows(Xdb), 1);

    for(i=0; i<initial; i++)
        setval(L, L0, i, 0);
    for(i=initial; i<nRows(Xdb); i++) {
        Xdbi = val(Xdb, i, 0);
        Xdbi = Xdbi;

        setval(L, L0-(X0-Xdbi)*Mdry/(rhow * length*width), i, 0);
    }

    return L;
}

/**
 * Diffusivity calculations
 */
matrix* DOswinVector(int initial, matrix *X, double T)
{
    double i;
    matrix *D;
    D = CreateMatrix(nRows(X), 1);
    for(i=0; i<nRows(X); i++)
        setval(D, DiffCh10(val(X, i, 0), T), i, 0);
    return D;
}


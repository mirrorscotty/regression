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
 * @param kF Vector of kF values.
 * @returns Row number
 */
int FindInitialPointkF(vector *kF)
{
    int i;
    for(i=len(kF); i>0; i--)
        if(valV(kF, i) < 0)
            return i+1;
    return 0;
}

/**
 * Determine the lowest number row in RH such that all rows after it fall within
 * the tolerance.
 * @param RH Vector of relative humidity values
 * @returns Row number
 */
int FindInitialPointRH(vector *RH)
{
    int i, /* Loop index */
        nrows = len(RH); /* Number of data points */
    double avg = 0, /* Set the average to zero initially */
           tol = 0.05; /* How close to the average we need to be */

    /* Calculate the average relative humidity */
    for(i=0; i<nrows; i++)
        avg += valV(RH, i);
    avg = avg/nrows;

    /* Find the lowest number row such that all rows after it fall within the
     * tolerance */
    for(i=nrows/2; i>0; i--)
        if(fabs(valV(RH, i) - avg) > tol)
            return i+1;

    return 0;
}

/**
 * Calculate the Deborah number. The characteristic length used to calculate
 * diffusion time scale is taken to be half the thickness of the slab.
 * @param initial Row number for the first data point to use.
 * @param point Row number for the current data point.
 * @param Xdb Vector of moisture contents [kg/kg db]
 * @param kF Vector of kF values [1/s]
 * @param L0 Initial slab thickness (full thickness) [m]
 * @param T Temperature (assumed constant) [K]
 * @param m Set of Maxwell parameters used for determining mean relaxation time.
 * @returns Deborah number [-]
 */
double DeborahNumber(int initial,
                     int point,
                     vector *Xdb,
                     vector *kF,
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
    Xdb0 = valV(Xdb, initial);
    kf0 = valV(kF, initial);

    /* Get the ones at the time we're interested in. */
    Xdbi = valV(Xdb, point);
    kFi = valV(kF, point);

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
 * Create a matrix of Deborah numbers, one of each data point.
 * @param Xdb Row matrix of moisture contents [kg/kg db]
 * @param kF Row matrix of kF values [1/s]
 * @param L0 Initial length [m]
 * @param T Drying temperature [K]
 * @param m Set of Maxwell material parameters.
 * @returns Vector of Deborah numbers.
 */
vector* DeborahMatrix(int initial,
                      vector *Xdb,
                      vector *kF,
                      double L0,
                      double T,
                      maxwell* m)
{
    int i;

    vector* De;

    De = CreateVector(len(kF));
    for(i=0; i<initial; i++)
        setvalV(De, i, 0);
    for(i=initial; i<len(De); i++)
        setvalV(De, i, DeborahNumber(initial, i, Xdb, kF, L0, T, m));

    return De;
}

/**
 * Diffusivity calculations
 */
vector* DOswinVector(int initial, vector *X, double T)
{
    int i;
    vector *D;
    D = CreateVector(len(X));
    for(i=0; i<len(X); i++)
        setvalV(D, i, DiffCh10(valV(X, i), T));
    return D;
}


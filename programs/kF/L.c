#include "kf.h"
#include "diffusivity.h"
#include <math.h>
#include "mechanical.h"
#include "matrix.h"
#include "choi-okos.h"
#include "constants.h"

/**
 * @file L.c
 * Functions for calculating shrinkage.
 */

/**
 * Determine the current thickness of the sample based on kF value and
 * calculated diffusivity. This function returns the full thickness of the slab.
 * @param initial Row number for the first data point to use.
 * @param point Row number for the current data point.
 * @param Xdb Vector of moisture contents [kg/kg db]
 * @param kF Vector of kF values [1/s]
 * @param L0 Initial length [m]
 * @param T Temperature (assumed constant) [K]
 * @returns Thickness [-]
 */
double NewLength(int initial,
                 int point,
                 vector *Xdb,
                 vector *kF,
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

    return Li;
}

/**
 * Create a matrix of lengths, one for each data point. Each length is
 * calculated from the kF and diffusivity values at that point in time.
 * @param Xdb Row matrix of moisture contents [kg/kg db]
 * @param kF Row matrix of kF values [1/s]
 * @param L0 Initial length [m]
 * @param T Drying temperature [K]
 * @returns Vector of lengths
 */
vector* LengthMatrix(int initial,
                     vector *Xdb,
                     vector *kF,
                     double L0,
                     double T)
{
    int i;

    vector *L;

    L = CreateVector(len(kF));
    for(i=0; i<initial; i++)
        setvalV(L, i, L0);
    for(i=initial; i<len(L); i++)
        setvalV(L, i, NewLength(initial, i, Xdb, kF, L0, T));

    return L;
}

/**
 * @param Xdb Vector of moisture contents [kg/kg db]
 * @param kF Vector of kF values [1/s]
 * @param L0 Initial length [m]
 * @param T Drying temperature [K]
 * @returns Vector matrix of lengths
 */
vector* LengthConstD(int initial,
                     vector *kF,
                     double L0,
                     double T)
{
    int i;
    double Dkf0,
           kf0,
           kfi,
           Li;
    vector *L;

    /* Grab the initial kF and moisture content */
    kf0 = valV(kF, initial);

    /* Calculate initial diffusivity */
    Dkf0 = kf0*L0*L0/(M_PI*M_PI);

    L = CreateVector(len(kF));
    for(i=0; i<initial; i++)
        setvalV(L, i, L0);
    for(i=initial; i<len(L); i++) {
        kfi = valV(kF, i);

        Li = sqrt(M_PI*M_PI*Dkf0/kfi);
        
        setvalV(L, i, Li);
    }

    return L;
}

/**
 * Calculate the length change due to water loss, assuming maximum shrinkage.
 * @param initial First row to look at
 * @param Xdb Vector of moisture contents [kg/kg db]
 * @param L0 Initial length [m]
 * @param Mdry Mass of the bone-dry sample
 * @param T Temperature at which the sample was dried.
 * @returns Vector of sample thicknesses
 */
vector* LengthWaterLoss(int initial,
                        vector *Xdb,
                        double L0,
                        double Mdry,
                        double T)
{
    double length = SLABLENGTH, /* Sample length [m] */
           width = SLABWIDTH, /* Sample width [m] */
           rhow, /* Density of water [kg/m^3] */
           Xdbi, /* Individual moisture content values */
           X0 = valV(Xdb, initial); /* Initial moisture content */
    int i; /* Loop index */
    choi_okos *co; /* Choi-okos values for water density */
    vector *L; /* Calculated matrix of thicknesses */

    /* Calculate the density of water at this temperature */
    co = CreateChoiOkos(WATERCOMP);
    rhow = rho(co, T);
    DestroyChoiOkos(co);

    L = CreateVector(len(Xdb));

    for(i=0; i<initial; i++)
        setvalV(L, i, L0);
    for(i=initial; i<len(Xdb); i++) {
        Xdbi = valV(Xdb, i);

        setvalV(L, i, L0-(X0-Xdbi)*Mdry/(rhow * length*width));
    }

    return L;
}


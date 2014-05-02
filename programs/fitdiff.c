/**
 * @file fitdiff.c
 * Simple program to calculate tortuosity from diffusivity data, assuming that
 * the diffusivity constant can be written in terms of porosity, tortuosity,
 * the self-diffusion coefficient of water, and the binding energy of water.
 */

#include <stdio.h>
#include <math.h>

#include "regress.h"

#include "matrix.h"

#include "diffusivity.h"
#include "constants.h"
#include "isotherms.h"

/**
 * Modified version of the diffusion model from the Handbook of Food Engineering.
 * This function is set up for use with fitnlm().
 * @param Xdb Moisture content [kg/kg db]
 * @param beta A 1x1 matrix containing the value of tau (tortuosity) [-]
 * @returns effective diffusivity
 */
double DiffModel(double Xdb, matrix *beta)
{
    oswin *dat;
    dat = OSWINDATA();

    double Deff,
           T = 55+273.15, /* Temperature [K] */
           Dself = SelfDiffWater(T),
           phi = POROSITY,
           tau = val(beta, 0, 0),
           K = 1032.558, /* Source: Xiong et al. (1991) */
           Eb = BindingEnergyOswin(dat, Xdb, T),
           R = 8.314; /* Gas Constant */


    /* Equation 13 from Ch10 of Handbook of Food Engineering, Second Edition */
    Deff = phi/tau * Dself
        * ( K*exp(-Eb/(R*T)) / (1+K*exp(-Eb/(R*T))) );

    return Deff;
}

/**
 * Calculate the X matrix for use in the regress function.
 * @param Xdb Moisture content [kg/kg db]
 * @param T Temperature [K]
 * @returns Matrix for input into regress
 */
double CalcX(double Xdb, double T)
{
    oswin *dat;
    dat = OSWINDATA();

    double X,
           Dself = SelfDiffWater(T),
           K = 1032.558, /* Source: Xiong et al. (1991) */
           Eb = BindingEnergyOswin(dat, Xdb, T),
           R = 8.314; /* Gas Constant */


    /* Equation 13 from Ch10 of Handbook of Food Engineering, Second Edition */
    X = Dself * ( K*exp(-Eb/(R*T)) / (1+K*exp(-Eb/(R*T))) );

    return X;
}

/**
 * Load in a data file and calculate tortuosity.
 */
int main(int argc, char *argv[])
{
    matrix *data, *Xdb, *D, *beta, *X;
    int dcol = 1, /* Column to get diffusivity from */
        xdbcol = 0, /* Column for Xdb */
        i;
    double T = 55+273.15;

    /* If a filename isn't supplied, spit out usage info and exit */
    if(argc != 2) {
        puts("Usage:");
        puts("diffusivity <datafile.csv>");
        return 0;
    }

    /* Load the csv file into a matrix */
    data = mtxloadcsv(argv[1], 0);

    /* Pull out the relevant data */
    Xdb = ExtractColumn(data, xdbcol);
    D = ExtractColumn(data, dcol);
    //beta0 = CreateMatrix(1, 1);
    
    X = CreateMatrix(nRows(D), 1);

    for(i=0; i<nRows(X); i++) 
        setval(X, CalcX(val(Xdb, i, 0), T), i, 0);

    //setval(beta0, 3, 0, 0);

    //printf("Dinit = %g\n", DiffModel(.05, beta0));
    //beta = fitnlm(&DiffModel, Xdb, D, beta0);
    beta = regress(D, X);
    printf("phi/tau = %g\n", val(beta, 0, 0));
    printf("tau = %g\n", POROSITY/val(beta, 0, 0));

    return 0;
}


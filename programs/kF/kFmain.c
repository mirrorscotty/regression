/**
 * @file kFmain.c
 * Program to analyze drying data (primarily from the IGASorp) and calculate
 * diffusivity and shrinkage based on the Crank equation. Also calculates
 * several other quantities such as Deborah number and mass/momentum flux at
 * the surface of the sample.
 */

#include "kf.h"
#include "matrix.h"
#include "regress.h"
#include "mechanical.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[])
{
    vector *t, /* Time vector [s] */
           *X, /* Moisture content [kg/kg db] */
           *RH, /* Relative humidity [%] */
           *kF, /* kF value [-] */
           *L, /* Thickness (from kF + diffusivity model) [m] */
           *De, /* Deborah Number [-] */
           *Lwat, /* Thickness (from volume of water lost) [m] */
           *Diff, /* Diffusivity [m/s^2] */
           *MFlux,
           *MomeFlux,
           *Lconst; /* Thickness (from kF, constant diffusivity) [m] */
    matrix *data; /* Matrix of data for output */
    int p0; /* Initial data point */
    double Xe, /* Equilibrium moisture content [kg/kg db]*/
           Mdry, /* Mass of dry sample [g] */
           L0 = 6.22e-4, /* Initial slab thickness [m] */
           T = 60+273.15; /* Drying temperature [K] */
    char *outfile; /* Filename to output data to */
    maxwell *m; /* Set of maxwell data for calculating Deborah number */
    m = CreateMaxwell();
    /* Allocate memory for the output filename */
    outfile = (char*) calloc(sizeof(char), 80);

    /* If a filename isn't supplied, spit out usage info and exit */
    if(argc < 3) {
        puts("Usage:");
        puts("kF <datafile.csv> <Mdry> <Xe>");
        puts("datafile.csv: The file to load data from.");
        puts("Mdry: The mass of the dry sample. (in g)");
        puts("L0: Initial thickness (in mm)");
        puts("Xe: Optionally supply the equilibrium moisture content.");
        puts("");
        puts("Output is saved to kF<datafile.csv>.");
        return 0;
    }

    /* Pull the dry mass from the command line arguments */
    Mdry = atof(argv[2]);
    L0 = atof(argv[3])/1000; /* Convert to meters for calculations */

    /* Load all the important information from the IGASorp file.*/
    t = LoadIGASorpTime(argv[1]);
    X = LoadIGASorpXdb(argv[1], Mdry);
    RH = LoadIGASorpRH(argv[1]);

    /* Determine the first point to use for equilibrium moisture
     * content and similar calculations. Values will be calculated
     * for rows before this, but they should be disregarded.
     */
    p0 = FindInitialPointRH(RH);
    //p0 = 0;
    printf("Starting calculations from row %d.\n", p0);

    /* If equilibrium moisture content is supplied, use that value.
     * Otherwise, calculate Xe iteratively. In either case, print out the
     * value. */
    if(argc == 5)
        Xe = atof(argv[4]);
    else
        Xe = CalcXeIt(p0, t, X, valV(X, len(X)-1)*.95);
    printf("Xe = %g\n", Xe);

    /* Create the filename for the output csv file. It is always
     * kF prepended to the supplied input filename. */
    sprintf(outfile, "kF%s", argv[1]);

    /* Calculate kF (ratio of diffusivity to length squared) */
    kF = calckf(t, X, Xe);
    /* Determine the length from the kF value and the diffusivity model */
    L = LengthMatrix(p0, X, kF, L0, T);
    /* Deborah Number */
    De = DeborahMatrix(p0, X, kF, L0, T, m);
    /* Length change due solely to volume of water lost */
    Lwat = LengthWaterLoss(p0, X, L0, Mdry*1e-6, T);
    /* Length change from kF assuming constant diffusivity */
    Lconst = LengthConstD(p0, kF, L0, T);
    /* Diffusivity (from model) */
    Diff = DOswinVector(p0, X, T);
    MFlux = MassFlux(p0, t, X, Mdry*1e-6);
    //MomeFlux = MomentumFlux(p0, t, X, L, T, m);
    MomeFlux = PastaMassFlux(p0, t, L, L0, T);

    /* Combine all of the vectors for output */
    data = CatColVector(10, t, X, kF, L, De, Lwat, Lconst, Diff, MFlux, MomeFlux);

    /* Write the calculated values to a csv file. */
    mtxprntfilehdr(data, outfile, "Time [s],Moisture Content [kg/kg db],kF,Thickness [m],Deborah Number,Shrinkage (Water Loss),,Diffusivity,Mass Flux,Momentum Flux\n");
    DestroyMatrix(data);
    //DestroyMaxwell(m);

    return 0;
}


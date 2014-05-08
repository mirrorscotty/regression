/**
 * @file main.c
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
    vector *t, *X, *RH, *kF, *L, *De, *Lwat, *Diff, *MFlux, *MomeFlux;
    matrix *data;
    int p0; /* Initial data point */
    double Xe,
           Mdry,
           L0 = 6.22e-4,
           T = 60+273.15;
    char *outfile;
    maxwell *m;
    m = CreateMaxwell();
    outfile = (char*) calloc(sizeof(char), 80);

    /* If a filename isn't supplied, spit out usage info and exit */
    if(argc < 3) {
        puts("Usage:");
        puts("kF <datafile.csv> <Xdry>");
        puts("datafile.csv is the file to load data from.");
        puts("Mdry is the moisture content of the dry sample.");
        puts("Output is saved to kF<datafile.csv>.");
        return 0;
    }

    Mdry = atof(argv[2]);

    t = LoadIGASorpTime(argv[1]);
    X = LoadIGASorpXdb(argv[1], Mdry);
    RH = LoadIGASorpRH(argv[1]);

    p0 = FindInitialPointRH(RH);
    //p0 = 0;
    
    if(argc == 4)
        Xe = atof(argv[3]);
    else
        Xe = CalcXeIt(p0, t, X, .0402);
    printf("Xe = %g\n", Xe);

    sprintf(outfile, "kF%s", argv[1]);

    kF = calckf(t, X, Xe);
    L = LengthMatrix(p0, X, kF, L0, T);
    De = DeborahMatrix(p0, X, kF, L0, T, m);
    Lwat = LengthWaterLoss(p0, X, L0, Mdry*1e-6, T);
    Diff = DOswinVector(p0, X, T);
    MFlux = MassFlux(p0, t, X, Mdry*1e-6);
    //MomeFlux = MomentumFlux(p0, t, X, L, T, m);
    MomeFlux = PastaMassFlux(p0, t, L, L0, T);

    data = CatColVector(9, t, X, kF, L, De, Lwat, Diff, MFlux, MomeFlux);

    mtxprntfile(data, outfile);
    DestroyMatrix(data);

    return 0;
}


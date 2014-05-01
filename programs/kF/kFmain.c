/**
 * @file main.c
 * Simple program to solve the Crank equation for kF.
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
    matrix *t, *X, *RH, *kF, *data, *L, *De, *tmp;
    int p0; /* Initial data point */
    double Xe,
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
        puts("Xdry is the moisture content of the dry sample.");
        puts("Output is saved to kF<datafile.csv>.");
        return 0;
    }

    t = LoadIGASorpTime(argv[1]);
    X = LoadIGASorpXdb(argv[1], atof(argv[2]));
    RH = LoadIGASorpRH(argv[1]);

    p0 = FindInitialPointRH(RH);
    //p0 = 0;
    
    if(argc == 4)
        Xe = atof(argv[3]);
    else
        Xe = CalcXeIt(p0, t, X, .0402);
    printf("Xe = %g\n", Xe);

    t = LoadIGASorpTime(argv[1]);
    sprintf(outfile, "kF%s", argv[1]);
    data = calckf(t, X, Xe);
    kF = ExtractColumn(data, 2);

    L = LengthMatrix(p0, X, kF, L0, T);
    De = DeborahMatrix(p0, X, kF, L0, T, m);

    tmp = AugmentMatrix(data, L);
    DestroyMatrix(data);
    data = tmp;
    tmp = AugmentMatrix(data, De);
    DestroyMatrix(data);

    mtxprntfile(tmp, outfile);
    DestroyMatrix(tmp);

    return 0;
}


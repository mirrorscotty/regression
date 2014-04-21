/**
 * @file main.c
 * Simple program to solve the Crank equation for kF.
 */

#include "kf.h"
#include "matrix.h"
#include "regress.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>


int main(int argc, char *argv[])
{
    matrix *t, *X;
    double Xe;
    char *outfile;
    outfile = (char*) calloc(sizeof(char), 80);

    /* If a filename isn't supplied, spit out usage info and exit */
    if(argc != 3) {
        puts("Usage:");
        puts("kF <datafile.csv> <Xdry>");
        puts("datafile.csv is the file to load data from.");
        puts("Xdry is the moisture content of the dry sample.");
        puts("Output is saved to kF<datafile.csv>.");
        return 0;
    }

    t = LoadIGASorpTime(argv[1]);
    X = LoadIGASorpXdb(argv[1], atof(argv[2]));
    
    Xe = CalcXe(t, X, 0);
    printf("Xe = %g\n", Xe);

    t = LoadIGASorpTime(argv[1]);
    sprintf(outfile, "kF%s", argv[1]);
    calckf(t, X, Xe, outfile);
    //calckfstep(t, X, Xe, "kFstep.csv");
    //fitkf(t, X);
}


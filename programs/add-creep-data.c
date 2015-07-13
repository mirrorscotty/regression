#include "matrix.h"
#include "material-data.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

int main(int argc, char *argv[])
{

    int i;
    matrix *input, *output;
    double tcol=0,
           xcol=1,
           ucol=2,
           j0col=3,
           j1col=4,
           j2col=5,
           tau1col=6,
           tau2col=7,
           T, ti, xi, ui;
    char *femdata, *creepdata, *outfile;

    if(argc != 5) {
        printf("Usage:\n"
              "add-creep-data <femdata.csv> <creepdata.csv> <T> <outfile.csv>\n"
               "<femdata.csv>\tFile containing average moisture contents\n"
             "<creepdata.csv>\tTable of creep data at the desired temperature\n"
               "<T>\t\tTemperature (K)\n"
               "<outfile.csv>\tSave the new data here.\n");
        exit(0);
    }

    femdata = argv[1];
    creepdata = argv[2];
    T = atof(argv[3]);
    outfile = argv[4];
    
    input = mtxloadcsv(femdata, 1);
    output = CreateMatrix(nRows(input), 8);
    for(i=0; i<nRows(input); i++) {
        ti = val(input, i, tcol);
        xi = val(input, i, xcol);
        ui = val(input, i, ucol);

        setval(output, ti, i, tcol);
        setval(output, xi, i, xcol);
        setval(output, ui, i, ucol);
        setval(output, CreepLookupJ0(creepdata, T, xi), i, j0col);
        setval(output, CreepLookupJ1(creepdata, T, xi), i, j1col);
        setval(output, CreepLookupJ2(creepdata, T, xi), i, j2col);
        setval(output, CreepLookupTau1(creepdata, T, xi), i, tau1col);
        setval(output, CreepLookupTau2(creepdata, T, xi), i, tau2col);
    }

    mtxprntfilehdr(output, outfile, "t,Xdb,u,J0,J1,J2,tau1,tau2\n");

    return 0;
}


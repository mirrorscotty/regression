#include "matrix.h"
#include "material-data.h"
#include "regress.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

double PronyModel(double t, matrix* beta)
{
    matrix *Ji, *taui;
    double J0, J, Jval, tauval;
    int n, i;

    /* Pull out all the parameter values */
    J0 = val(beta, 0, 0);
    n = (nRows(beta)-1)/2;
    Ji = CreateMatrix(n, 1);
    taui = CreateMatrix(n, 1);

    for(i=0; i<n; i++) {
        setval(Ji, val(beta, 2*i+1, 0), i, 0);
        setval(taui, val(beta, 2*i+2, 0), i, 0);
    }

    J = J0*J0;
    for(i=0; i<n; i++) {
        Jval = val(Ji, i, 0);
        Jval = Jval*Jval;
        tauval = val(taui, i, 0);
        tauval = tauval*tauval;
        J += Jval * (1-exp(-t/tauval));
    }
    DestroyMatrix(Ji);
    DestroyMatrix(taui);
    return J;
}

matrix* makedata(matrix *t, double T, double M)
{
    matrix *J;
    int i;
    double Ji, ti;
    J = CreateMatrix(nRows(t), 1);
    for(i=0; i<nRows(t); i++) {
        ti = val(t, i, 0);
        Ji = LLauraCreep(ti, T, M, 0);
        setval(J, Ji, i, 0);
    }
    return J;
}

matrix* fitdata(matrix *t, matrix *J)
{
    matrix *beta0, *beta;
    beta0 = CreateMatrix(5, 1);
    setval(beta0, sqrt(1e-8), 0, 0);
    setval(beta0, sqrt(1e-8), 1, 0);
    setval(beta0, sqrt(10), 2, 0);
    setval(beta0, sqrt(1e-8), 3, 0);
    setval(beta0, sqrt(100), 4, 0);
    beta = fitnlm(&PronyModel, t, J, beta0);
    DestroyMatrix(beta0);
    return beta;
}

int main(int argc, char *argv[])
{
    int i, j, k;
    double Ti, Mj, percent;
    vector *T, *M;
    matrix *t, *Jij, *betaij, *output, *ttmp;

    /*
    if(argc < 3) {
        printf("Usage:\n"
               "fitcreep: <file> <t1> <t2> ... <tn>\n"
               "<file>: Filename containing the creep function data.\n"
               "<t1>: First retardation time\n"
               "<t2>: Second retardation time\n"
               "...\n"
               "<tn>: Nth retardation time.\n");
        exit(0);
    }
    */

    T = linspaceV(293, 363, 10);
    M = linspaceV(0, .5, 10);

    ttmp = linspace(1, 1e3, 1000);
    t = mtxtrn(ttmp);
    DestroyMatrix(ttmp);

    output = CreateMatrix(len(T)*len(M), 2+5);

    for(i=0; i<len(T); i++) {
        Ti = valV(T, i);
        for(j=0; j<len(M); j++) {
            Mj = valV(M, j);
            Jij = makedata(t, Ti, Mj);
            betaij = fitdata(t, Jij);

            setval(output, Ti, i*len(M)+j, 0);
            setval(output, Mj, i*len(M)+j, 1);
            for(k=0; k<nRows(betaij); k++)
                setval(output, pow(val(betaij, k, 0), 2), i*len(T)+j, k+2);
            DestroyMatrix(Jij);
            DestroyMatrix(betaij);

            /* Print the percent done */
            percent = (1.*i*len(M)+j)/(len(M)*len(T))*100.;
            printf("%3.2f %%\r", percent);
            fflush(stdout);
        }
    }
    
    DestroyMatrix(t);
    DestroyVector(T);
    DestroyVector(M);
    mtxprntfilehdr(output, "output.csv", "T,M,J0,J1,tau1,J2,tau2\n");
    DestroyMatrix(output);
    return 0;
}


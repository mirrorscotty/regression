#include "matrix.h"
#include "fitnlm.h"

double gab(double aw, matrix* beta)
{
    double C, K, Xm, Xdb;

    C = val(beta, 0, 0);
    K = val(beta, 1, 0);
    Xm = val(beta, 2, 0);
    
    Xdb = C*K*Xm*aw/((1-K*aw)*(1-K*aw+C*K*aw));

    return Xdb;
}

int main(int argc, char *argv[])
{
    matrix *data, *aw, *Xdb, *tmp0, *tmp1, *beta0, *beta;
    data = mtxloadcsv("Andrieu.csv", 0);

    aw = ExtractColumn(data, 0);
    Xdb = ExtractColumn(data, 5);

    tmp0 = AugmentMatrix(aw, Xdb);
    tmp1 = DeleteNaNRows(tmp0);
    DestroyMatrix(aw);
    DestroyMatrix(Xdb);

    aw = ExtractColumn(tmp1, 0);
    Xdb = ExtractColumn(tmp1, 1);

    beta0 = CreateOnesMatrix(3, 1);
    setval(beta0, 6, 0, 0);
    setval(beta0, .5, 1, 0);
    setval(beta0, .04, 2, 0);

    beta = fitnlm(&gab, aw, Xdb, beta0);
    mtxprnt(beta);
    return 0;
}


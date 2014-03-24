#include <stdlib.h>
#include <math.h>

#include "regress.h"
#include "matrix.h"

/* Equivalent of the Matlab "regress" function. */
matrix* regress(matrix *y, matrix *X)
{
    matrix *Xt, *XtX, *XtXinv, *XtXinvXt, *beta;

    Xt = mtxtrn(X);
    XtX = mtxmul(Xt, X);
    XtXinv = CalcInv(XtX);
    XtXinvXt = mtxmul(XtXinv, Xt);

    beta = mtxmul(XtXinvXt, y);

    DestroyMatrix(Xt);
    DestroyMatrix(XtX);
    DestroyMatrix(XtXinv);
    DestroyMatrix(XtXinvXt);

    return beta;
}

/* Matlab "polyfit" function. */
matrix* polyfit(matrix* x, matrix* y, int order)
{
    matrix *Y, *X, *beta;
    int i, j, nelem;

    X = NULL;
    nelem = nRows(x); /* Row matrix */

    X = CreateMatrix(nelem, order+1);

    for(i=0; i<=nelem; i++) {
        for(j=0; j<=order; j++) {
            setval(X, pow(val(x, i, 0), j), i, j);
        }
    }
    //Y = mtxtrn(y);

    beta = regress(y, X);
    DestroyMatrix(X);
    //DestroyMatrix(Y);

    return beta;
}
 

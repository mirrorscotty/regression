#include "kf.h"
#include "matrix.h"

/**
 * Load the time data from an IGASorp data file. The file needs to be converted
 * to a CSV file before loading. The header at the top of the file is ignored,
 * but the values must be separated by commas.
 * @param file The name of the file to open.
 * @returns A column matrix of times [s]
 */
matrix* LoadIGASorpTime(char *file)
{
    int row0 = 17, /* First row that contains numbers */
        col = 0, /* Get mass data from column 1 */
        i; /* Loop index */
    matrix *data, /* Raw data loaded from the file */
           *min, /* Time in minutes */
           *t; /* Time data in seconds */

    /* Load the data file and extract the column containing time */
    data = mtxloadcsv(file, row0);
    min = ExtractColumn(data, col);

    /* Convert from minutes to seconds */
    t = mtxmulconst(min, 60);

    /* Clean up */
    DestroyMatrix(data);
    DestroyMatrix(min);
    
    return t;
}

/**
 * Load the moisture content data from an IGASorp data file. The bone dry mass
 * must be supplied separately, and the file itself needs to be converted to a
 * CSV file before loading. The header at the top of the file is ignored, but
 * the values must be separated by commas.
 * @param file The name of the file to open.
 * @param Mdry The bone dry mass of the sample [mg]
 * @returns A column matrix of moisture content values [kg/kg db]
 */
matrix* LoadIGASorpXdb(char *file, double Mdry)
{
    int row0 = 17, /* First row that contains numbers */
        col = 1, /* Get mass data from column 2 */
        i; /* Loop index */
    matrix *data, /* Raw data from CSV file */
           *M, /* Mass of sample [mg] */
           *Xdb; /* Calculate moisture content [kg/kg db] */

    /* Load the data and pull out the relevant column */
    data = mtxloadcsv(file, row0);
    M = ExtractColumn(data, col);
    /* Make a matrix to hold the moisture content values */
    Xdb = CreateMatrix(nRows(M), 1);

    /* Calculate moisture content based on mass and bone dry mass */
    for(i=0; i<nRows(M); i++)
        setval(Xdb, (val(M, i, 0)-Mdry)/Mdry, i, 0);

    /* Clean up */
    DestroyMatrix(M);
    DestroyMatrix(data);
    
    return Xdb;
}



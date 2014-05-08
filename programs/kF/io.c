#include "kf.h"
#include "matrix.h"

/**
 * Load the time data from an IGASorp data file. The file needs to be converted
 * to a CSV file before loading. The header at the top of the file is ignored,
 * but the values must be separated by commas.
 * @param file The name of the file to open.
 * @returns A vector of times [s]
 */
vector* LoadIGASorpTime(char *file)
{
    int row0 = 17, /* First row that contains numbers */
        col = 0; /* Get mass data from column 1 */
    matrix *data; /* Raw data loaded from the file */
    vector *min, /* Time in minutes */
           *t; /* Time data in seconds */

    /* Load the data file and extract the column containing time */
    data = mtxloadcsv(file, row0);
    min = ExtractColumnAsVector(data, col);

    /* Convert from minutes to seconds */
    t = scalarmultV(60, min);

    /* Clean up */
    DestroyMatrix(data);
    DestroyVector(min);
    
    return t;
}

/**
 * Load the moisture content data from an IGASorp data file. The bone dry mass
 * must be supplied separately, and the file itself needs to be converted to a
 * CSV file before loading. The header at the top of the file is ignored, but
 * the values must be separated by commas.
 * @param file The name of the file to open.
 * @param Mdry The bone dry mass of the sample [mg]
 * @returns A vector of moisture content values [kg/kg db]
 */
vector* LoadIGASorpXdb(char *file, double Mdry)
{
    int row0 = 17, /* First row that contains numbers */
        col = 1, /* Get mass data from column 2 */
        i; /* Loop index */
    matrix *data; /* Raw data from CSV file */
    vector *M, /* Mass of sample [mg] */
           *Xdb; /* Calculate moisture content [kg/kg db] */

    /* Load the data and pull out the relevant column */
    data = mtxloadcsv(file, row0);
    M = ExtractColumnAsVector(data, col);
    /* Make a matrix to hold the moisture content values */
    Xdb = CreateVector(len(M));

    /* Calculate moisture content based on mass and bone dry mass */
    for(i=0; i<len(M); i++)
        setvalV(Xdb, i, (valV(M, i)-Mdry)/Mdry);

    /* Clean up */
    DestroyVector(M);
    DestroyMatrix(data);
    
    return Xdb;
}

/**
 * Load the time data from an IGASorp data file. The file needs to be converted
 * to a CSV file before loading. The header at the top of the file is ignored,
 * but the values must be separated by commas.
 * @param file The name of the file to open.
 * @returns A column matrix of times [s]
 */
vector* LoadIGASorpRH(char *file)
{
    int row0 = 17, /* First row that contains numbers */
        col = 2; /* Get humidity data from column 3 */
    matrix *data; /* Raw data loaded from the file */
    vector *RH;

    /* Load the data file and extract the column containing humidity */
    data = mtxloadcsv(file, row0);
    RH = ExtractColumnAsVector(data, col);

    /* Clean up */
    DestroyMatrix(data);
    
    return RH;
}


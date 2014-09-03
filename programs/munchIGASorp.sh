#!/bin/bash

# Convert the default output from the IGASorp program data files into a standard
# CSV file with a header, suitable for loading into Excel, or whatever.
# Usage:
# ./munchIGASorp.sh [IN.DAT] > [out.csv]

if [[ $# < 1 ]]
then
    echo 'Convert the default output from the IGASorp program data files into'
    echo 'a standard CSV file with a header, suitable for loading into Excel,'
    echo 'or whatever.'
    echo 'Usage:'
    echo './munchIGASorp.sh [IN.DAT] > [out.csv]'
else
    sed -e '/^[0-9]/s/\s\s*/,/g' $1 | sed -e '/^[0-9]/s/,$//'
fi


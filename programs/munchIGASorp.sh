#!/bin/bash

# Convert the default output from the IGASorp program data files into a standard
# CSV file with a header, suitable for loading into Excel, or whatever.
# Usage:
# ./munchIGASorp.sh [IN.DAT] > [out.csv]

sed -e '/^[0-9]/s/\s\s*/,/g' $1 | sed -e '/^[0-9]/s/,$//'


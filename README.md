regression
==========
* Author: Alex Griessman (alex.griessman@gmail.com)
* Repository: https://github.com/mirrorscotty/regression
This repository contains the following functions:
* regress: Linear regression
* polyfit: Linear regression for polynomial functions
* fitnlm: Nonlinear regression

Additionally, it has several programs to fit pasta drying parameters from data.
* `gab` - Fit a set of water activity and moisture content data to the GAB equation
    using nonlinear regression.
* `fitdiff` - Simple program to calculate tortuosity from diffusivity data,
    assuming that the diffusivity constant can be written in terms of porosity,
    tortuosity, the self-diffusion constant of water, and the binding energy of
    water.
* `kF` - Program to analyze drying data (primarily from the IGASorp) and calculate
    diffusivity and shrinkage based on the Crank equation. Also calculates
    several other quantities such as Deborah number and mass/momentum flux at
    the surface of the sample.
* `modulus` - Calculate the storage and loss moduli of a viscoelastic material
    given a set of Maxwell material properties as well as an imposed strain
    magnitude and frequency.
* `creep-table` - Generate a table of creep data at a specified temperature based
    on data from Rozzi (2002).

Building
--------
To compile any program, type `make <program>`


Dependencies
------------
Compilation requires GCC and GNU Make. This code also requires the matrix library
found [here](https://github.com/mirrorscotty/matrix), and the material data library
found [here](https://github.com/mirrorscotty/material-data).

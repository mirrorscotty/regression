regression
==========
This repository contains the following functions:
* regress: Linear regression
* polyfit: Linear regression for polynomial functions
* fitnlm: Nonlinear regression

Additionally, it has several programs to fit pasta drying parameters from data.
* gab: Fit a set of water activity and moisture content data to the GAB equation
    using nonlinear regression.
* fitdiff: Simple program to calculate tortuosity from diffusivity data,
    assuming that the diffusivity constant can be written in terms of porosity,
    tortuosity, the self-diffusion constant of water, and the binding energy of
    water.
* kF: Program to analyze drying data (primarily from the IGASorp) and calculate
    diffusivity and shrinkage based on the Crank equation. Also calculates
    several other quantities such as Deborah number and mass/momentum flux at
    the surface of the sample.
* modulus: Calculate the storage and loss moduli of a viscoelastic material
    given a set of Maxwell material properties as well as an imposed strain
    magnitude and frequency.


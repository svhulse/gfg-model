# Code for Host-pathogen coevolution promotes the evolution of general, broad-spectrum resistance and reduces foreign pathogen spillover risk

This repository contains the code used data referenced in the manuscipt comes from the code provided here. There are four python files principally used.

## Preprint availible here: https://www.biorxiv.org/content/10.1101/2023.08.04.548430v1.full.pdf

model.py is used to create a class instance which stores all the parameter values for a particular model.
This class structure is also used to calculate the mating matrix and transmission matrix.

solve.py is where the ODE model is defined, and the numerical solving performed. It takes a model object as an input, as well as initial conditions, and outputs an equilibria.

gen_raster.py is used to create a 2D raster of simulations. The parameters varied along the x and y axes can be set to any parameter, as well as the range of values each parameter takes. Once computed, gen_raster saves the output to a .pkl file

utilities.py contains helper functions, mostly used to unpack and analyse the output from gen_raster. It also contains the code needed to calculate the transitivity slope for a particular equilibrium.

### The parameters used for each raster scenario are stored in rasters.json, the given scenarios are described below:

nocov_gs: No coevolution, fixed intermediate recombination, varied costs of general and specific resistance

cov_gs: Coevolution, fixed intermediate recombination, varied costs of general and specific resistance

cov_gv: Coevolution, fixed intermediate recombination, varied costs of general resistance and virulence costs

nocov_gs_rho0: No coevolution, variable recombination with allele fixed for linkage modifier, varied costs of general and specific resistance

cov_gs_rho0: Coevolution, variable recombination with allele fixed for linkage modifier, varied costs of general and specific resistance

cov_gv_rho0: Coevolution, variable recombination with allele fixed for linkage modifier, varied costs of general resistance and virulence costs

nocov_gs_rho1: No coevolution, variable recombination with loss of linkage modifier, varied costs of general and specific resistance

cov_gs_rho1: Coevolution, variable recombination with loss of linkage modifier, varied costs of general and specific resistance

cov_gv_rho1: Coevolution, variable recombination with loss of linkage modifier, varied costs of general resistance and virulence costs


nocov_gs_hr: No coevolution, fixed high recombination, varied costs of general and specific resistance

cov_gs_hr: Coevolution, fixed high recombination, varied costs of general and specific resistance

cov_gv_hr: Coevolution, fixed high recombination, varied costs of general resistance and virulence costs


nocov_gs_hs: No coevolution, fixed intermediate recombination, varied costs of general and specific resistance, hard selection model

cov_gs_hs: Coevolution, fixed intermediate recombination, varied costs of general and specific resistance, hard selection model

cov_gv_hs: Coevolution, fixed intermediate recombination, varied costs of general resistance and virulence costs, hard selection model


nocov_gs_sg: No coevolution, fixed intermediate recombination, varied costs of general and specific resistance, strong general resistance

cov_gs_sg: Coevolution, fixed intermediate recombination, varied costs of general and specific resistance, strong general resistance

cov_gv_sg: Coevolution, fixed intermediate recombination, varied costs of general resistance and virulence costs, strong general resistance
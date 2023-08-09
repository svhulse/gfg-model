# Code for Host-pathogen coevolution promotes the evolution of general, broad-spectrum resistance and reduces foreign pathogen spillover risk

Preprint availible here: https://www.biorxiv.org/content/10.1101/2023.08.04.548430v1.full.pdf

All the data referenced in the manuscipt comes from the code provided here. There are four python files principally used.

model.py is used to create a class instance which stores all the parameter values for a particular model.
This class structure is also used to calculate the mating matrix and transmission matrix.

solve.py is where the ODE model is defined, and the numerical solving performed. It takes a model object as an input, as well as initial conditions, and outputs an equilibria.

gen_raster.py is used to create a 2D raster of simulations. The parameters varied along the x and y axes can be set to any parameter, as well as the range of values each parameter takes. Once computed, gen_raster saves the output to a .pkl file

utilities.py contains helper functions, mostly used to unpack and analyse the output from gen_raster. It also contains the code needed to calculate the transitivity slope for a particular equilibrium.
# Host-pathogen coevolution promotes the evolution of general, broad-spectrum resistance and reduces foreign pathogen spillover risk

This repository contains the code used to generate the data referenced in our manuscipt. The code used to make each figure is also included, although most figures require several data rasters to be generated first. All of this code is written assuming the file structure of this repository.

## Code Organisation

The code used to generate data consists of model.py, solve.py, and gen_raster.py.

model.py is used to create a class instance which stores all the parameter values for a particular model.
This class structure is also used to calculate the mating matrix and transmission matrix.

solve.py is where the ODE model is defined, and the numerical solving performed. It takes a model object as an input, as well as initial conditions, and outputs solutions. This code contains two methods: get_sol, which returns equilibrium points, and run_sim, which returns trajectories. It takes in Model class objects to store the parameters for each simulation.

gen_raster.py is used to create a 2D raster of simulations. The parameters varied along the x and y axes can be set to any parameter, as well as the range of values each parameter takes. Once computed, gen_raster saves the output to a .p file

## Plotting Scripts

The other code in this repository is used for analysing and plotting, and includes utilities.py, style.py, fig_1.py,...,fig_S4.py.

utilities.py contains helper functions, mostly used to unpack and analyse the output from gen_raster. It also contains the code needed to calculate the transitivity slope for a particular equilibrium.

style.py defines the graphics themes used for all plots.

Finally, the figure scripts generate and save figures as .svg images. All figures that display data rasters require the data to be generated prior to running the script.

## Proceedure to Generate Plots

To recreate the data for the manuscript, first download this repository. You will then need to use the gen_raster script to compute solution rasters. The parameter values used for all rasters are stored in the rasters.json file. This file can also be used to define additional raster scenarios. An example scenario is given below.

```
"cov_gs" : {
    "filename"  : "./data/cov_gs.p",
    "var_1"     : "c_g",
    "var_2"     : "c_s",
    "S_init"    : [1, 0.1, 0.1],
    "I_init"    : 0.9,
    "params"    : {
        "k":0.001,
		"mu":0.2,
		"b":1,
		"beta":0.5,
		"nh":0.1,
		"g":0.3,
		"s":0.9,
		"rho":[0.05, 0.05],
		"c_g":0.1,
		"c_s":0.2,
		"v":0.2,						
		"sel":"soft"
    }
},
```

Here, the first line defines the name of the scenario, var_1 refers to the the parameter which will be varied along the x-axis, and var_2 refers to the variable that will be varied along the y-axis.

S_init and I_init define the initial allele frequencies for the host and pathogen respectively. For the host, the first number is the frequency of the linkage modifer allele, the second is the general resistance allele and the third is the specific resistance allele. I_init refers to the initial frequency of the Avr genotype.

Finally, the values in the params list correspond to the simulation parameters in the Model class. For the recombination rate, rho[0] is the recombination rate with the linkage modifier, and rho[1] is the recombination rate without.

To generate a raster, set the scenario value to the name of the raster scenario desired.

To make each figure, you will need to generate the following rasters:
**Figure 2**: nocov_gs, cov_gs, cov_gv
**Figure 3**: nocov_gs_rho0, cov_gs_rho0
**Figure 4**: cov_gs, cov_gs_rho0
**Figure 5**: cov_gs_rho0, cov_gs_rho1
**Figure S1**: nocov_gs, cov_gs, cov_gv
**Figure S2**: cov_gs_rho1, cov_gv_rho1, cov_gs_rho0, cov_gv_rho0
**Figure S3**: cov_gs, cov_gv
**Figure S4**: all scenarios listed below

Once all necesary rasters are made, you should be able to run the figure scripts, which will save the images to the figures folder.

### Predefined Scenarios
The parameters used for each raster scenario are stored in rasters.json, the given scenarios are described below:

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
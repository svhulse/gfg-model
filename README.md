# Host-pathogen coevolution promotes the evolution of general, broad-spectrum resistance and reduces foreign pathogen spillover risk

This repository contains the code used to generate the data referenced in our manuscript. The code used to make each figure is also included, although most figures require several data rasters to be generated first. All of this code is written assuming the file structure of this repository. Unfortunately, the data files themselves are too large to be stored in this repository, and only low-resolution data for main-text figures is provided here. Higher resolution raster will need to be generated on the user's end.

## Code Organization

The code used to generate data consists of `model.py`, `solve.py`, and `gen_raster.py`.

`model.py` is used to create a class instance which stores all the parameter values for a particular model.
This class structure is also used to calculate the mating matrix and transmission matrix.

`solve.py` is where the ODE model is defined, and the numerical solving performed. It takes a model object as an input, as well as initial conditions, and outputs solutions. This code contains two methods: get_sol, which returns equilibrium points, and run_sim, which returns trajectories. It takes in Model class objects to store the parameters for each simulation.

`gen_raster.py` is used to create a 2D raster of simulations. The parameters varied along the x and y axes can be set to any parameter, as well as the range of values each parameter takes. Once computed, gen_raster saves the output to a .p file

## Plotting Scripts

The other code in this repository is used for analyzing and plotting, and includes `utilities.py`, `style.py`, `fig_1.py`,...,`fig_S4.py`.

`utilities.py` contains helper functions, mostly used to unpack, and analyze the output from gen_raster. It also contains the code needed to calculate the transitivity slope for a particular equilibrium.

`style.py` defines the graphics themes used for all plots.

Finally, the figure scripts generate and save figures as .svg images. All figures that display data rasters require the data to be generated prior to running the script.

## Procedure to Generate Plots

To recreate the data for the manuscript, the procedure is:
* Download this repository
* Run the gen_raster.py script with the desired scenarios
* Run the figure desired figure scripts

While many figures used saved raster files, Fig. 1 only uses individual simulations and can be run without generating rasters. 

Note that gen_raster is set up to use multithreading, with a default value of 4 cores. This can be changed by altering the cores value at line 20. The parameter values used for all rasters are stored in the rasters.json file. This file can also be used to define additional raster scenarios. An example scenario is given below.

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

Here, the first line defines the name of the scenario, `var_1` refers to the parameter which will be varied along the x-axis, and `var_2` refers to the variable that will be varied along the y-axis.

`S_init` and `I_init` define the initial allele frequencies for the host and pathogen respectively. For the host, the first number is the frequency of the linkage modifier allele, the second is the general resistance allele and the third is the specific resistance allele. `I_init` refers to the initial frequency of the *Avr* genotype.

Finally, the values in the params list correspond to the simulation parameters in the Model class. For the recombination rate, `rho[0]` is the recombination rate with the linkage modifier, and `rho[1]` is the recombination rate without. For a detailed description of all parameters, please see the manuscript.

To generate a raster, set the scenario value to the name of the raster scenario desired. Raster scenarios, as defined in `rasters.json` and the raster size can be specified from the command line by running the code as follows:

```
python gen_raster.py cov_gs 200
```

Where the first argument is the name of the scenario and the second is the size. If no size is provided, a default value of 200 will be used.

To make each figure, you will need to generate the following rasters:

**Figure 2**: nocov_gs, cov_gs, cov_gv\
**Figure 3**: nocov_gs_rho0, cov_gs_rho0\
**Figure 4**: cov_gs, cov_gs_rho0\
**Figure 5**: cov_gs_rho0, cov_gs_rho1\
**Figure S1**: nocov_gs, cov_gs, cov_gv\
**Figure S2**: cov_gs_rho1, cov_gv_rho1, cov_gs_rho0, cov_gv_rho0\
**Figure S3**: cov_gs, cov_gv\
**Figure S4**: all scenarios listed below

Once all necessary rasters are made, you should be able to run the figure scripts, which will save the images to the figures folder.

### Predefined Scenarios
`nocov_gs`: no coevolution, fixed intermediate recombination, varied costs of general and specific resistance\
`cov_gs`: coevolution, fixed intermediate recombination, varied costs of general and specific resistance\
`cov_gv`: coevolution, fixed intermediate recombination, varied costs of general resistance and virulence costs\
`nocov_gs_rho0`: no coevolution, variable recombination with allele fixed for linkage modifier, varied costs of general and specific resistance\
`cov_gs_rho0`: coevolution, variable recombination with allele fixed for linkage modifier, varied costs of general and specific resistance\
`cov_gv_rho0`: coevolution, variable recombination with allele fixed for linkage modifier, varied costs of general resistance and virulence costs\
`nocov_gs_rho1`: no coevolution, variable recombination with loss of linkage modifier, varied costs of general and specific resistance\
`cov_gs_rho1`: coevolution, variable recombination with loss of linkage modifier, varied costs of general and specific resistance\
`cov_gv_rho1`: coevolution, variable recombination with loss of linkage modifier, varied costs of general resistance and virulence costs\
`nocov_gs_hr`: no coevolution, fixed high recombination, varied costs of general and specific resistance\
`cov_gs_hr`: coevolution, fixed high recombination, varied costs of general and specific resistance\
`cov_gv_hr`: coevolution, fixed high recombination, varied costs of general resistance and virulence costs\
`nocov_gs_hs`: no coevolution, fixed intermediate recombination, varied costs of general and specific resistance, hard selection model\
`cov_gs_hs`: coevolution, fixed intermediate recombination, varied costs of general and specific resistance, hard selection model\
`cov_gv_hs`: coevolution, fixed intermediate recombination, varied costs of general resistance and virulence costs, hard selection model\
`nocov_gs_sg`: no coevolution, fixed intermediate recombination, varied costs of general and specific resistance, strong general resistance\
`cov_gs_sg`: coevolution, fixed intermediate recombination, varied costs of general and specific resistance, strong general resistance\
`cov_gv_sg`: coevolution, fixed intermediate recombination, varied costs of general resistance and virulence costs, strong general resistance

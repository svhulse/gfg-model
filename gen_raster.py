import numpy as np
import tqdm
import json 
import pickle as pkl
import multiprocessing as mp

from model import Model
from solve import get_sol

scenario = 'cov_gv_rho0'					#Name of raster scenario
size = 200									#Raster dimension

with open('rasters.json', 'r') as data:
	param_set = json.load(data)[scenario]

output_path = param_set['filename']			#Output filename

var_1 = param_set['var_1']					#First parameter rastered
var_2 = param_set['var_2']					#Second parameter rastered

'''
Range of values for raster variables, too vary something other than v, c_g, or c_s, you 
will need to add a new range array and add it to the vars dict below
'''
V_costs = np.linspace(0, 0.3, size)			#Range of parameters for virulence costs
G_costs = np.linspace(0, 0.2, size)			#Range of parameters for general resistance costs
S_costs = np.linspace(0, 0.4, size)			#Range of parameters for speific resistance costs

S_init = param_set['S_init']				#Initial host allele frequencies (Recomb, General, Specific)
I_init = param_set['I_init']				#Initial proportion of the Avr pathogen genotype

params = param_set['params']

vars = {'c_g': G_costs, 'c_s': S_costs, 'v': V_costs}

def pass_to_sim(model):
	return get_sol(model, S_init, I_init)

if __name__ == '__main__':
	coords = []     #x, y coordinates of each simulation in raster
	models = []     #Empty tuple for model classes

	#Create raster of model classes for each parameter combination
	for i in range(size):
		for j in range(size):
			coords.append((i,j))

			params[var_1] = vars[var_1][i]
			params[var_2] = vars[var_2][j]
			new_model = Model(**params)

			models.append(new_model)
		
	#Run simluations for 4 core processor
	pool = mp.Pool(processes=4)	
	
	results = []
	for result in tqdm.tqdm(pool.imap(pass_to_sim, models), total=len(models)):
		results.append(result)

	raster = []
	for i in range(size):
		inds = [j for j in range(len(coords)) if coords[j][1] == i]
		raster.append([results[j] for j in inds])

	with open(output_path, 'wb') as f:
		pkl.dump([models, results], f)
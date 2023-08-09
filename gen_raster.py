import numpy as np
import tqdm
import pickle as pkl
import multiprocessing as mp

from model import Model
from solve import run_sim

output_path = './data/cov_gs.p'				#Output filename
size = 100									#Raster dimension

var_1 = 'c_g'								#First parameter rastered
var_2 = 'c_s'								#Second parameter rastered

S_init = [0, 0.1, 0.5]						#Initial host allele frequencies (Recomb, General, Specific)
I_init = 0.9								#Avr proportion

params = {  'k':0.001,						#Coefficient of density dependent growth
			'mu':0.2,						#Deathrate
			'b':1,							#Birthrate
			'beta':0.5,						#Baseline transmission rate
			'nh':0.1,						#Nonhost transmission rate
			'g':0.3,						#Strength of general resistance
			's':0.9,						#Strength of specific resistance
			'rho':[0.05, 0.05],				#Recombination rate for each allele
			'c_g':0.1,						#Cost of general resistance
			'c_s':0.2,						#Cost of specific resistance
			'v':0.2,
			'sel':'soft'}					#Cost of virulence      

V_costs = np.linspace(0, 0.3, size)			#Range of parameters for virulence costs
G_costs = np.linspace(0, 0.2, size)			#Range of parameters for general resistance costs
S_costs = np.linspace(0, 0.4, size)			#Range of parameters for speific resistance costs

vars = {'c_g': G_costs, 'c_s': S_costs, 'v': V_costs}

def pass_to_sim(model):
	return run_sim(model, S_init, I_init)

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
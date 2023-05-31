import itertools

import numpy as np
import pickle as pkl

def get_trans(model, S, I):
	'''
	Calculate the transitivity slope for a particular set of equilibrium 
	conditions

	Args:
		S: Vector of length 4 of equilibrium host abundances
		I: Vector of length 3 of equilibrium pathogen abundances
		**kwargs: Set of parameters defining a simulation

	Returns:
		full_sib_end: vector of full sibling familiy susceptibilities for the 
			endemic pathogen
		half_sib_freq: expected frequency of each half-sib family
	'''

	#Get the proability of each mating pair
	freq_matrix = np.outer(S, S)
	freq_matrix = freq_matrix / (np.sum(S))**2

	#Calculate the susceptibility of each genotype to the endemic pathogen 
	#by weighting their suseptibility to vir and avr by the relative abundances
	I_freq = I / np.sum(I)
	sus = np.dot(model.B, I_freq)

	full_sib = np.zeros((model.S_genotypes**2, 3))
	half_sib = np.zeros((model.S_genotypes, 3))

	full_sib[:, 0] = np.dot(model.M, sus)			#Full sibling susceptibility to endemic pathogen
	full_sib[:, 1] = np.dot(model.M, model.B[:, 2])	#Full sib susceptiblity to foreign pathogen
	full_sib[:, 2] = freq_matrix.flatten()			#Full sib family frequency

	half_sib[:, 0] = np.average(full_sib[:, 0].reshape(-1, model.S_genotypes), axis=1)
	half_sib[:, 1] = np.average(full_sib[:, 1].reshape(-1, model.S_genotypes), axis=1)
	half_sib[:, 2] = np.sum(full_sib[:, 2].reshape(-1, model.S_genotypes), axis=1)

	return (full_sib, half_sib)

def load_data(path, var_1, var_2):
	#Load simulation raster data
	with open(path, 'rb') as f:
		models, data = pkl.load(f)
	
	#Get all parameter values for focal parameters
	x_vals = np.sort(list(set([getattr(model, var_1) for model in models])))
	y_vals = np.sort(list(set([getattr(model, var_2) for model in models])))

	#Number of unique focal parameter values
	n_x = len(x_vals)
	n_y = len(y_vals)   

	allele_freq_raster = np.zeros((models[0].n_loci, n_x, n_y))
	V_raster, T_raster = [np.zeros((n_x, n_y)) for _ in range(2)]

	for i in range(len(data)):
		x_param = getattr(models[i], var_1)
		y_param = getattr(models[i], var_2)
		x_ind = np.where(x_vals == x_param)
		y_ind = np.where(y_vals == y_param) 

		t, sus, inf = data[i]
		N = np.sum(sus, axis=0)

		#Find the t indexes equal to the range over which to average solution
		avg_range = (t[-1] - 2000, t[-1])
		ind_1 = np.argmin(abs(t - avg_range[0]))
		ind_2 = np.argmin(abs(t - avg_range[1]))

		sus_freq = np.average((sus/N)[:, ind_1:ind_2], axis=1)
		avir_freq = np.average(inf[0, ind_1:ind_2] / np.sum(inf[:, ind_1:ind_2], axis=0))

		#Compute the average allele frequency for each locus
		for j in range(models[i].n_loci):
			allele_freq_raster[j, x_ind, y_ind] = np.sum(sus_freq[models[i].G[:, j] == 1])

		#Take the average solution value over range to account for oscillations
		S_avg = np.average(sus[:, ind_1:ind_2], axis=1)
		I_avg = np.average(inf[:, ind_1:ind_2], axis=1)
		full_sib, _ = get_trans(models[i], S_avg, I_avg)

		#Check if there is sufficient G polymorphism, and if so, append the transitivity slope
		if 1e-2 < allele_freq_raster[1, x_ind, y_ind] < 1 - 1e-2:
			T_raster[x_ind, y_ind] = np.polyfit(full_sib[:, 0], full_sib[:, 1], 1, w=full_sib[:, 2])[0]
		#Set the transitivity slope to 0 if there isn't sufficient polymorpism
		else:
			T_raster[x_ind, y_ind] = 0
	
		V_raster[x_ind, y_ind] = avir_freq

	return allele_freq_raster, V_raster, T_raster
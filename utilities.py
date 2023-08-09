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

	full_sib[:, 0] = np.dot(model.M, sus)				#Full sibling susceptibility to endemic pathogen
	full_sib[:, 1] = np.dot(model.M, model.B[:, 2])		#Full sib susceptiblity to foreign pathogen
	full_sib[:, 2] = freq_matrix.flatten()				#Full sib family frequency

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
	V_raster, T_raster, D_raster = [np.zeros((n_x, n_y)) for _ in range(3)]

	for i in range(len(data)):
		x_param = getattr(models[i], var_1)
		y_param = getattr(models[i], var_2)
		x_ind = np.where(x_vals == x_param)
		y_ind = np.where(y_vals == y_param) 

		sus, inf, _ = data[i]

		#Find the t indexes equal to the range over which to average solution
		N = np.sum(sus)
		sus_freq = sus/N
		avir_freq = inf[0] / np.sum(inf)

		#Compute the average allele frequency for each locus
		for j in range(models[i].n_loci):
			allele_freq_raster[j, x_ind, y_ind] = np.sum(sus_freq[models[i].G[:, j] == 1])

		G_freq = np.sum(sus[models[i].G[:,1] == 1]) / N
		S_freq = np.sum(sus[models[i].G[:,2] == 1]) / N
		GS_freq = np.sum(sus[np.logical_and(models[i].G[:,1], models[i].G[:,2] == 1)]) / N

		D = GS_freq - G_freq*S_freq

		D_max = np.max((-1*G_freq*S_freq, -1*(1-G_freq)*(1-S_freq)))
		D_min = np.min((G_freq*(1-S_freq), (1-G_freq)*S_freq))
		min_freq = np.min((G_freq, S_freq, 1-G_freq, 1-S_freq))

		if D < 0 and min_freq > 0.01:
			D_prime = D / D_max
		elif D > 0 and min_freq > 0.01:
			D_prime = D / D_min
		else:
			D_prime = 0

		#Take the average solution value over range to account for oscillations
		full_sib, _ = get_trans(models[i], sus, inf)

		#Check if there is sufficient G polymorphism, and if so, append the transitivity slope
		if 1e-2 < allele_freq_raster[1, x_ind, y_ind] < 1 - 1e-2:
			T_raster[x_ind, y_ind] = np.polyfit(full_sib[:, 0], full_sib[:, 1], 1, w=full_sib[:, 2])[0]
		#Set the transitivity slope to 0 if there isn't sufficient polymorpism
		else:
			T_raster[x_ind, y_ind] = 0
	
		V_raster[x_ind, y_ind] = avir_freq
		D_raster[x_ind, y_ind] = D_prime

	return allele_freq_raster, V_raster, T_raster, D_raster

def check_stab(path, var_1, var_2):
	#Load simulation raster data
	with open(path, 'rb') as f:
		models, data = pkl.load(f)
	
	#Get all parameter values for focal parameters
	x_vals = np.sort(list(set([getattr(model, var_1) for model in models])))
	y_vals = np.sort(list(set([getattr(model, var_2) for model in models])))

	#Number of unique focal parameter values
	n_x = len(x_vals)
	n_y = len(y_vals)   

	stab = np.zeros((n_x, n_y))

	for i in range(len(data)):
		x_param = getattr(models[i], var_1)
		y_param = getattr(models[i], var_2)
		x_ind = np.where(x_vals == x_param)
		y_ind = np.where(y_vals == y_param) 

		_, _, eigs = data[i]

		#Check if all eigenvalues are negative (except last which corresponds to foreign pathogen)
		if np.prod(np.real(eigs[:-1]) < 0):
			stab[x_ind, y_ind] = 1
		else:
			print(np.real(eigs))
			stab[x_ind, y_ind] = 0

	return stab

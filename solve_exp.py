import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root, approx_fprime

from model import Model

def run_sim(model, S_0, I_0, t=(0,5000)):
	'''
	Run ODE simulation for three locus, three pathogen genotype model

	Args:
		model: Model class instance
		af_S: Initial allele frequencies for each of the three host allele
		af_I: Initial frequency of the Avir pathogen genotype

	Returns:
		sol.t: Time points corresponding to the solution
		S: Solution for susceptible host abundances [genotype, time]
		I: Solution for infected host abundances [genotype, time]
	'''
	def df(t, X):
		#Seperate out uninfected and infected hosts
		S = X[:model.S_genotypes]
		I = X[model.S_genotypes:]

		N = np.sum(S) + np.sum(I)
 
		#Get the frequency of each genotype
		genotype_freq = S / np.sum(S)

		#Get parental pair frequencies and adjust by fecundity costs
		pair_freq = np.outer(model.C*genotype_freq, genotype_freq).flatten()

		dS = np.sum(S)*np.dot(pair_freq, model.M) - \
			S*(model.k*N + model.mu + np.dot(model.B, I)/N)
		dI = I*(np.dot(model.B.T, S)/N - model.mu)

		X_out = np.append(dS, dI)

		return X_out

	X_0 = np.append(S_0, I_0)
	sol = solve_ivp(df, t, X_0, method='DOP853', max_step=0.5)

	S = sol.y[:model.S_genotypes, :]
	I = sol.y[model.S_genotypes:, :]

	return sol.t, S, I
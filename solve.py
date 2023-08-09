import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root, approx_fprime

from model import Model

def get_sol(model, af_S, af_I, t=(0,5000), init_hosts=400, init_inf=10):
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

	def pass_to_df(X):
		return df(None, X)

	#Assign host genotype ICs based on allele frequencies
	S_0 = np.ones(model.S_genotypes)		
	for i in range(model.n_loci):
		S_0[model.G[:,i] == 0] = S_0[model.G[:,i] == 0] * (1 - af_S[i])
		S_0[model.G[:,i] == 1] = S_0[model.G[:,i] == 1] * (af_S[i])
	
	#Assign infected ICs based on Avr frequency
	I_0 = np.zeros(3)
	I_0[0] = af_I
	I_0[1] = 1 - af_I

	X_0 = np.append(S_0 * init_hosts, I_0 * init_inf)
	sol = solve_ivp(df, t, X_0, method='DOP853')
	eq = root(pass_to_df, sol.y[:,-1], tol=1e-10)

	S = eq.x[:model.S_genotypes]
	I = eq.x[model.S_genotypes:]
	eigs = np.linalg.eig(approx_fprime(eq.x, pass_to_df))[0]
	
	if not eq.success:
		print(eq.message)

	return S, I, eigs
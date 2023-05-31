import numpy as np
from scipy.integrate import solve_ivp

from model import Model

def run_sim(model, af_S, af_I, t=(0,10000), init_hosts=400, init_inf=10):

	def df(t, X):
		S = X[:model.S_genotypes]
		I = X[model.S_genotypes:]

		N = np.sum(S) + np.sum(I)
 
		genotype_freq = S / np.sum(S)
		pair_freq = np.outer(model.F*genotype_freq, genotype_freq).flatten()

		dS = np.sum(S)*np.dot(pair_freq, model.M) - \
			S*(model.k*N + model.mu + np.dot(model.B, I)/N)
		dI = I*(np.dot(model.B.T, S)/N - model.mu)

		X_out = np.append(dS, dI)

		return X_out

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

	S = sol.y[:model.S_genotypes, :]
	I = sol.y[model.S_genotypes:, :]

	return sol.t, S, I
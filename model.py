import itertools

import numpy as np

class Model:

	def __init__(self, **kwargs):
		self.n_loci = 3
		self.I_genotypes = 3
		self.S_genotypes = 2**self.n_loci

		self.G = np.array(list(itertools.product([0, 1], repeat=self.n_loci)))

		self.rho = [0.05, 0]

		self.b = 1
		self.mu = 0.2
		self.k = 0.001

		self.beta = 0.5
		self.g = 0.3
		self.s = 0.9
		self.v = 0.2

		self.c_g = 0.1
		self.c_s = 0.2

		for key, value in kwargs.items():
			setattr(self, key, value)

		costs = [0, self.c_g, self.c_s]

		self.F = 1 - np.dot(self.G, costs)
		self.B = self.transmission_matrix()
		self.M = self.mating_matrix()

	def transmission_matrix(self):
		B = np.ones((self.S_genotypes, self.I_genotypes)) * self.beta

		for i, host in enumerate(self.G):
			if host[1] == 1:
				B[i,:] = B[i,:] * (1-self.g)
			
			if host[2] == 1:
				B[i,0] = B[i,0] * (1-self.s)
			else:
				B[i,1] = B[i,1] * (1-self.v)

		return B
				
	def mating_matrix(self):
		paths = np.array(list(itertools.product([0, 1], repeat=self.n_loci)))
		p_paths = np.zeros((self.S_genotypes, 2))

		'''
		for i, path in enumerate(paths):
			p_locus = np.zeros(self.n_loci)
			p_locus[0] = 0.5

			for j in range(self.n_loci - 1):
				if path[j] == path[j+1]:
					p_locus[j+1] = 1 - self.rho if np.isscalar(self.rho) else 1 - self.rho[j]
				else:
					p_locus[j+1] = self.rho if np.isscalar(self.rho) else self.rho[j]
			
			p_paths[i] = np.product(p_locus)
		'''
		for i, path in enumerate(paths):
			p_locus = np.zeros((self.n_loci, 2))
			p_locus[0:2] = 0.5

			if path[2] == path[1]:
				p_locus[2, 0] = 1 - self.rho[0]
				p_locus[2, 1] = 1 - self.rho[1]
			else:
				p_locus[2, 0] = self.rho[0]
				p_locus[2, 1] = self.rho[1]
			
			p_paths[i] = np.product(p_locus, axis=0)

		M = np.zeros((self.S_genotypes, self.S_genotypes, self.S_genotypes))
		for i, parental in enumerate(self.G):
			for j, maternal in enumerate(self.G):
				offspring = np.zeros(paths.shape)
				offspring[paths==0] = np.tile(parental, (self.S_genotypes, 1))[paths==0]
				offspring[paths==1] = np.tile(maternal, (self.S_genotypes, 1))[paths==1]
				
				for k, progeny in enumerate(self.G):
					matches = np.product(offspring == progeny, axis=1, dtype=bool)
					
					if maternal[0] == 0:
						M[i,j,k]= np.dot(matches, p_paths[:,0])
					else:
						M[i,j,k]= np.dot(matches, p_paths[:,1])
		
		M = M.reshape((self.S_genotypes**2, self.S_genotypes))
		return M

def collapse_locus(model, S, locus):
	'''
	Reduce the number of genotypes in a simulation by undifferentiating a particular
	locus

	Args:
		model: Model class instance
		S: Simulation results for uninfected hosts
		locus: which locus to colapse

	Returns:
		S_prime: Transformed simulation results
		G_prime: Matrix of transformed genotypes

	'''
	G_prime = np.array(list(itertools.product([0, 1], repeat=model.n_loci-1)))
	S_prime = np.zeros((int(S.shape[0] / 2), S.shape[1]))

	for i, gtp in enumerate(G_prime):
		index = np.product(np.delete(model.G, locus, 1) == gtp, axis=1, dtype=bool)
		S_prime[i,:] = np.sum(S[index, :], axis=0)

	return S_prime, G_prime
import itertools
import numpy as np

class Model:

	def __init__(self, **kwargs):
		self.n_loci = 3	#Number of host loci, fixed at 3
		self.S_genotypes = 2**self.n_loci	#Number of host genotypes
		self.I_genotypes = 3	#Number of pathogen genotypes, fixed at 3

		#Matrix of all possible genotypes, G[i,j] is the jth allele of the ith genotype
		self.G = np.array(list(itertools.product([0, 1], repeat=self.n_loci)))

		#Recombination rate vector, rho[0] is the recombination rate for allele 0, rho[1] is the 
		#recombination rate for allele 1
		self.rho = [0.05, 0]

		self.b = 1				#Default birthrate
		self.mu = 0.2			#Default deathrate
		self.k = 0.001			#Coefficient of density-dependent growth

		self.beta = 0.5			#Baseline transmission rate
		self.g = 0.3			#Default strength of general resistance
		self.s = 0.9			#Default strength of specific resistance
		self.v = 0.2			#Default cost of virulence
		self.nh = 0				#Nonhost resistance

		self.c_g = 0.1			#Default cost of general resistance
		self.c_s = 0.2			#Default cost of specific resistance

		self.sel = 'hard'		#Form of virlence costs (hard vs soft)

		for key, value in kwargs.items():
			setattr(self, key, value)

		self.C = self.cost_vector()
		self.B = self.transmission_matrix()
		self.M = self.mating_matrix()

	def cost_vector(self):
		'''
		Define the cost vector, where C[i] is the total fecundity cost for genotype i,
		assuming multiplicitive interactions between loci

		Returns:
			C: Cost vector
		'''

		allele_costs = [0, self.c_g, self.c_s]

		costs = self.G * allele_costs

		C = np.ones(self.G.shape[0])
		for i in range(self.G.shape[1]):
			C = C * (1 - costs[:,i])

		return C

	def transmission_matrix(self):
		'''
		Define the transmission matrix B, where B[i,j] is the transmission rate of the
		jth pathogen genotype on the ith host genotype

		Returns:
			B: Transmission matrix
		'''
		
		#Define transmission matrix and apply nonhost resistance
		B = np.ones((self.S_genotypes, self.I_genotypes)) * self.beta
		B[:,2] = B[:,2]*(1-self.nh)

		#Soft selection
		if self.sel == 'soft':
			for i, host in enumerate(self.G):
				#Apply general resistance to all pathogen types
				if host[1] == 1:
					B[i,:] = B[i,:] * (1-self.g)
				
				#Apply specific resistance for resistant genotype
				if host[2] == 1:
					B[i,0] = B[i,0] * (1-self.s)

				#Apply virulence costs for susceptible genotype
				else:
					B[i,1] = B[i,1] * (1-self.v)

		#Hard selection
		if self.sel == 'hard':
			B[:,1] = B[:,1] * (1-self.v)

			for i, host in enumerate(self.G):
				#Apply general resistance to all pathogen types
				if host[1] == 1:
					B[i,:] = B[i,:] * (1-self.g)
				
				#Apply specific resistance for resistant genotype
				if host[2] == 1:
					B[i,0] = B[i,0] * (1-self.s)
					
		return B
				
	def mating_matrix(self):
		'''
		Define the matring matrix M, where M[i,j] is the probability of offspring genotype j given
		parental combination i (reduced from third order tensor)

		Returns:
			M: Mating matrix
		'''

		paths = np.array(list(itertools.product([0, 1], repeat=self.n_loci)))
		p_paths = np.zeros((self.S_genotypes, 2))

		for i, path in enumerate(paths):
			p_locus = np.zeros((self.n_loci, 2))
			p_locus[0:2] = 0.5

			if path[2] == path[1]:
				p_locus[2, 0] = 1 - self.rho[0]
				p_locus[2, 1] = 1 - self.rho[1]
			else:
				p_locus[2, 0] = self.rho[0]
				p_locus[2, 1] = self.rho[1]
			
			p_paths[i] = np.prod(p_locus, axis=0)

		M = np.zeros((self.S_genotypes, self.S_genotypes, self.S_genotypes))
		for i, parental in enumerate(self.G):
			for j, maternal in enumerate(self.G):
				offspring = np.zeros(paths.shape)
				offspring[paths==0] = np.tile(parental, (self.S_genotypes, 1))[paths==0]
				offspring[paths==1] = np.tile(maternal, (self.S_genotypes, 1))[paths==1]
				
				for k, progeny in enumerate(self.G):
					matches = np.prod(offspring == progeny, axis=1, dtype=bool)
					
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
		index = np.prod(np.delete(model.G, locus, 1) == gtp, axis=1, dtype=bool)
		S_prime[i,:] = np.sum(S[index, :], axis=0)

	return S_prime, G_prime
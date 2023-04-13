import numpy as np
import tqdm
import pickle as pkl

import multiprocessing as mp

from allelic_model import run_sim

#Default paramaters, some may be overwritten for individual simulations
def_params = {  'k':0.001,      #Coefficient of density dependent growth
                'mu':0.2,       #Deathrate
                'b':1,          #Birthrate
                'beta':0.5,     #Baseline transmission rate
                'nh':0.1,       #Nonhost transmission rate
                'g':0.3,        #Strength of general resistance
                's':0.9,        #Strength of specific resistance
                'rho':0.05,     #Recombination rate
                'c_g':0.1,     #Cost of general resistance
                'c_s':0.2,      #Cost of specific resistance
                'v':0.2,        #Cost of virulence
                'S_0':[100,100,100,100],
                'I_0':[10,1,0]}        

size = 100

output_path = 'Data/Mult_Costs/cov_gv_costs.p'

V_costs = np.linspace(0, 0.3, size) #Range of parameters for virulence costs
G_costs = np.linspace(0, 0.2, size) #Range of parameters for general resistance costs
S_costs = np.linspace(0, 0.4, size) #Range of parameters for speific resistance costs

def pass_to_sim(kwargs):
	return run_sim(**kwargs)

if __name__ == '__main__':
    coords = []
    params = []

    for i in range(size):
        for j in range(size):
            coords.append((i,j))

            new_param = def_params.copy()            
            new_param['c_g'] = G_costs[i]
            #new_param['c_s'] = S_costs[j]
            new_param['v'] = V_costs[j]

            params.append(new_param)
		
    #Run simluations for 4 core processor
    pool = mp.Pool(processes=4)	
    
    results = []
    for result in tqdm.tqdm(pool.imap(pass_to_sim, params), total=len(params)):
        results.append(result)

    raster = []
    for i in range(size):
        inds = [j for j in range(len(coords)) if coords[j][1] == i]
        raster.append([results[j] for j in inds])

    with open(output_path, 'wb') as f:
        pkl.dump([results, params], f)
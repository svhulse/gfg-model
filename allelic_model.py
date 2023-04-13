import numpy as np
from scipy.integrate import solve_ivp

from utilities import transmission_matrix, mating_matrix

def run_sim(**kwargs):
    params = {  'k':0.001,                  #Coefficient of density dependent growth
                'mu':0.2,                   #Deathrate
                'b':1,                      #Birthrate
                'beta':0.5,                 #Baseline transmission rate
                'nh':0.1,                   #Nonhost transmission rate
                'g':0.3,                    #Strength of general resistance
                's':0.9,                    #Strength of specific resistance
                'rho':0.05,                 #Recombination rate
                'c_g':0.1,                 #Cost of general resistance
                'c_s':0.2,                  #Cost of specific resistance
                'v':0.2,                    #Cost of virulence
                'S_0':[100,100,100,100],    #Initial S densities
                'I_0':[10,1,0]}             #Initial I densities

    params.update(kwargs)

    k = params['k']
    mu = params['mu']
    b = params['b']
    c_g = params['c_g']
    c_s = params['c_s']

    S_0 = np.array(params['S_0']) #Load initial conditions for hosts
    I_0 = np.array(params['I_0']) #Load initial conditions for infected hosts

    costs = np.zeros(4)
    costs[0] = 1-(1-c_g)*(1-c_s)    #Costs for GS
    costs[1] = c_g          #Costs for Gs
    costs[2] = c_s          #Costs for gS
    costs[3] = 0            #Costs for gs

    b = b - costs #Subtract costs from the baseline birthrate

    B = transmission_matrix(**params)
    M = mating_matrix(**params)

    def df(t, X):
        S = X[0:4]
        I = X[4:7]

        N = np.sum(S) + np.sum(I)
 
        genotype_freq = S / np.sum(S)
        pair_freq = np.outer(b*genotype_freq, genotype_freq).flatten()

        dS = np.sum(S)*np.dot(pair_freq, M) - S*(k*N + mu + np.dot(B, I)/N)
        dI = I*(np.dot(B.T, S)/N - mu)

        X_out = np.concatenate((dS, dI))

        return X_out
    
    t = (0, 10000)

    X_0 = np.concatenate((S_0, I_0))
    sol = solve_ivp(df, t, X_0, method='DOP853')

    S = sol.y[0:4, :]
    I = sol.y[4:7, :]

    return sol.t, S, I
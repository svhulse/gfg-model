import numpy as np
import pickle as pkl

def transmission_matrix(**kwargs):
    '''
	Calculate the transmission matrix based on a simulations parameters

	Args:
		**kwargs: Set of parameters defining a simulation, must include g, s, v, 
        nh, and beta

	Returns:
		B: Matrix defining the transmission rates between each host and each 
        pathogen 
	'''

    g = kwargs['g']
    s = kwargs['s']
    v = kwargs['v']
    nh = kwargs['nh']
    beta = kwargs['beta']

    #Matrix of transmission rates
    B = np.ones((4, 3))

    #QR transmission rates
    B[0,0] = (1-g)*(1-s)
    B[0,1] = 1-g
    B[0,2] = (1-g)*(1-nh)

    #Qr transmission rates
    B[1,0] = 1-g
    B[1,1] = (1-g)*(1-v)
    B[1,2] = (1-g)*(1-nh)

    #qR transmission rates
    B[2,0] = 1-s
    B[2,1] = 1
    B[2,2] = 1-nh

    #qr transmission rates
    B[3,0] = 1
    B[3,1] = 1-v
    B[3,2] = 1-nh

    B = B*beta

    return B

def mating_matrix(**kwargs):
    '''
	Calculate the mating matrix for a particular recombination rate

	Args:
		**kwargs: Set of parameters defining a simulation, must include rho 

	Returns:
		M: Mating matrix determining how offspring are distributed
	'''

    rho = kwargs['rho']

    M = np.zeros((4, 16))
    M[:, 0] = [1, 0, 0, 0]                                       #GS x GS
    M[:, 1] = [0.5, 0.5, 0, 0]                                   #GS x Gs
    M[:, 2] = [0.5, 0, 0.5, 0]                                   #GS x gS
    M[:, 3] = [0.5*(1-rho), 0.5*rho, 0.5*rho, 0.5*(1-rho)]       #GS x gs
    M[:, 4] = [0.5, 0.5, 0, 0]                                   #Gs x GS
    M[:, 5] = [0, 1, 0, 0]                                       #Gs x Gs
    M[:, 6] = [0.5*rho, 0.5*(1-rho), 0.5*(1-rho), 0.5*rho]       #Gs x gS
    M[:, 7] = [0, 0.5, 0, 0.5]                                   #Gs x gs
    M[:, 8] = [0.5, 0, 0.5, 0]                                   #gS x GS
    M[:, 9] = [0.5*rho, 0.5*(1-rho), 0.5*(1-rho), 0.5*rho]       #gS x Gs
    M[:, 10] = [0, 0, 1, 0]                                      #gS x gS
    M[:, 11] = [0, 0, 0.5, 0.5]                                  #gS x gs
    M[:, 12] = [0.5*(1-rho), 0.5*rho, 0.5*rho, 0.5*(1-rho)]      #gs X GS
    M[:, 13] = [0, 0.5, 0, 0.5]                                  #gs x Gs
    M[:, 14] = [0, 0, 0.5, 0.5]                                  #gs x gS
    M[:, 15] = [0, 0, 0, 1]                                      #gs x gs

    return M

def get_trans(S, I, **kwargs):
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
        full_sib_for: vector of full sibling famility susceptibilities for the
            foreign pathogen
        full_sib_freq: expected frequency of each full sibling family
        half_sib_end: vector of half sibling familiy susceptibilities for the 
            endemic pathogen
        half_sib_for: vector of half sibling familiy susceptibilities for the 
            foreign pathogen
        half_sib_freq: expected frequency of each half-sib family
	'''

    B = transmission_matrix(**kwargs)
    M = mating_matrix(**kwargs)

    #Get the proability of each mating pair
    freq_matrix = np.outer(S, S)
    freq_matrix = freq_matrix / (np.sum(S))**2

    #Calculate the susceptibility of each genotype to the endemic pathogen 
    #by weighting their suseptibility to vir and avr by the relative abundances
    I_freq = I / np.sum(I)
    sus = np.dot(B, I_freq)

    full_sib_end = np.dot(M.T, sus)
    full_sib_for = np.dot(M.T, B[:, 2])
    full_sib_freq = freq_matrix.flatten()

    half_sib_end = np.average(full_sib_end.reshape(-1, 4), axis=1)
    half_sib_for = np.average(full_sib_for.reshape(-1, 4), axis=1)
    half_sib_freq = np.sum(full_sib_freq.reshape(-1, 4), axis=1)

    return (full_sib_end, full_sib_for, full_sib_freq, half_sib_end, half_sib_for, half_sib_freq)

def load_data(path, var_1, var_2):
    '''
	Load raw pickle data for a simulation raster set

	Args:
		path: path of the pickle file to be loaded
        var_1: first simulation parameter that is varied in the raster
        var_2: second simulation parameter that is varied in the raster

	Returns:
		G_raster: equilibrum G allele frequencies
        S_raster: equilibrum S allele frequencies
        V_raster: proportion of virulent endemic pathogens
        T_raster: transitivity slope
        D_raster: linkage disequilibrum (D')
	'''

    #Load simulation raster data
    with open(path, 'rb') as f:
        data, params = pkl.load(f)
    
    #Get all parameter values for focal parameters
    x_vals = np.sort(list(set([param[var_1] for param in params])))
    y_vals = np.sort(list(set([param[var_2] for param in params])))

    #Number of unique focal parameter values
    n_x = len(x_vals)
    n_y = len(y_vals)   

    G_raster, S_raster, V_raster, T_raster, D_raster = [np.zeros((n_x, n_y)) for _ in range(5)]

    for i in range(len(data)):
        x_param = params[i][var_1]
        y_param = params[i][var_2]

        x_ind = np.where(x_vals == x_param)
        y_ind = np.where(y_vals == y_param) 

        t, sus, inf = data[i]
        N = np.sum(sus, axis=0)

        #Find the t indexes equal to the range over which to average solution
        avg_range = (6000, 8000)
        ind_1 = np.argmin(abs(t - avg_range[0]))
        ind_2 = np.argmin(abs(t - avg_range[1]))

        sus_freq = np.average((sus/N)[:, ind_1:ind_2], axis=1)
        avir_freq = np.average(inf[0, ind_1:ind_2] / np.sum(inf[:, ind_1:ind_2], axis=0))
                
        G_freq = sus_freq[0] + sus_freq[1]
        S_freq = sus_freq[0] + sus_freq[2]
        GS_freq = sus_freq[0]

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
        S_avg = np.average(sus[:, ind_1:ind_2], axis=1)
        I_avg = np.average(inf[:, ind_1:ind_2], axis=1)
        full_sib_end, full_sib_for, full_sib_freq, _, _, _ = get_trans(S_avg, I_avg, **params[i])

        if sorted(full_sib_freq)[-2] > 1e-2:
            T_raster[x_ind, y_ind] = np.polyfit(full_sib_end, full_sib_for, 1, w=full_sib_freq)[0]
        else:
            T_raster[x_ind, y_ind] = 0
    
        G_raster[x_ind, y_ind] = G_freq
        S_raster[x_ind, y_ind] = S_freq
        D_raster[x_ind, y_ind] = D_prime
        V_raster[x_ind, y_ind] = avir_freq

    return G_raster, S_raster, V_raster, T_raster, D_raster

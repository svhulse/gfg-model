import numpy as np
import os

from model import Model
from run_sim import run_sim
from utilities import get_trans, collapse_locus

import matplotlib
from matplotlib import pyplot as plt

SMALL_SIZE = 8
MEDIUM_SIZE = 10
BIGGER_SIZE = 14

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=SMALL_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

#plt.rcParams["font.family"] = "Helvetica Neue"
matplotlib.rcParams['figure.dpi'] = 300

params = {      'k':0.001, 
                'mu':0.2,  
                'b':1,     
                'beta':0.5,
                'nh':0.1,
                'g':0.3,
                's':0.9,
                'rho':[0.05, 0.0], 
                'c_g':0.1,     
                'c_s':0.2,      
                'v':0.2}

sim = Model(**params)

t, S, I = run_sim(sim, [0.1, 0.5, 0.5], 0.9, t=(0, 100000))

Sc, _ = collapse_locus(sim, S, 0)

labels = [r'$\it{gs}$', r'$\it{gS}$', r'$\it{Gs}$', r'$\it{GS}$']

fig, ax = plt.subplots(nrows = 1, ncols = 3, figsize=(8,3))
fig.tight_layout()
plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.25, hspace=0.25)

f_R = np.sum(S[sim.G[:,0] == 1], axis=0) / np.sum(S, axis=0)
f_G = np.sum(S[sim.G[:,1] == 1], axis=0) / np.sum(S, axis=0)
f_S = np.sum(S[sim.G[:,2] == 1], axis=0) / np.sum(S, axis=0)

ax[0].plot(t, f_R, label = 'Linkage Modifier')
ax[0].plot(t, f_G, label = 'General Resistance')
ax[0].plot(t, f_S, label = 'Specific Resistance')

ax[0].legend(loc='upper left')
ax[0].set_xlabel('Time')
ax[0].set_ylabel('Frequency')
ax[0].set_title('Host Genotypes', fontsize=BIGGER_SIZE, pad=25)
ax[0].ticklabel_format(axis='x', style='sci', scilimits=(0,0))

ax[1].plot(t, I[0,:]/(I[0,:] + I[1,:]), label=r'$\it{Avr}$')
ax[1].plot(t, I[1,:]/(I[0,:] + I[1,:]), label=r'$\it{Vir}$')

ax[1].set_ylim([0, 1])
ax[1].set_xlabel('Time')
ax[1].set_ylabel('Frequency')
ax[1].legend(loc='upper left')
ax[1].set_title('Pathogen Genotypes', fontsize=BIGGER_SIZE, pad=25)
ax[1].ticklabel_format(axis='x', style='sci', scilimits=(0,0))

trans_slope = np.zeros(len(t))

for i in range(len(t)):
    trans = get_trans(sim, S[:,i], I[:,i])
    trans_slope[i], _ = np.polyfit(trans[3], trans[4], 1, w=trans[5])

ax[2].plot(t, trans_slope)
ax[2].set_title('Resistance Transitivity', fontsize=BIGGER_SIZE, pad=25)
ax[2].ticklabel_format(axis='x', style='sci', scilimits=(0,0))

ax[0].annotate("A", xy=(-0.2, 1.05), xycoords="axes fraction", fontsize=BIGGER_SIZE)
ax[1].annotate("B", xy=(-0.2, 1.05), xycoords="axes fraction", fontsize=BIGGER_SIZE)
ax[2].annotate("C", xy=(-0.2, 1.05), xycoords="axes fraction", fontsize=BIGGER_SIZE)

fig.savefig('fig_5.svg', bbox_inches='tight', pad_inches=0.1)
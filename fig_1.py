import numpy as np

from model import Model
from solve_exp import run_sim
from utilities import get_trans
from model import collapse_locus

import style
from matplotlib import pyplot as plt

params = {      'k':0.001, 
                'mu':0.2,  
                'b':1,     
                'beta':0.5,
                'nh':0.1,
                'g':0.3,
                's':0.9,
                'rho':[0.05, 0.05], 
                'c_g':0.1,     
                'c_s':0.2,      
                'v':0.2,
                'sel':'soft'}

sim = Model(**params)

init_cond = [0, 0.1, 0.1]
S_0 = np.ones(sim.S_genotypes)		

for i in range(sim.n_loci):
	S_0[sim.G[:,i] == 0] = S_0[sim.G[:,i] == 0] * (1 - init_cond[i])
	S_0[sim.G[:,i] == 1] = S_0[sim.G[:,i] == 1] * (init_cond[i])

t, S, I = run_sim(sim, S_0, [0.9, 0.1, 0], t=(0,2000))

Sc, _ = collapse_locus(sim, S, 0)

trans = get_trans(sim, S[:,-1], I[:,-1])

host_labels = [r'$\it{gs}$', r'$\it{gS}$', r'$\it{Gs}$', r'$\it{GS}$']
path_labels = [r'$\it{Avr}$', r'$\it{vir}$']

fig, ax = plt.subplots(nrows = 1, ncols = 3, figsize=(8,3))
fig.tight_layout()
plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.25, hspace=0.25)

ax[0].plot(t[75:], Sc.T[75:,:], label = host_labels)

ax[0].legend(loc='upper right')
ax[0].set_xlabel('Time')
ax[0].set_ylabel('Abundance')

ax[1].plot(t[75:], I.T[75:,0:2], label = path_labels)

ax[1].legend(loc='upper right')
ax[1].set_xlabel('Time')
ax[1].set_ylabel('Abundance')

m_f, b_f = np.polyfit(trans[0][:,0], trans[0][:,1], 1, w=trans[0][:,2])

ax[2].scatter(trans[0][:,0], trans[0][:,1], s = trans[0][:,2]*1000, label = 'Full-sib', c='C1')
ax[2].plot(trans[0][:,0], b_f + m_f*trans[0][:,0])

ax[2].set_xlim([0, 0.6])
ax[2].set_ylim([0, 0.6])
ax[2].set_xlabel('Mean Susceptibility to Endemic')
ax[2].set_ylabel('Susceptibility to Foreign')

ax[0].annotate("A", xy=(-0.2, 1.05), xycoords="axes fraction", fontsize=style.bigger)
ax[1].annotate("B", xy=(-0.2, 1.05), xycoords="axes fraction", fontsize=style.bigger)
ax[2].annotate("C", xy=(-0.2, 1.05), xycoords="axes fraction", fontsize=style.bigger)

fig.savefig('./figures/fig_1.svg', bbox_inches='tight', pad_inches=0.1)
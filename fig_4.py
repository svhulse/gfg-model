import numpy as np

from model import Model
from solve_exp import run_sim
from utilities import get_trans, load_data
from model import collapse_locus

from matplotlib import pyplot as plt
from matplotlib import ticker

import style
from style import label_axes

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
                'v':0.2,
                'sel':'soft'}

sim = Model(**params)

af_S1 = [0, 0.1, 0.1]
af_S2 = [1, 0.1, 0.1]

#Assign host genotype ICs based on allele frequencies
S_0_recomb = np.ones(sim.S_genotypes)		
S_0_nonrecomb = np.ones(sim.S_genotypes)		

for i in range(sim.n_loci):
	S_0_recomb[sim.G[:,i] == 0] = S_0_recomb[sim.G[:,i] == 0] * (1 - af_S1[i])
	S_0_recomb[sim.G[:,i] == 1] = S_0_recomb[sim.G[:,i] == 1] * (af_S1[i])
	S_0_nonrecomb[sim.G[:,i] == 0] = S_0_nonrecomb[sim.G[:,i] == 0] * (1 - af_S2[i])
	S_0_nonrecomb[sim.G[:,i] == 1] = S_0_nonrecomb[sim.G[:,i] == 1] * (af_S2[i])

t1, S_1, I_1 = run_sim(sim, S_0_recomb, [0.9, 0.1, 0], t=(0,20000))
t2, S_2, I_2 = run_sim(sim, S_0_nonrecomb, [0.9, 0.1, 0], t=(0,20000))

S_1c, _ = collapse_locus(sim, S_1, 0)
S_2c, _ = collapse_locus(sim, S_2, 0)

trans_1 = get_trans(sim, S_1[:,-1], I_1[:,-1])
trans_2 = get_trans(sim, S_2[:,-1], I_2[:,-1])

def_path = './data/cov_gs.p'
fxt_path = './data/cov_gs_rho0.p'

_, _, def_slope, _ = load_data(def_path, 'c_s', 'c_g') 
_, _, fxt_slope, _ = load_data(fxt_path, 'c_s', 'c_g')

labels = [r'$\it{gs}$', r'$\it{gS}$', r'$\it{Gs}$', r'$\it{GS}$']

fig, ax = plt.subplots(nrows = 2, ncols = 3, figsize=(10,6))
fig.tight_layout()
plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.25, hspace=0.25)

ax[0,0].plot(t1, S_1c.T, label = labels)

ax[0,0].legend()
ax[0,0].set_xlabel('Time')
ax[0,0].set_ylabel('Abundance')

m_f, b_f = np.polyfit(trans_1[0][:,0], trans_1[0][:,1], 1, w=trans_1[0][:,2])

ax[0,1].scatter(trans_1[0][:,0], trans_1[0][:,1], s = trans_1[0][:,2]*1000, label = 'Full-sib', c='C1')
ax[0,1].plot(trans_1[0][:,0], b_f + m_f*trans_1[0][:,0])

ax[0,1].set_xlim([0, 0.6])
ax[0,1].set_ylim([0, 0.6])
ax[0,1].set_xlabel('Mean Susceptibility to Endemic')
ax[0,1].set_ylabel('Susceptibility to Foreign')

R_plot = ax[0,2].imshow(def_slope, extent=[0,1,0,1], norm=plt.Normalize(-1,1), cmap='coolwarm', aspect='auto')
label_axes(ax[0,2], 's')

ax[1,0].plot(t2, S_2c.T, label = labels)

ax[1,0].legend()
ax[1,0].set_xlabel('Time')
ax[1,0].set_ylabel('Abundance')

m_f, b_f = np.polyfit(trans_2[0][:,0], trans_2[0][:,1], 1, w=trans_2[0][:,2])

ax[1,1].scatter(trans_2[0][:,0], trans_2[0][:,1], s = trans_2[0][:,2]*1000, label = 'Full-sib', c='C1')
ax[1,1].plot(trans_2[0][:,0], b_f + m_f*trans_2[0][:,0])

ax[1,1].set_xlim([0, 0.6])
ax[1,1].set_ylim([0, 0.6])
ax[1,1].set_xlabel('Mean Susceptibility to Endemic')
ax[1,1].set_ylabel('Susceptibility to Foreign')

ax[1,2].imshow(fxt_slope, extent=[0,1,0,1], norm=plt.Normalize(-1,1), cmap='coolwarm', aspect='auto')
label_axes(ax[1,2], 's')

ax[0,0].annotate('Recombination', xy=(0, 0.5), xytext=(-50,0), xycoords='axes fraction', textcoords='offset points', ha='center', va='center', rotation=90, fontsize=style.bigger)
ax[1,0].annotate('No Recombination', xy=(0, 0.5), xytext=(-50,0), xycoords='axes fraction', textcoords='offset points', ha='center', va='center', rotation=90, fontsize=style.bigger)


ax[0,0].annotate("A", xy=(-0.2, 1.05), xycoords="axes fraction", fontsize=style.bigger)
ax[0,1].annotate("B", xy=(-0.2, 1.05), xycoords="axes fraction", fontsize=style.bigger)
ax[0,2].annotate("C", xy=(-0.2, 1.05), xycoords="axes fraction", fontsize=style.bigger)
ax[1,0].annotate("D", xy=(-0.2, 1.05), xycoords="axes fraction", fontsize=style.bigger)
ax[1,1].annotate("E", xy=(-0.2, 1.05), xycoords="axes fraction", fontsize=style.bigger)
ax[1,2].annotate("F", xy=(-0.2, 1.05), xycoords="axes fraction", fontsize=style.bigger)

cbar = fig.colorbar(R_plot, ax=ax.ravel().tolist(), location='right', aspect=40, pad=0.025)
tick_locator = ticker.MaxNLocator(nbins=5)
cbar.locator = tick_locator
cbar.update_ticks()

cbar.set_label('Transitivity Slope', fontsize=style.medium)

fig.savefig('./figures/fig_4.svg', bbox_inches='tight', pad_inches=0.1)
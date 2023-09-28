import numpy as np

from model import Model, collapse_locus
from solve import run_sim
from utilities import load_data, get_trans, check_stab

from matplotlib import pyplot as plt
from matplotlib import ticker

import style
from style import label_axes

params_0 = {    'k':0.001, 
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

gs_path = './data/cov_gs_rho0.p'
gs_mask = check_stab('./data/cov_gs_rho1.p', 'c_s', 'c_g') 

_, _, gs_slope, _ = load_data(gs_path, 'c_s', 'c_g') 

gs_masked = np.ma.masked_where(gs_mask == 1, gs_slope)

sim = Model(**params_0)
t, S, I = run_sim(sim, np.ones(8)*10, np.array([10,10,0]), t=(0, 10000))

S_0 = S[:,-1]
I_0 = I[:,-1]

sim = Model(**params)
t, S, I = run_sim(sim, S_0, I_0, t=(0, 150000))

Sc, _ = collapse_locus(sim, S, 0)

labels = [r'$\it{gs}$', r'$\it{gS}$', r'$\it{Gs}$', r'$\it{GS}$']

fig, ax = plt.subplots(nrows = 1, ncols = 3, figsize=(10,3))
fig.tight_layout()
plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.4, hspace=0.4)

f_R = np.sum(S[sim.G[:,0] == 1], axis=0) / np.sum(S, axis=0)
f_G = np.sum(S[sim.G[:,1] == 1], axis=0) / np.sum(S, axis=0)
f_S = np.sum(S[sim.G[:,2] == 1], axis=0) / np.sum(S, axis=0)

ax[0].plot(t, f_R, label = 'Linkage Modifier')
ax[0].plot(t, f_G, label = 'General Resistance')
ax[0].plot(t, f_S, label = 'Specific Resistance')

#ax[0].legend(loc='center right')
ax[0].set_xlabel('Time')
ax[0].set_ylabel('Frequency')
ax[0].ticklabel_format(axis='x', style='sci', scilimits=(0,0))

trans_slope = np.zeros(len(t))

for i in range(len(t)):
    f_sib, h_sib = get_trans(sim, S[:,i], I[:,i])
    trans_slope[i], _ = np.polyfit(f_sib[:,0], f_sib[:,1], 1, w=f_sib[:,2])

ax[1].plot(t, trans_slope)
ax[1].set_xlabel('Time')
ax[1].set_ylabel('Transitivity Slope')
ax[1].ticklabel_format(axis='x', style='sci', scilimits=(0,0))

R_plot = ax[2].imshow(gs_masked, extent=[0,1,0,1], norm=plt.Normalize(-1,1), cmap='coolwarm', aspect='auto')
label_axes(ax[2], 's')

cbar = fig.colorbar(R_plot, ax=ax.ravel().tolist(), location='right', aspect=40, pad=0.01)
tick_locator = ticker.MaxNLocator(nbins=5)
cbar.locator = tick_locator
cbar.update_ticks()

cbar.set_label('Transitivity Slope', fontsize=style.medium)

ax[0].annotate("A", xy=(-0.2, 1.05), xycoords="axes fraction", fontsize=style.bigger)
ax[1].annotate("B", xy=(-0.2, 1.05), xycoords="axes fraction", fontsize=style.bigger)
ax[2].annotate("C", xy=(-0.2, 1.05), xycoords="axes fraction", fontsize=style.bigger)

fig.savefig('./figures/fig_5.svg', bbox_inches='tight', pad_inches=0.1)
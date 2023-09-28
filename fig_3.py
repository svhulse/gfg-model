import numpy as np

from utilities import load_data
from matplotlib import pyplot as plt
import style 
from style import label_axes

nocov_gs_rho0 = './data/nocov_gs_rho0.p'
gs_cov_rho0 = './data/cov_gs_rho0.p'

allele_freqs_nc, _, _, _ = load_data(nocov_gs_rho0, 'c_s', 'c_g')
allele_freqs_gs, _, _, _ = load_data(gs_cov_rho0, 'c_s', 'c_g') 

fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(5,3))
fig.tight_layout()
plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.25)

NoCov_plot = ax[0].imshow(allele_freqs_nc[1,:,:], extent=[0,1,0,1], norm=plt.Normalize(0,1), cmap='plasma')
Cov_plot = ax[1].imshow(allele_freqs_gs[1,:,:], extent=[0,1,0,1], norm=plt.Normalize(0,1), cmap='plasma')

label_axes(ax[0], 's')
label_axes(ax[1], 's')

ax[0].annotate("A", xy=(-0.2, 1.05), xycoords="axes fraction", fontsize=style.bigger)
ax[1].annotate("B", xy=(-0.2, 1.05), xycoords="axes fraction", fontsize=style.bigger)

cbar = fig.colorbar(Cov_plot, ax=ax.ravel().tolist(), location='bottom', aspect=40, pad=0.2)
cbar.set_label('General Resistance Allele Frequency', fontsize=style.medium)

fig.savefig('./figures/fig_3.svg', bbox_inches='tight', pad_inches=0.1)
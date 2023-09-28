import numpy as np

from utilities import load_data
from matplotlib import pyplot as plt

import style
from style import label_axes

nc_path = './data/nocov_gs.p'
gs_path = './data/cov_gs.p'
gv_path = './data/cov_gv.p'

allele_freqs_nc, _, _, _ = load_data(nc_path, 'c_s', 'c_g')
allele_freqs_gs, _, _, _ = load_data(gs_path, 'c_s', 'c_g') 
allele_freqs_gv, _, _, _ = load_data(gv_path, 'v', 'c_g')

fig, ax = plt.subplots(nrows=2, ncols=3, figsize=(8,6))
fig.tight_layout()
plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.25)

Q_nc_plot = ax[0,0].imshow(allele_freqs_nc[1,:,:], extent=[0,1,0,1], norm=plt.Normalize(0,1), cmap='plasma')
Q_plot = ax[0,1].imshow(allele_freqs_gs[1,:,:], extent=[0,1,0,1], norm=plt.Normalize(0,1), cmap='plasma')
Q_vq_plot = ax[0,2].imshow(allele_freqs_gv[1,:,:], extent=[0,1,0,1], norm=plt.Normalize(0,1), cmap='plasma')

R_nc_plot = ax[1,0].imshow(allele_freqs_nc[2,:,:], extent=[0,1,0,1], norm=plt.Normalize(0,1), cmap='plasma')
R_plot = ax[1,1].imshow(allele_freqs_gs[2,:,:], extent=[0,1,0,1], norm=plt.Normalize(0,1), cmap='plasma')
R_vq_plot = ax[1,2].imshow(allele_freqs_gv[2,:,:], extent=[0,1,0,1], norm=plt.Normalize(0,1), cmap='plasma')

ax[0,0].set_title('No Coevolution \n Resistance Costs', fontsize=style.bigger, pad=25)
ax[0,1].set_title('Coevolution \n Resistance Costs', fontsize=style.bigger, pad=25)
ax[0,2].set_title('Coevolution \n Virulence Costs', fontsize=style.bigger, pad=25)

ax[0,0].annotate('General Resistance', xy=(0, 0.5), xytext=(-50,0), xycoords='axes fraction', textcoords='offset points', ha='center', va='center', rotation=90, fontsize=style.bigger)
ax[1,0].annotate('Specific Resistance', xy=(0, 0.5), xytext=(-50,0), xycoords='axes fraction', textcoords='offset points', ha='center', va='center', rotation=90, fontsize=style.bigger)

label_axes(ax[0,0], 's')
label_axes(ax[0,1], 's')
label_axes(ax[0,2], 'v')
label_axes(ax[1,0], 's')
label_axes(ax[1,1], 's')
label_axes(ax[1,2], 'v')

ax[0,0].annotate("A", xy=(-0.2, 1.05), xycoords="axes fraction", fontsize=style.bigger)
ax[0,1].annotate("B", xy=(-0.2, 1.05), xycoords="axes fraction", fontsize=style.bigger)
ax[0,2].annotate("C", xy=(-0.2, 1.05), xycoords="axes fraction", fontsize=style.bigger)
ax[1,0].annotate("D", xy=(-0.2, 1.05), xycoords="axes fraction", fontsize=style.bigger)
ax[1,1].annotate("E", xy=(-0.2, 1.05), xycoords="axes fraction", fontsize=style.bigger)
ax[1,2].annotate("F", xy=(-0.2, 1.05), xycoords="axes fraction", fontsize=style.bigger)

cbar = fig.colorbar(Q_plot, ax=ax.ravel().tolist(), location='bottom', aspect=40, pad=0.1)
cbar.set_label('Allele Frequency', fontsize=style.medium)

fig.savefig('./figures/fig_2.svg', bbox_inches='tight', pad_inches=0.1)
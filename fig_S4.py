import numpy as np

from utilities import load_data
from matplotlib import pyplot as plt

import style
from style import label_axes

nc_dd = './data/nocov_gs_dd.p'
gs_dd = './data/cov_gs_dd.p'
gv_dd = './data/cov_gv_dd.p'
allele_freqs_ncdd, _, _, _ = load_data(nc_dd, 'c_s', 'c_g')
allele_freqs_gsdd, _, _, _ = load_data(gs_dd, 'c_s', 'c_g') 
allele_freqs_gvdd, _, _, _ = load_data(gv_dd, 'v', 'c_g')

nc_hs = './data/nocov_gs_hs.p'
gs_hs = './data/cov_gs_hs.p'
gv_hs = './data/cov_gv_hs.p'
allele_freqs_nchs, _, _, _ = load_data(nc_hs, 'c_s', 'c_g')
allele_freqs_gshs, _, _, _ = load_data(gs_hs, 'c_s', 'c_g') 
allele_freqs_gvhs, _, _, _ = load_data(gv_hs, 'v', 'c_g')

nc_sg = './data/nocov_gs_sg.p'
gs_sg = './data/cov_gs_sg.p'
gv_sg = './data/cov_gv_sg.p'
allele_freqs_ncsg, _, _, _ = load_data(nc_sg, 'c_s', 'c_g')
allele_freqs_gssg, _, _, _ = load_data(gs_sg, 'c_s', 'c_g') 
allele_freqs_gvsg, _, _, _ = load_data(gv_sg, 'v', 'c_g')

nc_hr = './data/nocov_gs_hr.p'
gs_hr = './data/cov_gs_hr.p'
gv_hr = './data/cov_gv_hr.p'
allele_freqs_nchr, _, _, _ = load_data(nc_hr, 'c_s', 'c_g')
allele_freqs_gshr, _, _, _ = load_data(gs_hr, 'c_s', 'c_g') 
allele_freqs_gvhr, _, _, _ = load_data(gv_hr, 'v', 'c_g')

fig, ax = plt.subplots(nrows=4, ncols=3, figsize=(8,12))
fig.tight_layout()
plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.25)

Q_nc_plot = ax[0,0].imshow(allele_freqs_ncdd[1,:,:], extent=[0,1,0,1], norm=plt.Normalize(0,1), cmap='plasma')
Q_plot = ax[0,1].imshow(allele_freqs_gsdd[1,:,:], extent=[0,1,0,1], norm=plt.Normalize(0,1), cmap='plasma')
Q_vq_plot = ax[0,2].imshow(allele_freqs_gvdd[1,:,:], extent=[0,1,0,1], norm=plt.Normalize(0,1), cmap='plasma')

ax[1,0].imshow(allele_freqs_nchs[1,:,:], extent=[0,1,0,1], norm=plt.Normalize(0,1), cmap='plasma')
ax[1,1].imshow(allele_freqs_gshs[1,:,:], extent=[0,1,0,1], norm=plt.Normalize(0,1), cmap='plasma')
ax[1,2].imshow(allele_freqs_gvhs[1,:,:], extent=[0,1,0,1], norm=plt.Normalize(0,1), cmap='plasma')

ax[2,0].imshow(allele_freqs_ncsg[1,:,:], extent=[0,1,0,1], norm=plt.Normalize(0,1), cmap='plasma')
ax[2,1].imshow(allele_freqs_gssg[1,:,:], extent=[0,1,0,1], norm=plt.Normalize(0,1), cmap='plasma')
ax[2,2].imshow(allele_freqs_gvsg[1,:,:], extent=[0,1,0,1], norm=plt.Normalize(0,1), cmap='plasma')

ax[3,0].imshow(allele_freqs_nchr[1,:,:], extent=[0,1,0,1], norm=plt.Normalize(0,1), cmap='plasma')
ax[3,1].imshow(allele_freqs_gshr[1,:,:], extent=[0,1,0,1], norm=plt.Normalize(0,1), cmap='plasma')
ax[3,2].imshow(allele_freqs_gvhr[1,:,:], extent=[0,1,0,1], norm=plt.Normalize(0,1), cmap='plasma')

ax[0,0].set_title('No Coevolution \n Resistance Costs', fontsize=style.bigger, pad=25)
ax[0,1].set_title('Coevolution \n Resistance Costs', fontsize=style.bigger, pad=25)
ax[0,2].set_title('Coevolution \n Virulence Costs', fontsize=style.bigger, pad=25)

ax[0,0].annotate('Density Dependent', xy=(0, 0.5), xytext=(-50,0), xycoords='axes fraction', textcoords='offset points', ha='center', va='center', rotation=90, fontsize=style.bigger)
ax[1,0].annotate('Hard Selection', xy=(0, 0.5), xytext=(-50,0), xycoords='axes fraction', textcoords='offset points', ha='center', va='center', rotation=90, fontsize=style.bigger)
ax[2,0].annotate('Strong Resistance', xy=(0, 0.5), xytext=(-50,0), xycoords='axes fraction', textcoords='offset points', ha='center', va='center', rotation=90, fontsize=style.bigger)
ax[3,0].annotate('High Recombination', xy=(0, 0.5), xytext=(-50,0), xycoords='axes fraction', textcoords='offset points', ha='center', va='center', rotation=90, fontsize=style.bigger)

label_axes(ax[0,0], 's')
label_axes(ax[0,1], 's')
label_axes(ax[0,2], 'v')
label_axes(ax[1,0], 's')
label_axes(ax[1,1], 's')
label_axes(ax[1,2], 'v')
label_axes(ax[2,0], 's')
label_axes(ax[2,1], 's')
label_axes(ax[2,2], 'v')
label_axes(ax[3,0], 's')
label_axes(ax[3,1], 's')
label_axes(ax[3,2], 'v')

ax[0,0].annotate("A", xy=(-0.2, 1.05), xycoords="axes fraction", fontsize=style.bigger)
ax[0,1].annotate("B", xy=(-0.2, 1.05), xycoords="axes fraction", fontsize=style.bigger)
ax[0,2].annotate("C", xy=(-0.2, 1.05), xycoords="axes fraction", fontsize=style.bigger)
ax[1,0].annotate("D", xy=(-0.2, 1.05), xycoords="axes fraction", fontsize=style.bigger)
ax[1,1].annotate("E", xy=(-0.2, 1.05), xycoords="axes fraction", fontsize=style.bigger)
ax[1,2].annotate("F", xy=(-0.2, 1.05), xycoords="axes fraction", fontsize=style.bigger)
ax[2,0].annotate("G", xy=(-0.2, 1.05), xycoords="axes fraction", fontsize=style.bigger)
ax[2,1].annotate("H", xy=(-0.2, 1.05), xycoords="axes fraction", fontsize=style.bigger)
ax[2,2].annotate("I", xy=(-0.2, 1.05), xycoords="axes fraction", fontsize=style.bigger)
ax[3,0].annotate("J", xy=(-0.2, 1.05), xycoords="axes fraction", fontsize=style.bigger)
ax[3,1].annotate("K", xy=(-0.2, 1.05), xycoords="axes fraction", fontsize=style.bigger)
ax[3,2].annotate("L", xy=(-0.2, 1.05), xycoords="axes fraction", fontsize=style.bigger)

cbar = fig.colorbar(Q_plot, ax=ax.ravel().tolist(), location='bottom', aspect=40, pad=0.05)
cbar.set_label('Allele Frequency', fontsize=style.medium)

fig.savefig('./figures/fig_S4.svg', bbox_inches='tight', pad_inches=0.1)
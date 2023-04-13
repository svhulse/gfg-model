import numpy as np
import pickle as pkl

from utilities import load_data

import matplotlib
from matplotlib import pyplot as plt

SMALL_SIZE = 8
MEDIUM_SIZE = 10
BIGGER_SIZE = 14

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=SMALL_SIZE)     # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

plt.rcParams["font.family"] = "Helvetica Neue"
matplotlib.rcParams['figure.dpi'] = 300

nc_path = 'Data/High_G/nocov_gs_costs.p'
gs_path = 'Data/High_G/cov_gs_costs.p'
gv_path = 'Data/High_G/cov_gv_costs.p'

G_nc, S_nc, _, _, D_1 = load_data(nc_path, 'c_s', 'c_g')
G_gs, S_gs, _, _, D_2 = load_data(gs_path, 'c_s', 'c_g') 
G_gv, S_gv, _, _, D_3 = load_data(gv_path, 'v', 'c_g')

def label_axes(ax, type):
    ax.invert_yaxis()
    ax.set_xlabel(r'$\it{G}$ cost')

    n_ticks = 5
    x_range = [0, 0.2]

    if type == 's':
        ax.set_ylabel(r'$\it{S}$ costs')
        y_range = [0, 0.4]
    elif type == 'v':
        ax.set_ylabel(r'$\it{V}$ cost')
        y_range = [0, 0.3]

    x_ticks = np.linspace(x_range[0], x_range[1], n_ticks)
    y_ticks = np.linspace(y_range[1], y_range[0], n_ticks)
    
    ax.set_yticks(np.linspace(0, 1, n_ticks))
    ax.set_yticklabels([str(tick) for tick in np.round(y_ticks, 4)])
    ax.set_xticks(np.linspace(0, 1, n_ticks))
    ax.set_xticklabels([str(tick) for tick in np.round(x_ticks, 4)])

fig, ax = plt.subplots(nrows=2, ncols=3, figsize=(8,6))
fig.tight_layout()
plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.25)

Q_nc_plot = ax[0,0].imshow(G_nc, extent=[0,1,0,1], norm=plt.Normalize(0,1), cmap='plasma')
Q_plot = ax[0,1].imshow(G_gs, extent=[0,1,0,1], norm=plt.Normalize(0,1), cmap='plasma')
Q_vq_plot = ax[0,2].imshow(G_gv, extent=[0,1,0,1], norm=plt.Normalize(0,1), cmap='plasma')

R_nc_plot = ax[1,0].imshow(S_nc, extent=[0,1,0,1], norm=plt.Normalize(0,1), cmap='plasma')
R_plot = ax[1,1].imshow(S_gs, extent=[0,1,0,1], norm=plt.Normalize(0,1), cmap='plasma')
R_vq_plot = ax[1,2].imshow(S_gv, extent=[0,1,0,1], norm=plt.Normalize(0,1), cmap='plasma')

ax[0,0].set_title('No Coevolution \n Resistance Costs', fontsize=BIGGER_SIZE, pad=25)
ax[0,1].set_title('Coevolution \n Resistance Costs', fontsize=BIGGER_SIZE, pad=25)
ax[0,2].set_title('Coevolution \n Virulence Costs', fontsize=BIGGER_SIZE, pad=25)

ax[0,0].annotate('General Resistance', xy=(0, 0.5), xytext=(-50,0), xycoords='axes fraction', textcoords='offset points', ha='center', va='center', rotation=90, fontsize=BIGGER_SIZE)
ax[1,0].annotate('Specific Resistance', xy=(0, 0.5), xytext=(-50,0), xycoords='axes fraction', textcoords='offset points', ha='center', va='center', rotation=90, fontsize=BIGGER_SIZE)

label_axes(ax[0,0], 's')
label_axes(ax[0,1], 's')
label_axes(ax[0,2], 'v')
label_axes(ax[1,0], 's')
label_axes(ax[1,1], 's')
label_axes(ax[1,2], 'v')

ax[0,0].annotate("A", xy=(-0.2, 1.05), xycoords="axes fraction", fontsize=BIGGER_SIZE)
ax[0,1].annotate("B", xy=(-0.2, 1.05), xycoords="axes fraction", fontsize=BIGGER_SIZE)
ax[0,2].annotate("C", xy=(-0.2, 1.05), xycoords="axes fraction", fontsize=BIGGER_SIZE)
ax[1,0].annotate("D", xy=(-0.2, 1.05), xycoords="axes fraction", fontsize=BIGGER_SIZE)
ax[1,1].annotate("E", xy=(-0.2, 1.05), xycoords="axes fraction", fontsize=BIGGER_SIZE)
ax[1,2].annotate("F", xy=(-0.2, 1.05), xycoords="axes fraction", fontsize=BIGGER_SIZE)

fig.colorbar(Q_plot, ax=ax.ravel().tolist(), location='bottom', aspect=40, pad=0.1)

fig.savefig('fig_S4.svg', bbox_inches='tight', pad_inches=0.1)
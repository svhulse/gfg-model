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
plt.rc('axes', labelsize=SMALL_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

plt.rcParams["font.family"] = "Helvetica Neue"
matplotlib.rcParams['figure.dpi'] = 300

gs_path = 'Data/Default/cov_gs_costs.p'
gv_path = 'Data/Default/cov_gv_costs.p'

_, _, _, gs_slope, _ = load_data(gs_path, 'c_s', 'c_g') 
_, _, _, gv_slope, _ = load_data(gv_path, 'v', 'c_g')

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

fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(6,6))
fig.tight_layout()
plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.3, hspace=0.25)

ax[0].set_title('Effect of Resistance Costs', fontsize=BIGGER_SIZE, pad=25)
R_plot = ax[0].imshow(gs_slope, extent=[0,1,0,1], norm=plt.Normalize(0,1), cmap='plasma')
label_axes(ax[0], 's')

ax[1].set_title('Effect of Virulence Costs', fontsize=BIGGER_SIZE, pad=25)
T_plot = ax[1].imshow(gv_slope, extent=[0,1,0,1], norm=plt.Normalize(0,1), cmap='plasma')#, alpha=F_region)
label_axes(ax[1], 'v')

ax[0].annotate("A", xy=(-0.2, 1.05), xycoords="axes fraction", fontsize=BIGGER_SIZE)
ax[1].annotate("B", xy=(-0.2, 1.05), xycoords="axes fraction", fontsize=BIGGER_SIZE)

ax[0].annotate("*", xy=(0.5, 0.45), xycoords="axes fraction", fontsize=BIGGER_SIZE)

fig.colorbar(T_plot, ax=ax.ravel().tolist(), location='bottom', aspect=40, pad=0.1)

fig.savefig('fig_3.svg', bbox_inches='tight', pad_inches=0.1)
import numpy as np

from utilities import check_stab

from matplotlib import pyplot as plt
import matplotlib.patches as mpatches

import style
from style import label_axes

gs_path = './data/cov_gs_rho1.p'
gv_path = './data/cov_gv_rho1.p'
gs_path2 = './data/cov_gs_rho0.p'
gv_path2 = './data/cov_gv_rho0.p'

stab_gs = check_stab(gs_path, 'c_s', 'c_g') 
stab_gv = check_stab(gv_path, 'v', 'c_g')
stab_gs2 = check_stab(gs_path2, 'c_s', 'c_g') 
stab_gv2 = check_stab(gv_path2, 'v', 'c_g')

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

fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(5,5.75))
fig.tight_layout()
plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.3, hspace=0.1)

im=ax[0,0].imshow(stab_gs, extent=[0,1,0,1], norm=plt.Normalize(0,1), cmap='plasma')
label_axes(ax[0,0], 's')

ax[0,1].imshow(stab_gv, extent=[0,1,0,1], norm=plt.Normalize(0,1), cmap='plasma')
label_axes(ax[0,1], 'v')

ax[1,0].imshow(stab_gs2, extent=[0,1,0,1], norm=plt.Normalize(0,1), cmap='plasma')
label_axes(ax[1,0], 's')

ax[1,1].imshow(stab_gv2, extent=[0,1,0,1], norm=plt.Normalize(0,1), cmap='plasma')
label_axes(ax[1,1], 'v')

values = np.unique(stab_gs.ravel())
colors = [ im.cmap(im.norm(value)) for value in values]
patches = [ mpatches.Patch(color=colors[i], label="Level {l}".format(l=values[i]) ) for i in range(len(values)) ]
labels = ['Center or Unstable', 'Asymptotically Stable']

ax[0,0].annotate('Modifier Lost', xy=(0, 0.5), xytext=(-50,0), xycoords='axes fraction', textcoords='offset points', ha='center', va='center', rotation=90, fontsize=style.bigger)
ax[1,0].annotate('Modifier Fixed', xy=(0, 0.5), xytext=(-50,0), xycoords='axes fraction', textcoords='offset points', ha='center', va='center', rotation=90, fontsize=style.bigger)

ax[0,0].annotate("A", xy=(-0.2, 1.05), xycoords="axes fraction", fontsize=style.bigger)
ax[0,1].annotate("B", xy=(-0.2, 1.05), xycoords="axes fraction", fontsize=style.bigger)
ax[1,0].annotate("C", xy=(-0.2, 1.05), xycoords="axes fraction", fontsize=style.bigger)
ax[1,1].annotate("D", xy=(-0.2, 1.05), xycoords="axes fraction", fontsize=style.bigger)

fig.legend(handles=patches, labels=labels, loc="lower center", ncol=2)

fig.savefig('./figures/fig_S2.svg', bbox_inches='tight', pad_inches=0.1)
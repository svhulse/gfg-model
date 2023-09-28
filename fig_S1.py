import numpy as np
import matplotlib

from utilities import check_stab
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches

import style
from style import label_axes

nc_path = './data/nocov_gs.p'
gs_path = './data/cov_gs.p'
gv_path = './data/cov_gv.p'

stab_nc = check_stab(nc_path, 'c_s', 'c_g')
stab_gs = check_stab(gs_path, 'c_s', 'c_g') 
stab_gv = check_stab(gv_path, 'v', 'c_g')

fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(8,3.25))
fig.tight_layout()
plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.3, hspace=0.25)

im = ax[0].imshow(stab_nc, extent=[0,1,0,1], norm=plt.Normalize(0,1), cmap='plasma')
label_axes(ax[0], 's')

ax[1].imshow(stab_gs, extent=[0,1,0,1], norm=plt.Normalize(0,1), cmap='plasma')#, alpha=F_region)
label_axes(ax[1], 's')

ax[2].imshow(stab_gv, extent=[0,1,0,1], norm=plt.Normalize(0,1), cmap='plasma')#, alpha=F_region)
label_axes(ax[2], 'v')

values = np.unique(stab_nc.ravel())
colors = [ im.cmap(im.norm(value)) for value in values]
patches = [ mpatches.Patch(color=colors[i], label="Level {l}".format(l=values[i]) ) for i in range(len(values)) ]
labels = ['Center or Unstable', 'Asymptotically Stable']

ax[0].annotate("A", xy=(-0.2, 1.05), xycoords="axes fraction", fontsize=style.bigger)
ax[1].annotate("B", xy=(-0.2, 1.05), xycoords="axes fraction", fontsize=style.bigger)
ax[2].annotate("C", xy=(-0.2, 1.05), xycoords="axes fraction", fontsize=style.bigger)

fig.legend(handles=patches, labels=labels, loc="lower center", ncol=2)

fig.savefig('./figures/fig_S1.svg', bbox_inches='tight', pad_inches=0.1)
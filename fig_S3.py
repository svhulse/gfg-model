from utilities import load_data
from matplotlib import pyplot as plt

import style
from style import label_axes

gs_path = 'data/cov_gs.p'
gv_path = 'data/cov_gv.p'

_, _, _, D_gs = load_data(gs_path, 'c_s', 'c_g') 
_, _, _, D_gv = load_data(gv_path, 'v', 'c_g')


fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(5,3))
fig.tight_layout()
plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.25, hspace=0.25)

plot = ax[0].imshow(D_gs, extent=[0,1,0,1], norm=plt.Normalize(0,0.6), cmap='plasma')
ax[1].imshow(D_gv, extent=[0,1,0,1], norm=plt.Normalize(0,0.6), cmap='plasma')

label_axes(ax[0], 's')
label_axes(ax[1], 'v')

ax[0].annotate("A", xy=(-0.2, 1.05), xycoords="axes fraction", fontsize=style.bigger)
ax[1].annotate("B", xy=(-0.2, 1.05), xycoords="axes fraction", fontsize=style.bigger)

cbar = fig.colorbar(plot, ax=ax.ravel().tolist(), location='bottom', aspect=40, pad=0.2)
cbar.set_label('Linkage Disequilibrium ($D^\prime$)', fontsize=style.medium)

fig.savefig('./figures/fig_S3.svg', bbox_inches='tight', pad_inches=0.1)
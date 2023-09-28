import matplotlib
import numpy as np

from matplotlib import pyplot as plt

small = 8
medium = 10
bigger = 14

plt.rc('font', size=small)          # controls default text sizes
plt.rc('axes', titlesize=small)     # fontsize of the axes title
plt.rc('axes', labelsize=small)     # fontsize of the x and y labels
plt.rc('xtick', labelsize=small)    # fontsize of the tick labels
plt.rc('ytick', labelsize=small)    # fontsize of the tick labels
plt.rc('legend', fontsize=small)    # legend fontsize
plt.rc('figure', titlesize=bigger)  # fontsize of the figure title

matplotlib.rcParams['font.family'] = ['Helvetica Neue']
matplotlib.rcParams['figure.dpi'] = 600

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
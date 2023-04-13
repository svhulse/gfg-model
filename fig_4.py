import numpy as np

from allelic_model import run_sim
from utilities import get_trans

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

params_1 = {    'k':0.001, 
                'mu':0.2,  
                'b':1,     
                'beta':0.5,
                'nh':0.1,
                'g':0.3,
                's':0.9,
                'rho':0.05, 
                'c_g':0.1,     
                'c_s':0.2,      
                'v':0.2,
                'S_0':[100,100,100,100],
                'I_0':[10,1,0]}   

params_2 = {    'k':0.001, 
                'mu':0.2,  
                'b':1,     
                'beta':0.5,
                'nh':0.1,
                'g':0.3,
                's':0.9,
                'rho':0.0, 
                'c_g':0.1,     
                'c_s':0.2,      
                'v':0.2,
                'S_0':[100,100,100,100],
                'I_0':[10,1,0]}      

t1, S_1, I_1 = run_sim(**params_1)
trans_1 = get_trans(S_1[:,-1], I_1[:,-1], **params_1)

t2, S_2, I_2 = run_sim(**params_2)
trans_2 = get_trans(S_2[:,-1], I_2[:,-1], **params_2)

labels = [r'$\it{GS}$', r'$\it{Gs}$', r'$\it{gS}$', r'$\it{gs}$']

fig, ax = plt.subplots(nrows = 2, ncols = 3, figsize=(8,6))
fig.tight_layout()
plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.25, hspace=0.25)

ax[0,0].plot(t1, S_1.T, label = labels)

ax[0,0].legend()
ax[0,0].set_xlabel('Time')
ax[0,0].set_ylabel('Abundance')
ax[0,0].set_title('Host Genotypes', fontsize=BIGGER_SIZE, pad=25)

ax[0,1].plot(t1, I_1[0,:]/(I_1[0,:] + I_1[1,:]), label=r'$\it{Avr}$')
ax[0,1].plot(t1, I_1[1,:]/(I_1[0,:] + I_1[1,:]), label=r'$\it{Vir}$')

ax[0,1].set_ylim([0, 1])
ax[0,1].set_xlabel('Time')
ax[0,1].set_ylabel('Frequency')
ax[0,1].legend()
ax[0,1].set_title('Pathogen Genotypes', fontsize=BIGGER_SIZE, pad=25)

m_f, b_f = np.polyfit(trans_1[0], trans_1[1], 1, w=trans_1[2])
m_h, b_h = np.polyfit(trans_1[3], trans_1[4], 1, w=trans_1[5])

ax[0,2].scatter(trans_1[0], trans_1[1], s = trans_1[2]*1000, label = 'Full-sib')
ax[0,2].scatter(trans_1[3], trans_1[4], s = trans_1[5]*1000, marker = '^', label = 'Half-sib')
ax[0,2].plot(trans_1[0], b_f + m_f*trans_1[0])
ax[0,2].plot(trans_1[3], b_h + m_h*trans_1[3])

ax[0,2].set_xlim([0, 0.6])
ax[0,2].set_ylim([0, 0.6])
ax[0,2].set_xlabel('Mean Susceptibility to Endemic')
ax[0,2].set_ylabel('Susceptibility to Foreign')
ax[0,2].set_title('Resistance Transitivity', fontsize=BIGGER_SIZE, pad=25)
ax[0,2].legend(loc=4)

ax[1,0].plot(t2, S_2.T, label = labels)


ax[1,0].legend()
ax[1,0].set_xlabel('Time')
ax[1,0].set_ylabel('Abundance')

ax[1,1].plot(t2, I_2[0,:]/(I_2[0,:] + I_2[1,:]), label=r'$\it{Avr}$')
ax[1,1].plot(t2, I_2[1,:]/(I_2[0,:] + I_2[1,:]), label=r'$\it{Vir}$')

ax[1,1].set_ylim([0, 1])
ax[1,1].set_xlabel('Time')
ax[1,1].set_ylabel('Frequency')
ax[1,1].legend()

m_f, b_f = np.polyfit(trans_2[0], trans_2[1], 1, w=trans_2[2])
m_h, b_h = np.polyfit(trans_2[3], trans_2[4], 1, w=trans_2[5])

ax[1,2].scatter(trans_2[0], trans_2[1], s = trans_2[2]*1000, label = 'Full-sib')
ax[1,2].scatter(trans_2[3], trans_2[4], s = trans_2[5]*1000, marker = '^', label = 'Half-sib')
ax[1,2].plot(trans_2[0], b_f + m_f*trans_2[0])
ax[1,2].plot(trans_2[3], b_h + m_h*trans_2[3])

ax[1,2].set_xlim([0, 0.6])
ax[1,2].set_ylim([0, 0.6])
ax[1,2].set_xlabel('Mean Susceptibility to Endemic')
ax[1,2].set_ylabel('Susceptibility to Foreign')
ax[1,2].legend(loc=4)

ax[0,0].annotate('Recombination', xy=(0, 0.5), xytext=(-50,0), xycoords='axes fraction', textcoords='offset points', ha='center', va='center', rotation=90, fontsize=BIGGER_SIZE)
ax[1,0].annotate('No Recombination', xy=(0, 0.5), xytext=(-50,0), xycoords='axes fraction', textcoords='offset points', ha='center', va='center', rotation=90, fontsize=BIGGER_SIZE)

ax[0,0].annotate("A", xy=(-0.2, 1.05), xycoords="axes fraction", fontsize=BIGGER_SIZE)
ax[0,1].annotate("B", xy=(-0.2, 1.05), xycoords="axes fraction", fontsize=BIGGER_SIZE)
ax[0,2].annotate("C", xy=(-0.2, 1.05), xycoords="axes fraction", fontsize=BIGGER_SIZE)
ax[1,0].annotate("D", xy=(-0.2, 1.05), xycoords="axes fraction", fontsize=BIGGER_SIZE)
ax[1,1].annotate("E", xy=(-0.2, 1.05), xycoords="axes fraction", fontsize=BIGGER_SIZE)
ax[1,2].annotate("F", xy=(-0.2, 1.05), xycoords="axes fraction", fontsize=BIGGER_SIZE)

fig.savefig('fig_4.svg', bbox_inches='tight', pad_inches=0.1)
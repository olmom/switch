"""
Illustration of ultrasensitivity
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint
from utils.signal_response import SignalResponseCurve
import os
plt.rcParams['text.usetex'] = True
plt.rcParams['font.size'] = '13'
plt.rcParams['legend.fontsize'] = '12'
plt.rcParams['xtick.labelsize'] = '12'
plt.rcParams['ytick.labelsize'] = '12'
plt.rcParams['axes.spines.top'] = False
plt.rcParams['axes.spines.right'] = False

# FIGURE 1A-C: ILLUSTRATION OF LINEAR VS MICHAELIAN VS SIGMOIDAL
# SIGNAL RESPONSE CURVE
fig1 = plt.figure(figsize=(11,3.7))

linear_obj = SignalResponseCurve(k=0.1)
nonlin_obj = SignalResponseCurve(k=1, K=2.5, n=5.3)

signal = np.arange(0, 10, 0.01)

response_linear  = linear_obj.linear(signal)
response_mmenten = nonlin_obj.michaelis_menten(signal)
nonlin_obj.K=4
response_sigmoid = nonlin_obj.sigmoidal(signal)

ax1 = fig1.add_subplot(151)
ax2 = fig1.add_subplot(152)
ax3 = fig1.add_subplot(153)
ax1.plot(signal, response_linear, c='k')
ax2.plot(signal, response_mmenten, c='k')
ax3.plot(signal, response_sigmoid, c='k')
ax1.set_title('linear'); ax2.set_title('Michaelis\nMenten'); ax3.set_title('sigmoidal')
ax1.set_xlabel('stimulus'); ax2.set_xlabel('stimulus');ax3.set_xlabel('stimulus');
ax1.set_ylabel('response'); ax2.set_ylabel('response'); ax3.set_ylabel('response')
ax1.set_ylim([0,1]); ax2.set_ylim([0,1]); ax3.set_ylim([0,1]);
ax1.set_xlim([0,10]); ax2.set_xlim([0,10]); ax3.set_xlim([0,10]);
ax1.set_xticklabels([]); ax2.set_xticklabels([]); ax3.set_xticklabels([])
ax1.set_yticklabels([]); ax2.set_yticklabels([]); ax3.set_yticklabels([])
ax1.set_aspect(1.00/ax1.get_data_ratio(), adjustable='box')
ax2.set_aspect(1.00/ax2.get_data_ratio(), adjustable='box')
ax3.set_aspect(1.00/ax3.get_data_ratio(), adjustable='box')

##########################################

# FIGURE 1D,E: BISTABLE SWITCHES
stimulus = np.arange(0, 0.1, 0.0001)
strength_feedback = np.array([0.14, 0.20])

t = np.arange(0,1000,0.1)
colors = ['k', 'r', 'b', 'g']
titles = ['reversible\nswitch', 
          'irreversible\nswitch']
titles2 = ['D', 'E']

for f in range(len(strength_feedback)):
    ax = fig1.add_subplot(1,5,f+4)

    # Bistability with ODEs
    responses, responses_inv = [], []
    y0 = 0
    for s in stimulus:
        obj = SignalResponseCurve(stimulus=s, f=strength_feedback[f], 
                                  n=5, d=0.01, A_tot=1, K=1)
        sol = odeint(obj.bistable_ODE, y0, t)
        SS = sol[-1]
        y0 = SS
        responses.append(SS)
    ax.plot(stimulus, responses, 
            label='$f={}$'.format(format(strength_feedback[f], '.2f')), 
        c='k')

    # Bistability with ODEs -- inverse
    y0 = 0
    for s in stimulus[::-1]:
        obj = SignalResponseCurve(stimulus=s, f=strength_feedback[f], 
                                  n=5, d=0.01, A_tot=1, K=1)
        sol = odeint(obj.bistable_ODE, y0, t)
        SS = sol[-1]
        y0 = SS
        responses_inv.append(SS)
    ax.plot(stimulus, responses_inv[::-1], c='k')

    # Bistability with steady state equation from
    # (Xiong & Ferrel -- corrigendum)
    obj = SignalResponseCurve(f=strength_feedback[f], 
                              n=5, d=0.01, A_tot=1, K=1)
    responses = np.asarray(responses).flatten()
    responses = np.linspace(responses.min(), responses.max(), 1000)
    input = obj.bistable(responses)
    ax.plot(input, responses, '--', color='k')
    ax.set_xlabel('stimulus')
    ax.set_ylabel('response')   
    ax.set_xlim([0,.025/2.5]); ax.set_xticks([0,0.025/5, 0.025/2.5])
    ax.set_ylim([0,0.93]); ax.set_yticks([0,0.93/2, 0.93])
    ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')
    ax.set_title(titles[f])
    ax.set_title(titles2[f], loc='left')
    ax.set_xticklabels([]); ax.set_yticklabels([])

ax1.set_title('A', loc='left')
ax2.set_title('B', loc='left')
ax3.set_title('C', loc='left')

fig1.subplots_adjust(
top=0.88,
bottom=0.11,
left=0.07,
right=0.93,
hspace=0.59,
wspace=0.6
)

isExist = os.path.exists('./figures/')
if not isExist:  
    os.makedirs('./figures/')
#fig1.savefig('./figures/fig1.pdf', format='pdf')
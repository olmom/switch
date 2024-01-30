"""
dxdt = f(x) = -x(x-1)(x-2) + p 
illustration of steady states, bifurcation diagrams, hysteresis
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
plt.rcParams['lines.linewidth'] = 2

t = np.arange(0, 1000, .1)

stimulus, K = 0.0025, 0.95 
strength_feedback = np.array([0, 0.15])

fig3 = plt.figure(figsize=(11,5))

ax0 = fig3.add_subplot(251)
ax00 = fig3.add_subplot(256)
ax0.axis('off'); ax00.axis('off')

# FIGURE 3A-D: MONOSTABLE SYSTEM
# ------------------------------
for f in range(len(strength_feedback)):
    bistable_obj = SignalResponseCurve(stimulus=stimulus, 
                                       f=strength_feedback[f], 
                                       n=5, d=0.01, A_tot=1.0, K=K)
    sol1 = odeint(bistable_obj.bistable_ODE, .01, t)
    sol2 = odeint(bistable_obj.bistable_ODE, 2, t)
    sol3 = odeint(bistable_obj.bistable_ODE, 0.6, t)

    # time series: A* approaches the steady state over time
    idx_fig = 2 if f==0 else 7 
    ax1 = fig3.add_subplot(2,5,idx_fig)
    ax1.plot(t, sol1, c='#1B9E77'); ax1.plot(t, sol2, c='#D95F02'); 
    ax1.plot(t, sol3, c='#7570B3')
    ax1.set_ylabel('conc. $A^*$'); ax1.set_xlabel('time')
    ax1.set_xlim([0,500])
    ax1.set_xticklabels([]);ax1.set_yticklabels([]) 

    # calculate steady state graphically
    ax2 = fig3.add_subplot(2,5,idx_fig+1)
    x = np.linspace(0, 3, 1000)
    prod = bistable_obj.production_term(x) + bistable_obj.PFL_term(x)
    deg  = bistable_obj.degradation_term(x)
    ax2.plot(x, prod, label='$A^*$ activation', c='goldenrod')
    ax2.plot(x, deg, label='$A^*$ inactivation', c='dimgray')
    ax2.legend(framealpha=0, loc='center left', 
               bbox_to_anchor=(1, 0.5))
    ax2.set_ylabel('term');ax2.set_xlabel('conc. $A^*$')
    ax2.set_ylim([0, 0.01]); ax2.set_xlim([0,1])
    ax2.set_xticklabels([]);ax2.set_yticklabels([]) 

    if f == 1:
        ax1.set_xlim([0,250])
        ax1.set_title('F', loc='left')
        ax2.set_title('G', loc='left')
    else:
        ax1.set_title('B', loc='left')
        ax2.set_title('C', loc='left')

    ax1.set_aspect(1/ax1.get_data_ratio(), adjustable='box')
    ax2.set_aspect(1/ax2.get_data_ratio(), adjustable='box')


# FIGURE 3E-F: BISTABLE SYSTEM WITH HYSTERESIS (POSITIVE FEEDBACK)
# ----------------------------------------------------------------
ax5 = fig3.add_subplot(255)
ax7 = fig3.add_subplot(2,5,10)
stimulus = np.arange(0, .01, 1e-5)

y0 = 0.1
response, response_inv = [], []
# increase stimulus in one direction
for s in stimulus:
    bistable_obj = SignalResponseCurve(stimulus=s, 
                                       f=strength_feedback[1], 
                                       n=5, d=0.01, A_tot=1, K=1.)
    sol = odeint(bistable_obj.bistable_ODE, y0, t) 
    r = sol[-1]
    y0 = r
    response.append(r)
# decrease stimulus (in other direction, to see hysteresis curve)
for s in stimulus[::-1]:
    bistable_obj = SignalResponseCurve(stimulus=s, 
                                       f=strength_feedback[1], 
                                       n=5, d=0.01, A_tot=1, K=1.)
    sol = odeint(bistable_obj.bistable_ODE, y0, t) 
    r = sol[-1]
    y0 = r
    response_inv.append(r)    
ax7.plot(stimulus, response, c='k', alpha=0.8)
ax7.plot(stimulus, response_inv[::-1], c='k', alpha=0.8)

bistable_obj = SignalResponseCurve(f=strength_feedback[1], 
                                   n=5, d=0.01, A_tot=1, K=1.)
responses = np.asarray(response).flatten()
responses = np.linspace(responses.min(), responses.max(), 10000)
inputs = bistable_obj.bistable(responses)
idxs_UFP = np.where(np.diff(inputs)< 0)[0]
ax7.plot(inputs[idxs_UFP], responses[idxs_UFP], '--', c='k')
ax7.set_xlabel('stimulus'); ax7.set_ylabel('response $A^*$')
ax7.set_xlim([0,0.01]); ax7.set_ylim([0, 0.88])
ax7.set_xticklabels([]); ax7.set_yticklabels([])

response = []
for s in stimulus:
    bistable_obj = SignalResponseCurve(stimulus=s, 
                                       f=strength_feedback[0], 
                                       n=5, d=0.01, A_tot=1, K=1.)
    sol = odeint(bistable_obj.bistable_ODE, y0, t) 
    r = sol[-1]
    y0 = r
    response.append(r)

ax5.plot(stimulus, response, c='k', alpha=0.8)
responses = np.asarray(response).flatten()
responses = np.linspace(responses.min(), responses.max(), 10000)
inputs = bistable_obj.bistable(responses)
idxs_UFP = np.where(np.diff(inputs) < 0)[0]
ax5.plot(inputs[idxs_UFP], responses[idxs_UFP], '--', c='k')
ax5.set_xlabel('stimulus'); ax5.set_ylabel('response $A^*$')
ax5.set_xlim([0,0.01]); ax5.set_ylim([0, 0.88])
ax5.set_xticklabels([]); ax5.set_yticklabels([])

ax5.set_aspect(1/ax5.get_data_ratio(), adjustable='box')
ax7.set_aspect(1/ax7.get_data_ratio(), adjustable='box')
ax0.set_aspect(1/ax0.get_data_ratio(), adjustable='box')
ax00.set_aspect(1/ax00.get_data_ratio(), adjustable='box')

ax0.set_title('A', loc='left')
ax00.set_title('E', loc='left')
ax5.set_title('D', loc='left')
ax7.set_title('H', loc='left')

fig3.subplots_adjust(
top=0.935,
bottom=0.11,
left=0.075,
right=0.945,
hspace=0.1,
wspace=0.3
)

isExist = os.path.exists('./figures/')
if not isExist:  
    os.makedirs('./figures/')
#fig3.savefig('./figures/fig3_2.pdf', format='pdf')
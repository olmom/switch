"""
add panel with prod vs degradation where we see the SS
period as a function of delay
with DDEs -> better in R?
"""
import matplotlib.pyplot as plt
import numpy as np
from utils.delay_equation import DDE
import os
from scipy.signal import argrelextrema
from ddeint import ddeint
plt.rcParams['text.usetex'] = True
plt.rcParams['font.size'] = '13'
plt.rcParams['legend.fontsize'] = '12'
plt.rcParams['xtick.labelsize'] = '12'
plt.rcParams['ytick.labelsize'] = '12'
plt.rcParams['axes.spines.top'] = False
plt.rcParams['axes.spines.right'] = False


obj = DDE(beta=4, n=5, d=0.07, K=1)
g = lambda t : 1 #y(t) before integration interval
dt = 0.1
t = np.arange(0, 2000, dt)

delays = np.arange(3.8, 37.8, 1) 

fig4 = plt.figure(figsize=(11,7))
ax1 = fig4.add_subplot(241)
ax2 = fig4.add_subplot(242)
ax3 = fig4.add_subplot(245)
ax4 = fig4.add_subplot(246)
ax5 = fig4.add_subplot(247)
ax6 = fig4.add_subplot(248)

# PANEL B: oscillation period as a function of changing delay
periods = []
for tau in delays:
    #solve DDE
    obj.tau = tau
    sol = ddeint(obj.delay_model1, g, t)
    x  = sol
    #timeseries  
    x_idx_max = argrelextrema(x, np.greater)[0]
    per = np.diff(t[x_idx_max[-9:]]).mean()
    periods.append(per) 

ax2.plot(np.asarray(delays)/4, np.asarray(periods)/4)
ax2.set_yticks([0, 6, 12, 18, 24, 30])
ax2.set_xticks([0, 2, 4, 6, 8, 10])
ax2.set_xlabel(r'delay $\tau$ (h)'); ax2.set_ylabel('period (h)')
ax2.set_aspect(1/ax2.get_data_ratio(), adjustable='box')

# PANEL D: circadian timeseries
delay = 25.8
obj = DDE(beta=4, n=5, d=0.07, K=1, tau=delay)

sol = ddeint(obj.delay_model1, g, t)
ax4.plot(t[0:int(96.55*2/dt)], sol[-int(96.55*2/dt):], c='k')
ax4.set_xlabel('time (h)'); ax4.set_ylabel('conc. $x$')
ax4.set_xticks([0, 96, 192]); 
ax4.set_xticklabels([0, 24, 48]); ax4.set_yticklabels([])
ax4.set_xlim([0, 193.1]); 
ax4.set_title('for delay=${}$h'.format(format(delay/4, '.2f')))

# PANEL E: circadian phase portrait
last_5d = int(5*96.55/dt)
last_tau = int(obj.tau/dt)
x = sol[-last_5d:]
xd = sol[-last_5d-last_tau : (-last_tau)]
ax5.plot(x, xd, c='k')
ax5.set_xlabel('$x(t)$')
ax5.set_ylabel(r'$x(t-\tau)$')
ax5.set_xticklabels([]); ax5.set_yticklabels([])

# PANEL F: Hopf bifurcation for changing Hill exponent
obj = DDE(beta=4, tau=6.4, d=0.1, K=1)

g = lambda t : 1 #y(t) before integration interval
dt = 0.1
t = np.arange(0, 1000, dt)

hill_exponents = np.arange(1, 8, .10) 

periods, amplitudes, magnitudes = [], [], []
x_minima, x_maxima = [], []
std_devs_per = []
for n in hill_exponents:
    #solve DDE
    obj.n = n 
    sol = ddeint(obj.delay_model1, g, t)
    x  = sol

    #timeseries  
    x_idx_max = argrelextrema(x, np.greater)[0]
    x_idx_min = argrelextrema(x, np.less)[0]
    per = np.diff(t[x_idx_max[-9:]]).mean()
    per_sd = np.diff(t[x_idx_max[-9:]]).std()

    mean = x[-int(per*5/dt):].mean()
    x_min, x_max= x[x_idx_min[-9:]].mean(), x[x_idx_max[-9:]].mean()
    amp = (x_max - x_min)/mean
    periods.append(per); amplitudes.append(amp)
    x_maxima.append(x_max); x_minima.append(x_min)
    magnitudes.append(mean); std_devs_per.append(per_sd)

periods, magnitudes = np.asarray(periods), np.asarray(magnitudes)
x_minima, x_maxima = np.asarray(x_minima), np.asarray(x_maxima)
std_devs_per = np.asarray(std_devs_per)
amplitudes = np.asarray(amplitudes).flatten()

ax6.plot(hill_exponents, x_minima/magnitudes, 
         c='k', markersize=4)
ax6.plot(hill_exponents, x_maxima/magnitudes, 
         c='k', markersize=4)
ax6.set_xlabel('Hill exponent $n$'); 
ax6.set_ylabel('$x$min, max (norm.)')

###############################

ax1.axis('off'); ax3.axis('off')
ax1.set_aspect(1/ax1.get_data_ratio(), adjustable='box')
ax2.set_aspect(1/ax2.get_data_ratio(), adjustable='box')
ax3.set_aspect(1/ax3.get_data_ratio(), adjustable='box')
ax4.set_aspect(1/ax4.get_data_ratio(), adjustable='box')
ax5.set_aspect(1/ax5.get_data_ratio(), adjustable='box')
ax6.set_aspect(1/ax6.get_data_ratio(), adjustable='box')

ax1.set_title('A',loc='left')
ax2.set_title('B',loc='left')
ax3.set_title('C',loc='left')
ax4.set_title('D',loc='left')
ax5.set_title('E',loc='left')
ax6.set_title('F',loc='left')

fig4.subplots_adjust(top=0.88,
bottom=0.1,
left=0.045,
right=0.955,
hspace=0.8,
wspace=0.2)

fig4.show()
isExist = os.path.exists('./figures/')
if not isExist:  
    os.makedirs('./figures/')
#fig4.savefig('./figures/fig4.pdf', format='pdf')
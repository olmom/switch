"""
Gonze 2013: equations, time series, phase portrait
small d for 24h oscillations, n=10
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint
from utils.goodwin import GoodwinOscillator
import pandas as pd
import os
from scipy.signal import argrelextrema
plt.rcParams['text.usetex'] = True
plt.rcParams['font.size'] = '13'
plt.rcParams['legend.fontsize'] = '12'
plt.rcParams['xtick.labelsize'] = '12'
plt.rcParams['ytick.labelsize'] = '12'
plt.rcParams['axes.spines.top'] = False
plt.rcParams['axes.spines.right'] = False

# SOLVE GOODWIN SYSTEM
y0 = [.1, .1, .1]
dt = .01
t = np.arange(0, 10000, dt)

obj = GoodwinOscillator(n=9.5)

sol = odeint(obj.dynamics, y0, t)
x_max = argrelextrema(sol[:,0], np.greater)
period = np.diff(t[x_max[-9:]]).mean()

sol_mask = sol[-(int(4*period/dt)):,:]
sol = sol_mask/sol_mask.mean(axis=0)
t = t[0:(int(4*period/dt))]

# PLOT TIMESERIES AND PHASE PORTRAIT
fig6 = plt.figure(figsize=(11,7))
ax2 = fig6.add_subplot(233) # timeseries
ax3 = fig6.add_subplot(234, projection='3d') # phase portrait
ax2.set_title('C', loc='left')

labels = ['$x$', '$y$', '$z$']
for i in range(len(y0)):
    ax2.plot(t, sol[:,i], label=labels[i])
ax3.plot(sol[:,0], sol[:,1], sol[:,2], c='k')
ax2.legend(framealpha=0, bbox_to_anchor=(1.0, 0.5), loc='center left')
ax2.set_xlabel('time $t$ (h)'); ax2.set_ylabel('rel. conc.')
ax2.set_xticks([0,24,48,72,96])
ax3.set_xlabel('$x$'); ax3.set_ylabel('$y$'); ax3.set_zlabel('$z$')
print(period)

ax2.set_aspect(0.75/ax2.get_data_ratio(), adjustable='box')
ax3.set_facecolor((0,0,0,0))

# PLOT BIFURCATION DIAGRAM: load XPPAUTO results of bifurcation 
data = pd.read_csv("./results/goodwin_k2_allinfo.dat", 
                   sep=" ", header=None)
print(data.keys())

ax3 = fig6.add_subplot(235)
ax4 = fig6.add_subplot(236)
ax3.set_title('E', loc='left')
ax4.set_title('F', loc='left')

subset = data.loc[data[0]==1]
ax3.plot(subset[3], subset[6], c="k")
subset2 = data.loc[data[0]==2]
ax3.plot(subset2[3].values[70:(-100)], subset2[6].values[70:(-100)], 
         c="white", lw=6)
subset3 = data.loc[data[0]==3]
ax3.plot(subset3[3], subset3[6], color='k')
ax3.plot(subset3[3], subset3[9], color='k')
idx = np.where(subset3[3].values[1:] < obj.k2)[0][0]
ax3.plot(subset3[3].values[idx], subset3[6].values[idx], marker='*',
         color="crimson", markersize=15)
ax3.plot(subset3[3].values[idx], subset3[9].values[idx], marker='*', 
         color="crimson", markersize=15)
ax3.set_xlim([0, 0.4]); ax3.set_ylim([0, 0.085])
ax3.set_aspect(0.75/ax3.get_data_ratio(), adjustable='box')
ax3.set_ylabel("$x$ min, max")
ax3.set_xlabel("$x$ degradation rate (1/h)")

ax4.plot(subset3[3],  subset3[5], color='k')
ax4.plot(subset3[3].values[idx], subset3[5].values[idx], marker='*',
         color="crimson", markersize=15)
ax4.set_xlabel("$x$ degradation rate (1/h)")
ax4.set_ylabel(r"period (h)")
ax4.set_xlim([0, 0.4]); ax4.set_ylim([20, 40])
ax4.set_aspect(1/ax4.get_data_ratio(), adjustable='box')

ax0 = fig6.add_subplot(231)
ax00 = fig6.add_subplot(234)
ax0.set_axis_off(); ax00.set_axis_off()
ax0.set_title('A', loc='left'); ax00.set_title('D', loc='left')
ax0.set_aspect(1/ax0.get_data_ratio(), adjustable='box')
ax00.set_aspect(1/ax00.get_data_ratio(), adjustable='box')


# PLOT TIMESERIES WITH "BROKEN" AXIS, UNNORMALIZED SOLUTION (PANEL B)
y0 = [.1, .1, .1]
dt = .01
t = np.arange(0, 10000, dt)

obj = GoodwinOscillator(n=9.5)

sol = odeint(obj.dynamics, y0, t)
x_max = argrelextrema(sol[:,0], np.greater)
period = np.diff(t[x_max[-9:]]).mean()

sol_mask = sol[-(int(4*period/dt)):,:]
sol = sol_mask
t = t[0:(int(4*period/dt))]

ax_bot = fig6.add_subplot(435)
ax_top = fig6.add_subplot(432, sharex=ax_bot)
ax_top.set_title('B', loc='left')

ax_bot.plot(t, sol[:,0], label='$x$')
ax_bot.plot(t, sol[:,1], label='$y$')
ax_top.plot(t, sol[:,2], label='$z$', c='#2ca02c')
ax_bot.set_xlabel('time $t$ (h)'); ax_bot.set_ylabel('conc.')
ax_bot.set_xticks([0,24,48,72,96])

# hide the spines between ax and ax2
ax_top.spines['bottom'].set_visible(False)
ax_bot.spines['top'].set_visible(False)
ax_top.spines['top'].set_visible(False)
ax_top.spines['right'].set_visible(False)
ax_bot.spines['right'].set_visible(False)
ax_top.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
ax_bot.xaxis.tick_bottom()

d = .5  # proportion of vertical to horizontal extent of the slanted line
kwargs = dict(marker=[(-1, -d), (1, d)], markersize=12,
              linestyle="none", color='k', mec='k', mew=1, clip_on=False)
ax_top.plot([0, 1], [0, 0], transform=ax_top.transAxes, **kwargs)
ax_bot.plot([0, 1], [1, 1], transform=ax_bot.transAxes, **kwargs)

ax_bot.set_aspect((0.75/2)/ax_bot.get_data_ratio(), adjustable='box')
ax_top.set_aspect((0.75/2)/ax_top.get_data_ratio(), adjustable='box')

fig6.subplots_adjust(
top=0.88,
bottom=0.11,
left=0.005,
right=0.925,
hspace=0.8,
wspace=0.3,
)

fig6.patch.set_alpha(0)
fig6.show()
isExist = os.path.exists('./figures/')
if not isExist:  
    os.makedirs('./figures/')
#fig6.savefig('./figures/fig6bis.pdf', format='pdf')

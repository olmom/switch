"""
tetramer formation and multisite phosphosylation
can lead to ultrasensitivity that can be quantified
with the Hill function
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint
from utils.multisite_phosphorylation import MultisitePhosphorylation
from utils.fits import Fits
from scipy.optimize import curve_fit
import os
plt.rcParams['text.usetex'] = True
plt.rcParams['font.size'] = '13'
plt.rcParams['legend.fontsize'] = '12'
plt.rcParams['xtick.labelsize'] = '12'
plt.rcParams['ytick.labelsize'] = '12'
plt.rcParams['axes.spines.top'] = False
plt.rcParams['axes.spines.right'] = False


# FIGURE 2: ULTRASENSITIVITY THROUGH COOPERATIVE MULTISITE 
# PHOSPHORYLATION (OR BINDING OF O2 TO HEMOGLOBIN)
t = np.arange(0, 10000, 1)
y0 = [0.1, 0.1, 0.1, 0.1, 0.1]

colors = ['#1B9E77', '#D95F02', '#7570B3']
fig2 = plt.figure(figsize=(11,4))
ax = fig2.add_subplot(122)

kinase = np.arange(0,5,0.01)

# Plot how the steady state of xPPPP changes with increasing kinase
# concentrations, for different values of k4
k4_values = np.array([1, 10, 100])
responses = []
# iterate through different k4 values
for k in range(len(k4_values)):
    obj = MultisitePhosphorylation(
        k1=1, k2=1, k3=1, k4=k4_values[k],
        k1r=1, k2r=1, k3r=1, k4r=1, pase=1,
        kinase=kinase 
    )  
    sol = odeint(obj.dynamics_4P, np.tile(y0, len(kinase)), t)
    xPPPP = sol[:,4::5]
    xTot = sol[:,0::5] + sol[:,1::5] + sol[:,2::5] + \
        sol[:,3::5] + sol[:,4::5]
    rel_xPPPP = xPPPP[-1]/xTot[-1]

    relation = obj.k4 / obj.k1
    ax.plot(kinase, rel_xPPPP, c=colors[k],
            label='$k_4={}k_1$'.format(format(relation, '.0f')))

    # fit Hill curve and calculate n
    x_data, y_data = kinase, rel_xPPPP
    popt, pcov = curve_fit(Fits.hill_curve, x_data, y_data, p0=(1,1,1))
    print('Hill coeff. of fit for k4={}k1 => {}'.format(
        int(k4_values[k]), format(popt[1], '.2f')))

ax.set_xlabel('concentration of kinase (or O$_2$)'); ax.set_ylabel('fraction of xPPPP (or HbO$_2$)')
ax.set_ylim([0,1.05]); ax.set_xlim([0,5])
ax.set_yticks([0,0.5,1]); ax.set_xticklabels([])
ax.legend(framealpha=0)
ax.set_aspect(1/ax.get_data_ratio(), adjustable='box')
ax.set_title('B', loc='left')

ax0 = fig2.add_subplot(121)
ax0.axis('off')
ax0.set_title('A', loc='left')
ax0.set_aspect(1/ax0.get_data_ratio(), adjustable='box')

fig2.subplots_adjust(
    top=0.88,
    bottom=0.14,
    left=0.11,
    right=0.9,
    hspace=0.2,
    wspace=0.2
)

isExist = os.path.exists('./figures/')
if not isExist:  
    os.makedirs('./figures/')
#fig2.savefig('./figures/fig2.pdf', format='pdf')



import numpy as np
from scipy.interpolate import splrep 
from scipy.interpolate import splev
from scipy.integrate import odeint

class DDE:

    def __init__ (self, tau=10, beta=4, n=4, d=0.1, K=1):
                  #k=0.1, K=1, n=2, f=0, d=.1, stimulus=.1, 
                  #A_tot=1):
        """
        Parameters for different signal-response curves
        """
        self.tau = tau
        self.beta = beta
        self.n = n
        self.d = d
        self.K = K

        self.dimension = 1


    def production_term (self, y, t):
        xd = y
        prod = self.beta / (self.K + (xd)**self.n) 
        return prod

    def degradation_term (self, y, t):
        x = y
        degr = self.d*x 
        return degr        

    def delay_model1 (self, y, t):
        """
        http://be150.caltech.edu/2019/handouts/09_delay_oscillators.html
        """
        x = y(t)
        xd = y(t-self.tau)
        dxdt = self.beta / (self.K + (xd)**self.n) - self.d*x
        return dxdt      

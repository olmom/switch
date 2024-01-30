import numpy as np

class SignalResponseCurve:

    def __init__ (self, k=0.1, K=1, n=2, f=0, d=.1, stimulus=.1, 
                  A_tot=1):
        """
        Parameters for different signal-response curves
        """
        self.k = k
        self.d = d
        self.K = K
        self.n = n
        self.f = f 
        self.stimulus = stimulus
        self.A_tot = A_tot 

        self.dimension = 1

    def linear (self, input):
        output = self.k * input
        return output

    def michaelis_menten (self, input):
        output = self.k*input / (self.K+input)
        return output

    def sigmoidal (self, input):
        output = self.k*input**self.n / (self.K**self.n+input**self.n)
        return output

    def bistable (self, output):
        A_star = output
        stimulus = (self.f * A_star**self.n * self.A_tot -\
                   self.d * self.K**self.n * A_star -\
                   (self.f + self.d) * A_star**(self.n+1)) /\
                   ((A_star - self.A_tot) * \
                    (A_star**self.n + self.K**self.n))
        return stimulus 
    
    def production_term (self, y0):
        A_star = y0
        prod = self.stimulus*(self.A_tot - A_star)
        return prod

    def PFL_term (self, y0):
        A_star = y0
        PFL  = (self.A_tot - A_star) * self.f * (A_star**self.n) / \
            (self.K**self.n + A_star**self.n)
        return PFL 

    def degradation_term (self, y0):
        A_star = y0
        degr = self.d * A_star
        return degr 

    def bistable_ODE (self, y0, t):
        """
        Equation from Xiong and Ferrell, 2003
        """
        prod = self.production_term(y0)
        PFL  = self.PFL_term(y0)
        degr = self.degradation_term(y0)
        dydt = prod + PFL - degr
        return dydt
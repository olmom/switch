import numpy as np


class GoodwinOscillator:

    def __init__ (self, k1=1.0, k2=0.20, k3=1.0, k4=0.15, k5=1.0, k6=0.10, 
                  K1=1, n=10, 
                  ):
        """
        Parameters of a Goodwin oscillator with linear degradation
        Reference: Goodwin BC. Oscillatory behavior in enzymatic control 
        processes
        Advances in enzyme regulation. 1965;3:425-37.
        Parameters from:
        https://journals.sagepub.com/doi/pdf/10.1177/074873099129001037
        https://www.sciencedirect.com/science/article/pii/S0022519396900673?via%3Dihub
        but modified because default parameters produced damped rhythms
        """
        self.k1 = k1
        self.k2 = k2
        self.k3 = k3
        self.k4 = k4
        self.k5 = k5
        self.k6 = k6
        self.K1 = K1
        self.n = n
        self.dimension = 3

    def dynamics (self, y0, t): 
        """
        Goodwin oscillator with linear degradation
        Ref: Goodwin BC. Oscillatory behavior in enzymatic control processes. 
        Advances in enzyme regulation. 1965;3:425-37.
        """ 
        f = np.zeros_like(y0)
        x = y0[0::3]
        y = y0[1::3]
        z = y0[2::3]

        fx = f[0::3]
        fy = f[1::3]
        fz = f[2::3]

        # odes
        fx[:] = self.k1*((self.K1**self.n)/(self.K1**self.n + z**self.n)) - \
                self.k2*x 
        fy[:] = self.k3*x - self.k4*y
        fz[:] = self.k5*y - self.k6*z
        
        return f       

import numpy as np

class MultisitePhosphorylation:

    def __init__ (self, k1=.1, k1r=.1, k2=.1, k2r=.2, 
                  k3=.1, k3r=.2, k4=.1, k4r=.2,
                  kinase=.1, pase=.1):
        """
        Parameters for different signal-response curves
        """
        self.k1  = k1 
        self.k1r = k1r
        self.k2  = k2 
        self.k2r = k2r
        self.k3  = k3 
        self.k3r = k3r    
        self.k4  = k4 
        self.k4r = k4r
        self.kinase = kinase
        self.pase = pase

        self.dimension = 3

    def dynamics (self, y0, t):
        """
        Ferrell Ultrasensitivity Part 2
        """
        f = np.zeros_like(y0) 
        x, xP, xPP = y0[0::3], y0[1::3], y0[2::3]

        fx, fxP, fxPP   = f[0::3], f[1::3], f[2::3]
        fx[:] = -self.k1*self.kinase*x + self.k1r*self.pase*xP
        fxP[:] = self.k1*self.kinase*x - self.k1r*self.pase*xP -\
                 self.k2*self.kinase*xP + self.k2r*self.pase*xPP
        fxPP[:] = self.k2*self.kinase*xP - self.k2r*self.pase*xPP 
        return f

    def dynamics_3P (self, y0, t):
        """
        Ferrell Ultrasensitivity Part 2
        """
        f = np.zeros_like(y0) 
        x, xP, xPP, xPPP = y0[0::4], y0[1::4], y0[2::4], y0[3::4]

        fx, fxP, fxPP, fxPPP   = f[0::4], f[1::4], f[2::4], f[3::4]

        fx[:] = -self.k1*self.kinase*x + self.k1r*self.pase*xP
        fxP[:] = self.k1*self.kinase*x - self.k1r*self.pase*xP -\
                 self.k2*self.kinase*xP + self.k2r*self.pase*xPP
        fxPP[:] = self.k2*self.kinase*xP - self.k2r*self.pase*xPP -\
                  self.k3*self.kinase*xPP + self.k3r*self.pase*xPPP
        fxPPP[:] = self.k3*self.kinase*xPP - self.k3r*self.pase*xPPP 
        return f         

    def dynamics_4P (self, y0, t):
        """
        Ferrell Ultrasensitivity Part 2
        """
        f = np.zeros_like(y0) 
        x, xP, xPP, xPPP, xPPPP = y0[0::5], y0[1::5], y0[2::5], \
            y0[3::5], y0[4::5]

        fx, fxP, fxPP, fxPPP, fxPPPP = f[0::5], f[1::5], f[2::5], \
            f[3::5], f[4::5]

        fx[:] = -self.k1*self.kinase*x + self.k1r*self.pase*xP
        fxP[:] = self.k1*self.kinase*x - self.k1r*self.pase*xP -\
                 self.k2*self.kinase*xP + self.k2r*self.pase*xPP
        fxPP[:] = self.k2*self.kinase*xP - self.k2r*self.pase*xPP -\
                  self.k3*self.kinase*xPP + self.k3r*self.pase*xPPP
        fxPPP[:] = self.k3*self.kinase*xPP - self.k3r*self.pase*xPPP -\
                   self.k4*self.kinase*xPPP + self.k4r*self.pase*xPPPP 
        fxPPPP[:] = self.k4*self.kinase*xPPP - self.k4r*self.pase*xPPPP  
        return f    

class CascadeAmplification:

    def __init__ (self, K=.5, n=1, 
                  ):
        """
        Parameters for different signal-response curves
        """
        self.K  = K 
        self.n = n
        self.dimension = 1

    def input_output (self, input):
        return input**self.n / (self.K**self.n + input**self.n)




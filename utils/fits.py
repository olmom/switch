import numpy as np

class Fits:

    def hill_curve (x, k, n, K):
        y = k * x**n / (K**n + x**n)
        return y 
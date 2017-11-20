# -*- coding: utf-8 -*-
"""
Created on Sun Aug 27 19:35:34 2017

@author: dwi
"""

import numpy as np
from numpy import cos, pi, exp
import matplotlib.pyplot as pl

t = np.linspace (0, 50e-6, 300)
x0 = np.array([1, 0.9, 0.8, 0.7, 0.6, 5, 20e-6, 0.0]) #Initial coefficient

def func (t, A, B, C, D, E, tR, T2, c):
    return A * exp(-((t)/T2)**4) + B * exp(-((t-tR)/T2)**4) + C * exp(-((t-2*tR)/T2)**4) + D * exp(-((t-3*tR)/T2)**4) + E * exp(-((t-4*tR)/T2)**4) + c

print(x0)
pl.plot(t*1e6, func(t, *x0))
pl.show()
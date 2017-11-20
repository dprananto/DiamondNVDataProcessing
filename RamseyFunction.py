# -*- coding: utf-8 -*-
"""
Created on Tue Aug  1 16:12:35 2017

@author: dwi
"""
import numpy as np
from numpy import cos, pi, exp
import matplotlib.pyplot as pl

d =10e6
f = 2.2e6
f1 = d - f
f2 = d
f3 = d + f
T2s = 0.2e-6
p1 = 2e6
p2 = 2e6
p3 = 2e6
x0 = np.array([f1, f2, f3, p1, p2, p3, T2s, 1]) #Initial coefficient

t = np.linspace (0, 1e-6, 300)
def func (t, f1, f2, f3, p1, p2, p3, T, c):
    return (cos(2*pi*f1*t + p1) + cos(2*pi*f2*t + p2) + cos(2*pi*f3*t + p3)) * exp(-(t/T)**1) + c

print(f1, f2, f3)
pl.plot(t*1e6, func(t, *x0))
pl.show()

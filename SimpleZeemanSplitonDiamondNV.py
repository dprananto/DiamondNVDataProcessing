# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 21:34:37 2017

@author: Dwi Prananto
"""

from physicsconstants import h_bar, G_e, k_B, u0, g_e, u_B
import numpy as np
from numpy import exp, cos, sin, pi, sqrt, mat, dot, linalg
import matplotlib.pyplot as pl

t = 60.0 * pi /180.0 #\theta
p = 00.0 * pi / 180.0 #\varphi
D = 2.873e9 #Zero-field splitting
B = np.linspace(0.0, 100e-3)
I = 1j
f_p = D + G_e * B * cos(t)
f_m = D - G_e * B * cos(t)

print(I**2)
pl.plot (B, f_p)
pl.plot (B, f_m)
pl.show()
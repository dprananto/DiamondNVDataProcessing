# -*- coding: utf-8 -*-
"""
Created on Tue Jan  3 17:42:52 2017

@author: dwi
"""

from physicsconstants import h_bar, G_e, k_B, u0, g_e, u_B
import numpy as np
from numpy import exp, cos, sin, pi, sqrt
import matplotlib.pyplot as pl


D = 2.873e9 #* 2*pi #zero-field splitting
B = 10.0e-3
r = 40e-9
S = 1
J = 0.5
w = -G_e * B
T_m = np.linspace(298, 350, 100) #magnon temperature
#phi = np.linspace(0, 180, 100) 
phi = 0 * pi /180.0
n = 1 / (exp(h_bar * w / k_B * T_m) - 1.0) #magnons population - Bose-Einstein
b_jk = - (u0 * g_e**2 * u_B**2) / (r**3.0 * h_bar * 4.0 * pi)#dipolar strength
y1 = - (B/3.0) * cos(phi) * n - 2 * n * b_jk * 1.0/3.0 * h_bar
#y2 = D - (B/3) * cos(phi) * n #+ (sin(phi))**2.0) * 0.5 - bjk * ((3.0 * (r**2.0 / 3.0)) - 1.0/3.0) * n * 0.5
#print (y2)
pl.plot(T_m, y1)#, T_m, y2)
print(b_jk)
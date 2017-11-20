# -*- coding: utf-8 -*-
"""
Created on Mon Jan  9 19:48:24 2017

@author: dwi
"""
#Calculation of magnons population based on Bose-Einstein statistics

from physicsconstants import h_bar, G_e, k_B
import numpy as np
from numpy import exp
import matplotlib.pyplot as pl

B = 10.0e-3
w = G_e * B

T_m = np.linspace(300.0, 360.0, 100)
n = 1 / (exp(h_bar * w / k_B * T_m) - 1.0)

pl.xlabel('Magnon temperature (K)')
pl.ylabel('Magnon population')
pl.plot (T_m, n)

print (n)
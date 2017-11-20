# -*- coding: utf-8 -*-
"""
Created on Fri Jan 20 18:04:17 2017

@author: dwi
"""

######## CURVE FIT ##########################
import numpy as np
from numpy import exp
from scipy.optimize import curve_fit
import matplotlib.pyplot as pl
import pandas as pd
from physicsconstants import h_bar, G_e, k_B

path = r'/home/dwi/ownCloud/Anlabshared/Dwi/researchdata/'
filename = r'Calibration_Big_Electromagnet_for_SSE-NVC_Setup'
Data = pd.read_csv (path + filename + r'.csv')
x = np.array (Data['V[V]'])
y = np.array (Data['B[Oe]'])
B = 10.0e-3
w = G_e * B
#Initial guess
a = 1
b = 1

x0 = np.array([a, b])

def func (x, a, b):
    return a + b * x
    
    
fit = curve_fit (func, x, y, x0)
f = fit[0]
print (a)
X = np.linspace(0, 4, 100)
fitfunc = f[0] + f[1] * X 

pl.xlabel (r'$V$ [V]', fontsize = 20.0)
pl.ylabel ('$H$ [Oe]', fontsize = 20.0)
pl.rc('text', usetex=True)
pl.rc('font', family='serif', size = 20)
pl.plot (X, fitfunc)
pl.plot (x, y, 'ro', markersize = 10.0)#, label = r'$\downarrow$')
print (f[0], f[1])

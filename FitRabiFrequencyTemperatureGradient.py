# -*- coding: utf-8 -*-
"""
Created on Mon Sep 11 16:49:20 2017

@author: dwi
"""

######## CURVE FIT ##########################
import numpy as np
from numpy import exp, cos, sin, pi
from scipy.optimize import curve_fit
import matplotlib.pyplot as pl
import pandas as pd
#from physicsconstants import h_bar, G_e, k_B
path = r'/home/dwi/ownCloud/researchdata/171112SSE/'
filename = r'SummaryRabiFrequencyH33V0dBmTdia310K'
Data = pd.read_csv ( path + filename + '.csv')
x = np.array (Data['DT'])
y = np.array (Data['f'])
#Initial guess

a = 1.0
b = 1.0
x0 = np.array([a, b])

def func (x, a, b):
    return a + b * x

popt, pcov = curve_fit (func, x, y, p0=x0)    
    
DT = np.linspace(-15, 15, 100)
#fitfunc = a[0] * exp(-a[1]*T)

pl.xlabel (r'$\Delta T$ (K)', fontsize = 20.0)
pl.ylabel ('$\Omega_R/ 2\pi$ (MHz)', fontsize = 20.0)
pl.rc('text', usetex=True)
pl.rc('font', family='serif', size = 20)
#pl.ticklabel_format(style='sci', axis='x', scilimits=(00,00))
#pl.xlim(-16, 16)
#pl.ylim(2.7, 3.6)
#pl.plot (T, func(T, *popt2), lw = 2.0, color = 'g', label = r'2.8 GHz 37 dBm 100 ns, T$_1$ ='+ str(round((1/popt2[1]) * 1e9, 3)) + 'ns')
#pl.plot (T, func(T, *popt3), lw = 2.0, color = 'r', label = r'4 GHz 37 dBm 10 ns, T$_1$ ='+ str(round((1/popt3[1]), 20)) + 's')
pl.plot (x, y, 'ko',markersize = 15.0, markerfacecolor='None', markeredgewidth=2 ,markeredgecolor='k')#, label = r'No MW')
pl.plot (DT, func(DT, *popt), lw = 3, color = 'r')
#pl.plot (t, I2, 'go', markersize = 10.0)#, label = r'3 GHz 32 dBm')
#pl.plot (t, I3, 'ro', markersize = 10.0)#, label = r'$\Delta$T = 30 K')
#pl.legend(fontsize = 14)
pl.savefig(path + filename + '.png', format = 'png', bbox_inches= 'tight')
print('a =' + str(round(popt[0],2)) + 'MHz')
print('b =' + str(round(popt[1], 9)) + 'MHz/K')

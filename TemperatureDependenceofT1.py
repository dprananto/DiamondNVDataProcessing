# -*- coding: utf-8 -*-
"""
Created on Sat Apr 15 19:25:02 2017

@author: dwi
Temperature Dependence of T_1 based on Jarmola et al. PRL 108, 197601 (2012)
"""
import matplotlib.pyplot as pl
import numpy as np
from numpy import exp
from physicsconstants import k_B
path = r'/home/dwi/ownCloud/Anlabshared/Dwi/researchdata/170613SpinSeebeck110Diamond/'
T = np.linspace (298, 310, 100) #Temperature
A1 = 21
A2 = 2.1e3
D = 73
A3 = 2.2e-11
G = A1 + (A2 /(exp(D/(k_B*T)) - 1)) + A3*T**5 #Relaxation rate as a function of temperature
T1 = 1/G #Relaxation time
pl.plot (T, T1 * 1e3, lw = 3, color = 'r')
#pl.plot (T, G, lw = 2.5, color = 'r') #Plot relaxation rate
pl.xlabel(r'T (K)')
pl.ylabel(r'T$_1$ (ms)')
#pl.ticklabel_format(style='sci', axis='y', scilimits=(3,4))
pl.savefig(path + r'T1vsT_PhononProcess.png', format = 'png', bbox_inches = 'tight', frameon = False)
pl.show
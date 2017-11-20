# -*- coding: utf-8 -*-
"""
Created on Sat Apr  8 14:43:27 2017

@author: dwi
"""

######## CURVE FIT ##########################
import numpy as np
from numpy import exp
from scipy.optimize import curve_fit
import matplotlib.pyplot as pl
import pandas as pd
#from physicsconstants import h_bar, G_e, k_B

Data = pd.read_csv (r'/home/dwi/ownCloud/Anlabshared/Dwi/researchdata/170425YIGFilmNanodiamond/T1ND2PulsedMW100ns/200ns2.8GHz-27dBm/200ns2.8GHz-27dBm.csv', header=None)
Data.columns = ['Intensity']
Data2 = pd.read_csv (r'/home/dwi/ownCloud/Anlabshared/Dwi/researchdata/170425YIGFilmNanodiamond/T1ND2PulsedMW100ns/200ns2.8GHz10dBm/200ns2.8GHz10dBm.csv', header=None)
Data2.columns = ['Intensity']
#Data3 = pd.read_csv (r'/home/dwi/ownCloud/Anlabshared/Dwi/researchdata/170424YIGFilm//Nanodiamond/T1/200ns4GHz10dBmPulsed10ns/200ns4GHz10dBmPulsed10ns.csv', header=None)
#Data3.columns = ['Intensity']
t = np.linspace(0.0, 200.0e-9,11) #Delay time
I = np.array (Data['Intensity']/Data['Intensity'].max()) #Intensity
I2 = np.array (Data2['Intensity']/Data2['Intensity'].max())
#I3 = np.array (Data3['Intensity']/Data3['Intensity'].max())
#B = 10.0e-3
#w = G_e * B
#Initial guess
x0 = np.array([1.0, 200.0, 2.0]) #Initial coefficient

def func (t, A, R, c):
    return A * exp(-R * t) + c
    
    
popt, pcov = curve_fit (func, t, I, p0=x0)
popt2, pcov2 = curve_fit (func, t, I2, p0=x0)
#popt3, pcov3 = curve_fit (func, t, I3, p0=x0)
#a = fit[0]
#print (a)
T = np.linspace(0.0, 200.0e-9, 100)
#fitfunc = a[0] * exp(-a[1]*T)

pl.xlabel (r'$\tau$ (s)', fontsize = 20.0)
pl.ylabel ('Intensity (norm.)', fontsize = 20.0)
pl.rc('text', usetex=True)
pl.rc('font', family='serif', size = 20)
pl.ticklabel_format(style='sci', axis='x', scilimits=(0,00))
pl.plot (T, func(T, *popt), lw = 2.0, color = 'b', label = r'2.8 GHz 0 dBm, T$_1$ ='+ str(round((1/popt[1]) * 1e9, 3)) + 'ns')
pl.plot (T, func(T, *popt2), lw = 2.0, color = 'g', label = r'2.8 GHz 37 dBm 100 ns, T$_1$ ='+ str(round((1/popt2[1]) * 1e9, 3)) + 'ns')
#pl.plot (T, func(T, *popt3), lw = 2.0, color = 'r', label = r'4 GHz 37 dBm 10 ns, T$_1$ ='+ str(round((1/popt3[1]), 20)) + 's')
pl.plot (t, I, 'bo', markersize = 10.0)#, label = r'No MW')
pl.plot (t, I2, 'go', markersize = 10.0)#, label = r'3 GHz 32 dBm')
#pl.plot (t, I3, 'ro', markersize = 10.0)#, label = r'$\Delta$T = 30 K')
pl.legend(fontsize = 14)
print('Spin relaxation rate =' + str(round(popt[1],2)))
print('Spin relaxation time (T1) =' + str(round(1/popt[1], 5)) + 's')
print('A coefficient =' + str(popt[0]))
print('c coefficient =' + str(popt[2]))
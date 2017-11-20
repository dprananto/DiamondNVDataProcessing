# -*- coding: utf-8 -*-
"""
Created on Mon Sep  4 14:42:39 2017

@author: dwi
"""

######## CURVE FIT ##########################
import numpy as np
from numpy import exp, cos, sin, pi
from scipy.optimize import curve_fit
import matplotlib.pyplot as pl
import pandas as pd
#from physicsconstants import h_bar, G_e, k_B
path = r'/home/dwi/ownCloud/researchdata/171112SSE/171109SpinSeebeckonDiamond110/Rabi/'
filename = r'Rabi_FMR_0to1us_DT-12K_H33V_2708MHz_10dBm_atten'
Data = pd.read_csv (path + filename + r'.csv', header=None)
Data.columns = ['Intensity']
#Data2 = pd.read_csv (r'/home/dwi/ownCloud/Anlabshared/Dwi/researchdata/170425YIGFilmNanodiamond/T1ND2PulsedMW100ns/200ns2.8GHz10dBm/200ns2.8GHz10dBm.csv', header=None)
#Data2.columns = ['Intensity']
#Data3 = pd.read_csv (r'/home/dwi/ownCloud/Anlabshared/Dwi/researchdata/170424YIGFilm//Nanodiamond/T1/200ns4GHz10dBmPulsed10ns/200ns4GHz10dBmPulsed10ns.csv', header=None)
#Data3.columns = ['Intensity']
t = np.linspace(0, 1e-6, 51) #Delay time
I = np.array (Data['Intensity']/Data['Intensity'].max()) #Intensity
#I = I[0:64]
#I2 = np.array (Data2['Intensity']/Data2['Intensity'].max())
#I3 = np.array (Data3['Intensity']/Data3['Intensity'].max())
#B = 10.0e-3
#w = G_e * B

#Initial guesses
f = 3e6
p = 1
T2s = 0.3e-6
n = 1
d = 0.0000

x0 = np.array([1, f, p, T2s, n, 1, d]) #Initial coefficient

def func (t, A, f, p, T2s, n, c, d):
    return A * (0.5 * (cos(2*pi*f*t + p))**1) * exp(-(t/T2s)**n) + c + d*t
    
    
popt, pcov = curve_fit (func, t, I, p0=x0)
#popt2, pcov2 = curve_fit (func, t, I2, p0=x0)
#popt3, pcov3 = curve_fit (func, t, I3, p0=x0)
#a = fit[0]
#print (a)
T = np.linspace(0, 1e-6, 300)
#fitfunc = a[0] * exp(-a[1]*T)

pl.xlabel (r'MW duration ($\mu$s)', fontsize = 20.0)
pl.ylabel ('Intensity (norm.)', fontsize = 20.0)
pl.rc('text', usetex=True)
pl.rc('font', family='serif', size = 20)
pl.ticklabel_format(style='sci', axis='x', scilimits=(00,00))
#pl.ylim(0.958, 1)
#pl.plot (T, func(T, *popt2), lw = 2.0, color = 'g', label = r'2.8 GHz 37 dBm 100 ns, T$_1$ ='+ str(round((1/popt2[1]) * 1e9, 3)) + 'ns')
#pl.plot (T, func(T, *popt3), lw = 2.0, color = 'r', label = r'4 GHz 37 dBm 10 ns, T$_1$ ='+ str(round((1/popt3[1]), 20)) + 's')
pl.plot (t, I, 'o', markersize = 12.0, markerfacecolor='None', markeredgewidth=3 ,markeredgecolor='r',)#, label = r'No MW')
pl.plot (T, func(T, *popt), lw = 4, color = 'k', label = r'$\Delta T = -12$ K, $\Omega_R/2 \pi$ =' + str(round(popt[1] * 1e-6, 3) ) + r'MHz')
#pl.plot (t, I2, 'go', markersize = 10.0)#, label = r'3 GHz 32 dBm')
#pl.plot (t, I3, 'ro', markersize = 10.0)#, label = r'$\Delta$T = 30 K')
pl.legend(fontsize = 15)
pl.savefig(path + filename + r'.png', format = 'png', bbox_inches= 'tight')
#print('Rabi frequency =' + str(round(popt[1],2)))
print('Decoherence time (T2*) =' + str(round(popt[3], 12)) + 's')
#print('A coefficient =' + str(popt[0]))
print('n =' + str(popt[4]))
#print('c coefficient =' + str(popt[3]))
#print('p coefficient =' + str(popt[2]))
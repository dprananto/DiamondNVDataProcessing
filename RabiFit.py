# -*- coding: utf-8 -*-
"""
Created on Tue May 16 19:51:56 2017

@author: dwi
"""

######## CURVE FIT ##########################
import numpy as np
from numpy import exp, cos, sin, pi
from scipy.optimize import curve_fit
import matplotlib.pyplot as pl
import pandas as pd
#from physicsconstants import h_bar, G_e, k_B
path = r'/home/dwi/ownCloud/Anlabshared/Dwi/researchdata/170923SpinSeebeck110Diamond/Rabi/'
filename = r'Rabi_0to2us_2606MHz_10dBm_DT14.5K_Tdia310K_190G_distance70um_hyperfine_Anti'
Data = pd.read_csv (path + filename + r'.csv', header=None)
Data.columns = ['Intensity']
#Data2 = pd.read_csv (r'/home/dwi/ownCloud/Anlabshared/Dwi/researchdata/170425YIGFilmNanodiamond/T1ND2PulsedMW100ns/200ns2.8GHz10dBm/200ns2.8GHz10dBm.csv', header=None)
#Data2.columns = ['Intensity']
#Data3 = pd.read_csv (r'/home/dwi/ownCloud/Anlabshared/Dwi/researchdata/170424YIGFilm//Nanodiamond/T1/200ns4GHz10dBmPulsed10ns/200ns4GHz10dBmPulsed10ns.csv', header=None)
#Data3.columns = ['Intensity']
t = np.linspace(0.0, 2e-6,51) #Delay time
I = np.array (Data['Intensity']/Data['Intensity'].max()) #Intensity
I = I[0:64]
#I2 = np.array (Data2['Intensity']/Data2['Intensity'].max())
#I3 = np.array (Data3['Intensity']/Data3['Intensity'].max())
#B = 10.0e-3
#w = G_e * B
#Initial guess
f1 = 0.5e6
f2 = 0.8e6
p = 1.0
T2s = 0.01e-6
n = 0.1

x0 = np.array([1, f1, f2, p, T2s, n, 0]) #Initial coefficient

def func (t, A, f1, f2, p, T2s, n, c):
    return A**2 * ((0.5 * (cos(2*pi*f1*t + p))**1) + (0.5 * (cos(2*pi*f2*t + p))**1)) * exp(-(t/T2s)**n) + c
    
    
popt, pcov = curve_fit (func, t, I, p0=x0) 
#popt2, pcov2 = curve_fit (func, t, I2, p0=x0)
#popt3, pcov3 = curve_fit (func, t, I3, p0=x0)
#a = fit[0]
#print (a)
T = np.linspace(0.0, 5e-6, 300)
#fitfunc = a[0] * exp(-a[1]*T)

pl.xlabel (r'MW duration ($\mu$s)', fontsize = 20.0)
pl.ylabel ('Intensity (norm.)', fontsize = 20.0)
pl.rc('text', usetex=True)
pl.rc('font', family='serif', size = 20)
pl.ticklabel_format(style='sci', axis='x', scilimits=(00,00))
#pl.ylim(0.960, 1.005)
#pl.plot (T, func(T, *popt2), lw = 2.0, color = 'g', label = r'2.8 GHz 37 dBm 100 ns, T$_1$ ='+ str(round((1/popt2[1]) * 1e9, 3)) + 'ns')
#pl.plot (T, func(T, *popt3), lw = 2.0, color = 'r', label = r'4 GHz 37 dBm 10 ns, T$_1$ ='+ str(round((1/popt3[1]), 20)) + 's')
pl.plot (t, I, 'bo', markersize = 12.0)#, label = r'No MW')
pl.plot (T, func(T, *popt), lw = 3, color = 'r', label = r'fit, $\Omega_{R1} = 2\pi \times$' + str(round(popt[1] * 1e-6, 3) ) + r'MHz, $\Omega_{R2} = 2\pi \times$' + str(round(popt[2] * 1e-6, 3) ) + 'MHz')
#pl.plot (t, I2, 'go', markersize = 10.0)#, label = r'3 GHz 32 dBm')
#pl.plot (t, I3, 'ro', markersize = 10.0)#, label = r'$\Delta$T = 30 K')
pl.legend(fontsize = 14)
pl.savefig(path + filename + r'.png', format = 'png', bbox_inches= 'tight')
#print('Rabi frequency =' + str(round(popt[1],2)))
#print('Decoherence time (T2*) =' + str(round(popt[4], 9)) + 's')
#print('A coefficient =' + str(popt[0]))
#print('n =' + str(popt[5]))
#print('c coefficient =' + str(popt[3]))
#print('p coefficient =' + str(popt[3]))
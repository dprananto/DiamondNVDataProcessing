# -*- coding: utf-8 -*-
"""
Created on Fri Aug 25 17:54:15 2017

@author: dwi
"""

######## CURVE FIT ##########################
import numpy as np
from numpy import exp, cos, sin, pi
from scipy.optimize import curve_fit
from scipy.signal import hilbert
from scipy.signal import find_peaks_cwt
import matplotlib.pyplot as pl
import pandas as pd
#from physicsconstants import h_bar, G_e, k_B
path = r'/home/dwi/ownCloud/Anlabshared/Dwi/researchdata/170923SpinSeebeck110Diamond/SpinEcho/'
filename = r'SpinEcho_0to50us_2621MHz_10dBm_DT-12K_Tdia310K_190G_distance70um_hyperfine'
Data = pd.read_csv (path + filename + r'.csv', header=None)
Data.columns = ['Intensity']
#Data2 = pd.read_csv (r'/home/dwi/ownCloud/Anlabshared/Dwi/researchdata/170425YIGFilmNanodiamond/T1ND2PulsedMW100ns/200ns2.8GHz10dBm/200ns2.8GHz10dBm.csv', header=None)
#Data2.columns = ['Intensity']
#Data3 = pd.read_csv (r'/home/dwi/ownCloud/Anlabshared/Dwi/researchdata/170424YIGFilm//Nanodiamond/T1/200ns4GHz10dBmPulsed10ns/200ns4GHz10dBmPulsed10ns.csv', header=None)
#Data3.columns = ['Intensity']
t = np.linspace(0.0, 50e-6, 51) 
I = np.array (Data['Intensity']/Data['Intensity'].max()) #Intensity
#I = I[0:80]
#I2 = np.array (Data2['Intensity']/Data2['Intensity'].max())
#I3 = np.array (Data3['Intensity']/Data3['Intensity'].max())
#B = 10.0e-3
#w = G_e * B
#Initial guess
#analytic_signal = hilbert(I)
#envelope = np.abs(analytic_signal)
#peakind = find_peaks_cwt(I, np.arange(0.1e-6, 100e-6, 100))
#peakdata = I[peakind]

A = 1
f = 90e3
p = 0.2
T2 = 20e-6
n = 1
c = 1

x0 = np.array([A, f, T2, n, c]) #Initial coefficient

def func (t, A, f, T2, n, c):
    return A * cos (2*pi*f*t + p) * exp(-(t/T2)**n) + c 
    
    
popt, pcov = curve_fit (func, t, I, p0=x0)
#popt2, pcov2 = curve_fit (func, t, I2, p0=x0)
#popt3, pcov3 = curve_fit (func, t, I3, p0=x0)
#a = fit[0]
#print (a)
T = np.linspace(0.0, 50e-6, 100)
#tpeak = np.linspace (0.0, 50e-6, 9)
#fitfunc = a[0] * exp(-a[1]*T)

pl.xlabel (r'$\tau$(s)', fontsize = 20.0)
pl.ylabel ('Intensity (norm.)', fontsize = 20.0)
pl.rc('text', usetex=True)
pl.rc('font', family='serif', size = 20)
pl.ticklabel_format(style='sci', axis='x', scilimits=(00,00))
pl.ylim(0.955, 1.01)
#pl.plot (t, envelope, lw = 3, color = 'r', label = r'envelop')
#pl.plot(tpeak, peakdata, 'ro', markersize = 12.0)
#pl.plot (T, func(T, *popt2), lw = 2.0, color = 'g', label = r'2.8 GHz 37 dBm 100 ns, T$_1$ ='+ str(round((1/popt2[1]) * 1e9, 3)) + 'ns')
#pl.plot (T, func(T, *popt3), lw = 2.0, color = 'r', label = r'4 GHz 37 dBm 10 ns, T$_1$ ='+ str(round((1/popt3[1]), 20)) + 's')
pl.plot (t, I, 'bo-', markersize = 12.0)#, label = r'No MW')
pl.plot (T, func(T, *popt), lw = 3, color = 'r', label = r'fit, $T_2$ =' + str(round(popt[2] * 1e6, 3)) + '$\mu$s' + r', $\omega/2\pi$ = ' + str(round(popt[1] * 1e-3, 3)) + r'kHz')
#pl.plot (t, I2, 'go', markersize = 10.0)#, label = r'3 GHz 32 dBm')
#pl.plot (t, I3, 'ro', markersize = 10.0)#, label = r'$\Delta$T = 30 K')
pl.legend(fontsize = 14)
pl.savefig(path + filename + r'.png', format = 'png', bbox_inches= 'tight')
print('Decoherence time (T_2) =' + str(round(popt[2], 7)) + 's')
print('n =' + str(round(popt[3], 7)))
print('frequency =' + str(round(popt[1], 7)) + 'Hz')
#print('n =' + str(round(popt[3], 7)) + 'Hz')
#print(popt[0], popt[1], popt[2], popt[3], popt[4], popt[5])
#print(np.arange(0.01e-6, 100e-6, 0.1e-6))
#print(t)
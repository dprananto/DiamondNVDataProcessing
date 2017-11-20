# -*- coding: utf-8 -*-
"""
Created on Sat Jul 22 23:18:06 2017

@author: dwi
"""

######## CURVE FIT ##########################
import numpy as np
from numpy import exp, cos, pi
from scipy.optimize import curve_fit
import matplotlib.pyplot as pl
import pandas as pd
from scipy import stats
#from physicsconstants import h_bar, G_e, k_B
path = r'/home/dwi/ownCloud/Anlabshared/Dwi/researchdata/170826SpinSeebeck110Diamond/Ramsey/'
filename = r'Ramsey_0to1us_2602MHz_15dBm_DT6K_180G'
Data = pd.read_csv (path + filename + r'.csv', header=None)
Data.columns = ['Intensity']
#Data2 = pd.read_csv (r'/home/dwi/ownCloud/Anlabshared/Dwi/researchdata/170425YIGFilmNanodiamond/T1ND2PulsedMW100ns/200ns2.8GHz10dBm/200ns2.8GHz10dBm.csv', header=None)
#Data2.columns = ['Intensity']
#Data3 = pd.read_csv (r'/home/dwi/ownCloud/Anlabshared/Dwi/researchdata/170424YIGFilm//Nanodiamond/T1/200ns4GHz10dBmPulsed10ns/200ns4GHz10dBmPulsed10ns.csv', header=None)
#Data3.columns = ['Intensity']
t = np.linspace(0.0, 1.0e-6,128) #Delay time
I = np.array (Data['Intensity']/Data['Intensity'].max()) #Intensity
#I = I[0:80]
#I2 = np.array (Data2['Intensity']/Data2['Intensity'].max())
#I3 = np.array (Data3['Intensity']/Data3['Intensity'].max())
#B = 10.0e-3
#w = G_e * B

d = 12e6 #Detuning
f = 2.2e6 #Hyperfine split
p1 = 1
p2 = 1
p3 = 1
T2s = 0.2e-6
n = 2

#Initial guess
x0 = np.array([5, d, p1, p2, p3, T2s, n, 0.001]) #Initial coefficient

def func (t, A, d, p1, p2, p3, T2s, n, c):
    return A * (cos(2*pi*(d - f)*t + p1) + cos(2*pi*(d)*t + p2) + cos(2*pi*(d + f)*t + p3)) * exp(-(t/T2s)**n) + c
    
    
popt, pcov = curve_fit (func, t, I, p0=x0)
#popt2, pcov2 = curve_fit (func, t, I2, p0=x0)
#popt3, pcov3 = curve_fit (func, t, I3, p0=x0)
#a = fit[0]
#print (a)
tau = np.linspace(0.0, 1.0e-6, 500)
KS = stats.chisquare (I, func(t, *popt)) #Kolmogorov-Smirnov test


pl.xlabel (r'$\tau$ ($\mu$s)', fontsize = 20.0)
pl.ylabel ('Intensity (norm.)', fontsize = 20.0)
pl.rc('text', usetex=True)
pl.rc('font', family='serif', size = 20)
pl.ticklabel_format(style='sci', axis='x', scilimits=(0,00))
#pl.ylim(0.985,1.005)
pl.plot (t*1e6, I, 'bo', markersize = 10.0)#, label = r'No MW')
pl.plot (tau*1e6, func(tau, *popt), lw = 3, color = 'r', label = r'$T_2^*$ =' + str(round(popt[5]* 1e6, 3)) + r"$\mu$s" + '\n' + r'$\delta = $' + str(round(popt[1] * 1e-6,2)) + 'MHz')
pl.legend(fontsize = 14)
pl.savefig(path + filename + r'.png', format = 'png', bbox_inches= 'tight')
print('phase1 =' + str(round(popt[2],3))) + 'rad'
print('phase2 =' + str(round(popt[3],3))) + 'rad'
print('phase3 =' + str(round(popt[4],3))) + 'rad'
print('f_d =' + str(round(popt[1] * 1e-6,2)) + 'MHz') #Detuning frequency
#print('Dephasing time (T2*) =' + str(round(popt[4], 10)) + 's')
print('n =' + str(popt[6]))
#print(popt)
#print(pcov)
print(KS)
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 26 11:12:15 2017

@author: dwi
"""

######## CURVE FIT ##########################
import numpy as np
from numpy import exp
from scipy.optimize import curve_fit
import matplotlib.pyplot as pl
import pandas as pd
#from physicsconstants import h_bar, G_e, k_B

path = r"/home/dwi/ownCloud/Anlabshared/Dwi/researchdata/170826SpinSeebeck110Diamond/T1/"
Data = pd.read_csv (r'/home/dwi/ownCloud/Anlabshared/Dwi/researchdata/170826SpinSeebeck110Diamond/T1/T1_0to10ms_DT10K_180G.csv', header=None)
Data2 = pd.read_csv (r'/home/dwi/ownCloud/Anlabshared/Dwi/researchdata/170826SpinSeebeck110Diamond/T1/T1withPi_0to10ms_2612MHz_15dBm_DT10K_180G.csv', header=None)
Data.columns = ['Intensity']
Data2.columns = ['Intensity']
#Data2 = pd.read_csv (r'/home/dwi/ownCloud/Anlabshared/Dwi/researchdata/170613SpinSeebeck110Diamond/T1/20K_300G_0.5to20ms_49Ave.csv', header=None)
#Data2.columns = ['Intensity']
#Data3 = pd.read_csv (r'C:\Users\Dwi Prananto\ownCloud\Anlabshared\Dwi\researchdata\170429Diamond111onYIGFilm\T1\20ms3GHz10dBm10msPulsePos1-3\20ms3GHz10dBm10msPulsePos1-3.csv', header=None)
#Data3.columns = ['Intensity']
t = np.linspace(0.0e-3, 10.0e-3,21) #Delay time
I = np.array (Data['Intensity']) #IntensityT1
I2 = np.array (Data2['Intensity']) #IntensityT1withPi
Isub = I - I2 #Subtraction
Isubnorm = Isub/Isub.max()
#B = 10.0e-3
#w = G_e * B
#Initial guess
x0 = np.array([1.0, 2.0e-3, 1.0]) #Initial coefficient

def func (t, A, T1, c):
    return A * exp(-t/T1) + c
    
    
popt, pcov = curve_fit (func, t, Isubnorm, p0=x0)
#popt2, pcov2 = curve_fit (func, t, I2, p0=x0)
#popt3, pcov3 = curve_fit (func, t, I3, p0=x0)
#a = fit[0]
#print (a)
T = np.linspace(0.0e-3, 10.0e-3, 100)
#fitfunc = a[0] * exp(-a[1]*T)

pl.xlabel (r'$\tau$ (ms)', fontsize = 20.0)
pl.ylabel ('Intensity (norm.)', fontsize = 20.0)
pl.rc('text', usetex=True)
pl.rc('font', family='serif', size = 20)
pl.ticklabel_format(style='plain', axis='x', scilimits=(0,00))
pl.xlim(0, 11)
#pl.ylim(0.91, 1.01)
#pl.plot (T*1000, func(T, *popt2), lw = 3.0, color = 'b', label = r'$\Delta$T = 20K, T$_1$ ='+ str(round((1/popt2[1]) * 1e3, 3)) + 'ms')
#pl.plot (T, func(T, *popt3), lw = 2.0, color = 'b', label = r'3 GHz 37 dBm Pos 1-3, T$_1$ ='+ str(round((1/popt3[1]) * 1e3, 3)) + 'ms')
pl.plot (t*1000, Isubnorm, 'o', markersize = 15.0, markerfacecolor='None', markeredgewidth=2 ,markeredgecolor='b')#, label = r'No MW')
pl.plot (T*1000, func(T, *popt), lw = 3.0, color = 'r', label = r'fit, $T_1$ ='+ str(round((popt[1] * 1e3), 3)) + 'ms')
#pl.plot (t*1000, I2, '^', markersize = 15.0, markerfacecolor='None', markeredgewidth=2 ,markeredgecolor='k')#, label = r'3 GHz 32 dBm')
#pl.plot (t, I3, 'bo', markersize = 10.0)#, label = r'$\Delta$T = 30 K')
pl.legend(fontsize = 14)
pl.savefig(path + r"T1Sub_0to10ms_2612MHz_15dBm_DT10K_180G.png", format = 'png', bbox_inches = 'tight', frameon = False)
print('Spin relaxation rate =' + str(round(1/popt[1],2)))
print('Spin relaxation time (T1) =' + str(round(popt[1], 5)) + 's')
print('A coefficient =' + str(popt[0]))
print('c coefficient =' + str(popt[2]))
#print(I)
#print(I2)
#print(Isubnorm)
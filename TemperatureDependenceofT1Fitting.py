# -*- coding: utf-8 -*-
"""
Created on Sun Apr 16 15:14:22 2017

@author: dwi
"""

######## CURVE FIT ##########################
import numpy as np
from numpy import exp
from scipy.optimize import curve_fit
import matplotlib.pyplot as pl
import pandas as pd
from physicsconstants import k_B

Data = pd.read_csv ('/home/dwi/ownCloud/Anlabshared/Dwi/researchdata/170416SpinSeebeck/T1/SummarizedT1DataFunctionofTemperature.csv')
T = np.array (Data['T']) #Temperature
T1 = np.array (Data['T1']) #T1
#Initial guess
x0 = np.array([21.0, 2000.0, 2e-11, 73.0]) #Initial coefficient

def func (T, A1, A2, A3, D):
    return 1/ (A1 + (A2 /(exp(D/(k_B*T)) - 1)) + A3*T**5)
    
    
popt, pcov = curve_fit (func, T, T1, p0=x0)
Temp = np.linspace(300.0, 320, 100)
#fitfunc = a[0] * exp(-a[1]*T)

pl.xlabel (r'Temperature (K)', fontsize = 22.0)
pl.ylabel (r'T$_1$ (ms)', fontsize = 22.0)
pl.rc('text', usetex=True)
pl.rc('font', family='serif', size = 22)
#pl.ticklabel_format(style='sci', axis='x', scilimits=(0,00))
pl.plot (Temp, func(Temp, *popt), lw = 2.5, color ='r', label = r'Data fit')
pl.plot (T, T1, 'bo', markersize = 10.0, label = r'Data')
pl.legend(fontsize = 20)
print(r'A1 =' + str(popt[0]) + r's$^{-1}$')
print(r'A2 =' + str(popt[1]) + r's$^{-1}$')
print(r'A3 =' + str(popt[2]) + r'K$^{-5}$s$^{-1}$')
print(r'D =' + str(popt[3]) + 'meV')
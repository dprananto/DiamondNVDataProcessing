# -*- coding: utf-8 -*-
"""
Created on Fri Sep 22 15:00:36 2017

@author: dwi
"""

import pandas as pd
import numpy as np
from numpy import pi
from scipy.optimize import curve_fit
import matplotlib.pyplot as pl

#fi = input ("File path to open: ")
path = '/home/dwi/ownCloud/Anlabshared/Dwi/researchdata/170826SpinSeebeck110Diamond/ESR/'
filename = r'ESR_2595to2610MHz_DT-15K_-20dBm_190G_distance30um.csv'
Data = pd.read_csv (path + filename, header=None)
#Data2 = pd.read_csv ("/home/dwi/ownCloud/Anlabshared/Dwi/researchdata/170502Diamond111onYIGFilm/ESRYIG3+2um/ESRYIG3+2um.csv", header=None)
#ttl = input (r"Graph title: " )
#pl.title (ttl)
Data.columns = ['Intensity']
#Data.columns = ['Intensity']
#Data2.columns = ['Intensity']
fstart = 2.595e9
fstop = 2.610e9
x = np.linspace(fstart, fstop, 201)
y = np.array(Data['Intensity']/Data['Intensity'].max())
offset = 1.0 - np.mean (y[0:30])
#y2 = np.array(Data2['Intensity']/Data2['Intensity'].max())
yo = y + offset
#offset2 = 1.0 - np.mean (y2[0:30])
#yo2 = y2 + offset2
W = 1.000e6
f1 = 2.589e9
f2 = 2.591e9
f3 = 2.593e9

x0 = np.array([W, f1, f2, f3])

def func (x, W, f1, f2, f3):
    num1 = (0.5*W)
    denum1 = pi *(x - f1)**2 + pi*(0.5*W)**2
    num2 = (0.5*W)
    denum2 = pi *(x - f2)**2 + pi*(0.5*W)**2
    num3 = (0.5*W)
    denum3 = pi *(x - f3)**2 + pi*(0.5*W)**2 
    return 1 - ((num1/denum1) + (num2/denum2) + (num3/denum3))
    
    
popt, pcov = curve_fit (func, x, yo, p0=x0)


X = np.linspace(fstart, fstop, 501)
pl.rc('text', usetex=True)
pl.rc('font', family='serif', size = 20)
#pl.plot (X, func(X, *popt), lw = 3.0, color = 'r')
pl.plot (x, yo, 'bo-', color = 'b')
#pl.plot (x, yo2, lw = 2.5, label = 'pos 1 + 2 $\mu$m', color = 'b')
pl.ylim(0.988, 1.002)
pl.xlabel(r"$f$ (GHz)")
pl.ylabel('Intensity (norm.)')

#pl.legend(loc=0, ncol=1, fontsize = 12)
pl.savefig(path + filename + r'.png', format = 'png', bbox_inches = 'tight', frameon = False)
pl.show()
#print('FWHM =' + str(round(popt[0],2)))
#print('Center frequency 1=' + str(round(popt[1],2)))
#print('Center frequency 2=' + str(round(popt[2],2)))
#print('Center frequency 2=' + str(round(popt[3],2)))
#print('contrast = '+ str((1.0 - np.min(yo))*100))
#print (1 - np.min(yo))
#print (x)
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 18 18:37:46 2017

@author: dwi
"""

###############################################################################
#THIS PROGRAM PLOT ESR SPECTRUM AND LORENTZIAN FIT ############################
###############################################################################
import pandas as pd
import numpy as np
from numpy import pi
from scipy.optimize import curve_fit
import matplotlib.pyplot as pl

#fi = input ("File path to open: ")
path = '/home/dwi/ownCloud/Anlabshared/Dwi/researchdata/170725SpinSeebeck110Diamond/ESR/'
Data = pd.read_csv (r'/home/dwi/ownCloud/Anlabshared/Dwi/researchdata/170725SpinSeebeck110Diamond/ESR/DT-10K/2017_07_30_DT-10K_mT_Mf=18.2_Theta=0.0_Phi=90.0.csv', header=None)
#Data2 = pd.read_csv ("/home/dwi/ownCloud/Anlabshared/Dwi/researchdata/170502Diamond111onYIGFilm/ESRYIG3+2um/ESRYIG3+2um.csv", header=None)
#ttl = input (r"Graph title: " )
#pl.title (ttl)
Data.columns = ['Microwave frequency', 'Intensity']
#Data.columns = ['Intensity']
#Data2.columns = ['Intensity']
x = np.array(Data['Microwave frequency'])
x = x[116:136]
X = np.linspace(2.58e9,2.66e9,300)
y = np.array(Data['Intensity']/Data['Intensity'].max())
offset = 1.0 - np.mean (y[0:30])
y = y[116:136]
#y2 = np.array(Data2['Intensity']/Data2['Intensity'].max())
yo = y + offset
#offset2 = 1.0 - np.mean (y2[0:30])
#yo2 = y2 + offset2

x0 = np.array([8.0e6, 2.62e9])

def func (x, W, f):
    num = W**2
    denum = pi *(x - f)**2 + pi*(W)**2 
    return 1 - (num/denum)
    
    
popt, pcov = curve_fit (func, x, yo, p0=x0)



pl.rc('text', usetex=True)
pl.rc('font', family='serif', size = 20)
#pl.plot (X, func(X, *popt), lw = 3.0, color = 'r')
pl.plot (x, yo, 'bo-', color = 'b')
#pl.plot (x, yo2, lw = 2.5, label = 'pos 1 + 2 $\mu$m', color = 'b')
#pl.ylim(0.95, 1.005)
pl.xlabel(r"$f$ (GHz)")
pl.ylabel('Intensity (norm.)')

#pl.legend(loc=0, ncol=1, fontsize = 12)
pl.savefig(path + 'ESRfit_H183Oe_DT-10K.png', format = 'png', bbox_inches = 'tight', frameon = False)
pl.show()
print('FWHM =' + str(round(popt[0],2)))
print('Center frequency =' + str(round(popt[1],2)))
print('contrast = '+ str((1.0 - np.min(yo))*100))
print (1 - np.min(yo))
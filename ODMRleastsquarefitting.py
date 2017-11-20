# -*- coding: utf-8 -*-
"""
Created on Thu Jan 19 12:52:22 2017

@author: dwi
"""


import numpy as np
from numpy import pi
import matplotlib.pyplot as pl
from scipy.optimize import leastsq
import pandas as pd

Data = pd.read_csv("/home/dwi/ownCloud/Anlabshared/Dwi/researchdata/170118_temperature_callibration/2017_01_18_350K.csv")
Data.columns = ['frequency', 'intensity']
x = np.array(Data['frequency'])
y = np.array(Data['intensity'] / Data['intensity'].max())
offset = 1.0 - np.mean (y[0:20])
yo = y + offset
c = 0.003
H = 1.5e6
D = 2.87e9
E = 6.0e6
f1 = 2.862e9
f2 = 2.870e9
p = [c, f1, f2]

def lorentzian (x, p):
     num1 = H 
     denum1 = (x - (p[1]))**2.0 + H ** 2.0
     num2 = H 
     denum2 = (x - (p[2]))**2.0 + H ** 2.0
     lor1 = p[0] * (num1 / denum1) / pi
     lor2 = p[0] * (num2 / denum2) / pi
     return (1.0 - (lor1 + lor2))
     
def residuals(p, yo, x):
    return yo - lorentzian(x, p)

pbest = leastsq(residuals, p, args=(yo, x), full_output=1)
best_parameters = pbest[0]

Lor1 = 1.0 - ( H / ((x - best_parameters[1]) ** 2 + H ** 2) * best_parameters[0] / pi )
Lor2 = 1.0 - ( H / ((x - best_parameters[2]) ** 2 + H ** 2) * best_parameters[0] / pi )
D = (best_parameters[2] + best_parameters[1]) / 2.0

fit = lorentzian(x, best_parameters)
print (best_parameters)
print (D)
#print (fit)
pl.plot(x, fit)
pl.plot(x, yo)
pl.plot (x, Lor1, x, Lor2)
pl.show()

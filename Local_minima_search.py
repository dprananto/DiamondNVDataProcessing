import pandas as pd
import numpy as np
import matplotlib.pyplot as pl
from scipy.optimize import leastsq
from scipy.signal import argrelmin 
import sys
#from numpy import NaN, Inf, arange, isscalar, asarray, array

Data = pd.read_csv ('/home/dwi/ownCloud/Anlabshared/Dwi/Data and Analysis/160627/2016_06_27_MagPar_pos13.csv')
Data.columns = ['Microwave frequency', 'Intensity']
x = np.array(Data['Microwave frequency'])
y = np.array(Data['Intensity']/Data['Intensity'].max())


dipind = argrelmin(y, order=10, mode='clip')
ind = np.array(dipind[0])

xdip = x[dipind]
ydip = y[dipind]

p = [3550000, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03]  # [B, HWHM, dip]
z = [1*10**(-4), 10, 5] #B, Theta, Phi

def zeeman(z):
    f = 2.87*10**9 + (28*10**9 * z[0] * np.sqrt((np.sin (z[1]) * np.cos(z[2]))**2 + (np.sin(z[1]) * np.sin(z[2]))**2 + (np.cos(z[2]))**2))
    return f

def Lorentzian(x, p):
    num = p[0]**2
    denum = ( x - x[ind[1]])**2 + p[0]**2
    denum2 = ( x - x[ind[2]])**2 + p[0]**2
    denum3 = ( x - x[ind[3]])**2 + p[0]**2
    denum4 = ( x - x[ind[4]])**2 + p[0]**2
    denum5 = ( x - x[ind[5]])**2 + p[0]**2
    denum6 = ( x - x[ind[6]])**2 + p[0]**2
    denum7 = ( x - x[ind[7]])**2 + p[0]**2
    denum8 = ( x - x[ind[8]])**2 + p[0]**2
    Lor = p[1]*(num/denum)
    Lor2 = p[2]*(num/denum2)
    Lor3 = p[3]*(num/denum3)
    Lor4 = p[4]*(num/denum4)
    Lor5 = p[5]*(num/denum5)
    Lor6 = p[6]*(num/denum6)
    Lor7 = p[7]*(num/denum7)
    Lor8 = p[8]*(num/denum8)
    return 1-(Lor + Lor2 + Lor3 + Lor4 + Lor5 + Lor6 + Lor7 + Lor8)

def residuals(p,y,x):
    err = y - Lorentzian(x,p)
    return err

def reszeeman(z, xdip):
    err = xdip - zeeman(z)
    return err

pbest = leastsq(residuals, p ,args=(y,x),full_output=1)
best_parameters1 = pbest[0]

pbest2 = leastsq(reszeeman, z, args=(xdip), full_output=1)
best_parameters2 = pbest2[0] 

fit = Lorentzian(x,best_parameters1)
fit2 = zeeman(z)

pl.rc('text', usetex=True)
pl.rc('font', family='serif')
pl.plot (x, y)
pl.plot (x[dipind],y[dipind], 'ro')
pl.plot(x, fit, lw=2, color = 'r') 
#pl.plot(fit2, ydip)
pl.show()
#print dipind
#print ind[1]
#print x[ind[1]]
#print p[1]
#print Lor
print best_parameters2
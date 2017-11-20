import pandas as pd
import matplotlib.pyplot as pl
import numpy as np
from scipy.optimize import leastsq
Data = pd.read_csv('/home/dwi/ownCloud/Anlabshared/Dwi/Data and Analysis/160525/2016_05_25_magparticle_.post1.csv')
Data.columns = ['Microwave Frequency [GHz]', 'Intensity']
#x = Data['Microwave Frequency [GHz]']
#y = Data['Intensity']

v = [0.5*10**-4, 90, 30] #B, tetha, phi
#fpeak = 2.87*10**9
gamma= 3550000 #2.869*10**9

f=np.linspace(2801500000, 2950000000, 400)

def f1(v):
    return 2.87*10**9 - (28*10**9 * v[0] * np.sqrt((np.sin (v[1]) * np.cos(v[2]))**2 + (np.sin(v[1]) * np.sin(v[2]))**2 + (np.cos(v[2]))**2))

def f2(v):
    return 2.87*10**9 + (28*10**9 * v[0] * np.sqrt((np.sin (v[1]) * np.cos(v[2]))**2 + (np.sin(v[1]) * np.sin(v[2]))**2 + (np.cos(v[2]))**2))

def L1(f, f1, gamma, v):
    """ Return Lorentzian line shape at x with HWHM gamma """
    return 10**-5.8 - gamma / np.pi / ((f-f1(v))**2 + gamma**2)
    
def L2(f, f1, gamma, v):
    """ Return Lorentzian line shape at x with HWHM gamma """
    return 10**-5.8 - gamma / np.pi / ((f-f2(v))**2 + gamma**2) 

# optimization # 
#fbest = leastsq(residuals,f,full_output=1)
#best_parameters = fbest[0]

# fit to data #
#fit = L1(f, f1,gamma)

#fit = np.polyfit(Data['Microwave Frequency [GHz]'], Data['Intensity'], 30)

p1=L1(f, f1, gamma, v)
p2=L2(f, f2, gamma, v)
pl.plot (Data['Microwave Frequency [GHz]'], Data['Intensity']/Data['Intensity'].max())
pl.plot(f, p1/p1.max(), linewidth = 1.5, color="Red")
pl.plot(f, p2/p2.max(), linewidth = 1.5, color="Red")
pl.title ('ESR spectrum of diamond ensemble with magnetic particle')
pl.xlabel('Microwave Frequency [Ghz]')
pl.ylabel ('Intensity [cps.]')
pl.show()
#print Data
print f1(v)
print f2(v)
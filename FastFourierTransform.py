# -*- coding: utf-8 -*-
"""
Created on Thu Aug  3 18:36:51 2017

@author: dwi
"""

import numpy as np 
from numpy import fft
import matplotlib.pyplot as pl
import pandas as pd

path = r'/home/dwi/ownCloud/Anlabshared/Dwi/researchdata/170826SpinSeebeck110Diamond/Ramsey/'
filename = r'Ramsey_0to1us_2602MHz_15dBm_DT6K_180G'
Data = pd.read_csv(path + filename + r'.csv', header=None)
Data.columns = ['Intensity']
I = np.array (Data['Intensity']/Data['Intensity'].max())
t = np.linspace(0, 1e-6, 128)


n = len(I)
s = fft.fft(I)/n
timestep = t[1] - t[0]
freq = np.fft.fftfreq(n, d=timestep)
S = np.abs(s)

pl.rc('text', usetex=True)
pl.rc('font', family='serif', size = 22)
pl.xlabel('f (Hz)')
pl.ylabel('Intensity(arb.)')
pl.xlim(0, 6e7)
pl.plot(freq, s.imag, 'bo-')
pl.savefig(path + filename + r'.png', format = 'png', bbox_inches = 'tight', frameon = False)
pl.show()
print(n)
#print(np.abs(s))

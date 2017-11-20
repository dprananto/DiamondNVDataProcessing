# -*- coding: utf-8 -*-
"""
Created on Mon Sep  4 11:37:26 2017

@author: dwi
"""

import numpy as np 
from numpy import fft
import matplotlib.pyplot as pl
import pandas as pd

path = r'/home/dwi/ownCloud/Anlabshared/Dwi/researchdata/170923SpinSeebeck110Diamond/Rabi/'
filename = r'Rabi_0to2us_2611MHz_10dBm_DT6K_Tdia310K_190G_distance80um_hyperfine'
Data = pd.read_csv(path + filename + r'.csv', header=None)
Data.columns = ['Intensity']
I = np.array (Data['Intensity']/Data['Intensity'].max())
t = np.linspace(0, 2e-6, 51)

s = fft.fft(I)
n = len(I)
timestep = t[1]-t[0]
freq = np.fft.fftfreq(n, d=timestep)

pl.rc('text', usetex=True)
pl.rc('font', family='serif', size = 22)
pl.xlabel('f (Hz)')
pl.ylabel('Intensity(arb.)')
pl.xlim(0, 0.6e7)
pl.plot(freq, s.imag, 'bo-')
pl.savefig(path + r'FFT'+ filename + r'.png', format = 'png', bbox_inches = 'tight', frameon = False)
pl.show()
#print(t)
#print(s.imag.max())
#print(freq)
#print(s.imag)

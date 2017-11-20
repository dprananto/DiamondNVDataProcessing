# -*- coding: utf-8 -*-
"""
Created on Wed Jul 12 11:22:24 2017

@author: dwi
"""

import numpy as np
import pandas as pd
import glob
import matplotlib.pyplot as pl
import re

numbers = re.compile(r'(\d+)')
def numericalSort(value):
    parts = numbers.split(value)
    parts[1::2] = map(int, parts[1::2])
    return parts

path = '/home/dwi/ownCloud/Anlabshared/Dwi/researchdata/171101SSENVC_NVA/171101_SSE_NVC/'
files = sorted(glob.glob(path + r"/number0*00.csv"), key=numericalSort)

#Reference spectra
rfile = r'number0_I=0.0000'
rData = pd.read_csv(path + rfile + r'.csv', header=None)
rData.columns = ['freq', 'S11', 'phase', 'SWR', 'Imaginary' ]
ri = np.array(rData['S11']) #reference data

ilist = [] 
for f in files:
    Data = pd.read_csv(f, header=None)
    Data.columns = ['freq', 'S11', 'phase', 'SWR', 'Imaginary' ]
    i = np.array(Data['S11'])
    si = i - ri #Subtracting data with reference data
    ilist.append(si)

I = np.vstack(ilist)
F = np.linspace(1e9, 6e9, 5001)#
y = np.array([F])
Vdc = np.linspace(0, 1.01, 102) #Voltage of Advantest Source Meter
V = 5.311823977 * Vdc - 0.0039008842 #Voltage of Power Source Instrument
H = 27.816023953826136 + 125.41195086834962 * V #Conversion to External Fields
x = np.array([H])

Y, X = np.meshgrid(y, x)

pl.rc('text', usetex=True)
pl.rc('font', family='serif', size = 22)
pl.xlabel('$H$ (Oe)', fontsize=22)
pl.ylabel(r'$f$ (GHz)', fontsize = 22)
pl.xlim(50, 233)
pl.ylim(2.4e9,3.3e9)
pl.pcolormesh(X, Y, I)
pl.colorbar(label = '$S_{11}$ (dB)')
#pl.title(r'Temperature Gradient Dependence of NV ESR on YIG with No External Field')
pl.savefig(path + r"FMRSpectogram.png", format = 'png', bbox_inches = 'tight', frameon = False)
pl.show()

print(ri)


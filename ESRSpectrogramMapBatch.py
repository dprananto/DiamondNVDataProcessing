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

path = '/home/dwi/ownCloud/Anlabshared/Dwi/researchdata/170826SpinSeebeck110Diamond/ESR/ESRSweepH5to25/'
files = sorted(glob.glob(r"/home/dwi/ownCloud/Anlabshared/Dwi/researchdata/170826SpinSeebeck110Diamond/ESR/ESRSweepH5to25/*90.0.csv"), key=numericalSort)


ilist = [] 
for f in files:
    Data = pd.read_csv(f)
    Data.columns = ['Frequency', 'Intensity']
    i = np.array(Data['Intensity']/Data['Intensity'].max())
    offset = 1.0 - np.mean (i[0:20])
    io = i + offset
    ilist.append(io)

I = np.vstack(ilist)
f = np.linspace(2.4e9, 3.3e9, 501)#
y = np.array([f])
H = np.linspace(50, 233, 64)
x = np.array([H])

Y, X = np.meshgrid(y, x)

pl.rc('text', usetex=True)
pl.rc('font', family='serif', size = 22)
pl.xlabel('H (Oe)', fontsize=22)
pl.ylabel(r'f (GHz)', fontsize = 22)
pl.xlim(50, 233)
pl.ylim(2.4e9,3.3e9)
pl.pcolormesh(X, Y, I, vmax = 1.0, vmin = 0.945)
pl.colorbar(label = 'Intensity (norm.)')
#pl.title(r'Temperature Gradient Dependence of NV ESR on YIG with No External Field')
pl.savefig(path + r"ESR_110_DT18K_batch.png", format = 'png', bbox_inches = 'tight', frameon = False)
pl.show()

#print(f)


# -*- coding: utf-8 -*-
"""
Created on Wed Aug 24 09:06:21 2016

@author: dwi
"""

###############################################################################
#THIS PROGRAM PLOT ESR SPECTRUM ###############################################
###############################################################################
import pandas as pd
import numpy as np
import matplotlib.pyplot as pl
from scipy.signal import argrelmin
#from detect_peaks import detect_peaks

fi = input ('File path to open: ')
Data = pd.read_csv (fi)
ttl = input (r"Graph title: " )
pl.title (ttl)
Data.columns = ['Microwave frequency', 'Intensity']
x = np.array(Data['Microwave frequency'])
y = np.array(Data['Intensity']/Data['Intensity'].max())
offset = 1.0 - np.mean (y[0:30])
yo = y + offset

dip = argrelmin (yo, order = 12)
#index = np.array[dip[0]]

xdip = x[dip]
ydip = yo[dip]

pl.rc('text', usetex=True)
pl.rc('font', family='serif', size = 20)
pl.plot (x, yo, lw = 2, label = 'Experiment', color = 'g')
pl.xlabel('Microwave frequency (GHz)', fontsize = 20)
pl.ylabel('Intensity (norm.)', fontsize = 20)
pl.plot (xdip, ydip, 'ro')
#pl.plot (x, Y)

#print (dip)
print (xdip)
print (ydip)
#print (yo)
#pl.legend(loc=9, ncol=3)
pl.show()

# -*- coding: utf-8 -*-
"""
Created on Thu Dec  8 11:09:53 2016

@author: dwi
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as pl


#fi = input ('File path to open: ')
Data = pd.read_csv ('/home/dwi/ownCloud/Anlabshared/Dwi/researchdata/160830 Magnetic particle stray mapping/scan 25 x 25 um/Analysis/analysis_selected_2.csv')

x = np.array([-25.0, -20.0, -15.0, -10.0, -5.0, 0.0, 5.0, 10.0, 15.0, 20.0, 25.0])
y = np.array([-25.0, -20.0, -15.0, -10.0, -5.0, 0.0, 5.0, 10.0, 15.0, 20.0, 25.0])
B = np.array(Data['B'])
X, Y = np.meshgrid(x, y)
b = B.reshape((11, 11))

pl.rc('text', usetex=True)
pl.rc('font', family='serif', size = 20)
pl.xlabel(r'x ($\mu$m)', fontsize=20)
pl.ylabel(r'y ($\mu$m)', fontsize = 20)
pl.pcolormesh(X, Y, b)
pl.colorbar()
#pl.title(r'')
pl.show()
#print (x)
#print (y)
#print (I)
print (np.shape(X))
print (np.shape(Y))
print (np.shape(b))
#print(Im)
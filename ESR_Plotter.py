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

#fi = input ("File path to open: ")
path = '/home/dwi/ownCloud/researchdata/171112SSE/171109SpinSeebeckonDiamond110/ESR/'
filename = r'ESR_160Oe_2500to3200MHz_0dBm'
Data = pd.read_csv (path + filename + r'.csv', header=None)
#Data2 = pd.read_csv ("/home/dwi/ownCloud/Anlabshared/Dwi/researchdata/170502Diamond111onYIGFilm/ESRYIG3+2um/ESRYIG3+2um.csv", header=None)
#ttl = input (r"Graph title: " )
#pl.title (ttl)
#Data.columns = ['Microwave frequency', 'Intensity']
Data.columns = ['Intensity']
#Data2.columns = ['Intensity']
#x = np.array(Data['Microwave frequency'])
x = np.linspace(2.50e9, 3.2e9,201) #Start freq, Stop freq, N points
y = np.array(Data['Intensity']/Data['Intensity'].max())
#y2 = np.array(Data2['Intensity']/Data2['Intensity'].max())
offset = 1.0 - np.mean (y[0:30])
yo = y + offset
#offset2 = 1.0 - np.mean (y2[0:30])
#yo2 = y2 + offset2

pl.rc('text', usetex=True)
pl.rc('font', family='serif', size = 20)
pl.plot (x, yo, lw = 2.0, label = 'pos 1', color = 'b')
#pl.plot (x, yo2, lw = 2.5, label = 'pos 1 + 2 $\mu$m', color = 'b')
#pl.ylim(0.955, 1.005)
pl.xlabel(r"$f$ (GHz)")
pl.ylabel('Intensity (norm.)')

#pl.legend(loc=0, ncol=1, fontsize = 12)
pl.savefig(path + filename + '.png', format = 'png', bbox_inches = 'tight', frameon = False)
#pl.show()
#print(x)
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 24 13:56:41 2017

@author: dwi
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as pl
import cmocean
from matplotlib import ticker

path = "/home/dwi/ownCloud/Anlabshared/Dwi/Projects/Spin-waveHeatConveyer/LockinThermographyData/20170222/AllCSVAndFigures/20170222_YIGdisc_H+9.0A(gap80mm)_MW4.7GHz24(-6)dBm_(12.5Hz)_cal_ampl"
Data = pd.read_csv(path+".csv", usecols = range(1, 640) , skiprows = 4)
I = np.array(Data) #converting data into array
x = np.arange(0, 640, 1)
y = np.arange(0, 512, 1)
X, Y = np.meshgrid(x, y)

pl.rc('text', usetex=True)
pl.rc('font', family='serif', size = 20)
pl.xlim(0, 640)
pl.ylim(0, 512)
pl.gca().invert_yaxis() #Inverting Y axis
#pl.imshow(I, cmap = 'jet')
pl.pcolormesh(X, Y, I, cmap = cmocean.cm.thermal, vmax = 1.0, vmin = 0.0) #Amplitude
#pl.pcolormesh(X, Y, I, cmap = 'hsv_r', vmax = 360.0, vmin = 0.0) #Phase
cb = pl.colorbar(label = 'A (mK)')
#cb = pl.colorbar(label = r'\phi (deg)')
tick_locator = ticker.MaxNLocator(nbins=6)
cb.locator = tick_locator
cb.update_ticks()
frame1 = pl.gca()
frame1.axes.get_xaxis().set_visible(False)
frame1.axes.get_yaxis().set_visible(False)

pl.savefig(path+".png", format = 'png')
pl.show()

#print (I)
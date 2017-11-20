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

path = "/home/dwi/ownCloud/Anlabshared/Dwi/researchdata/171101SSENVC_NVA/171101_SSE_NVC/diff_by38Oe/"
filename = r'H_dependence'
Data = pd.read_csv(path+ filename + r".csv", usecols = range(1, 239))
I = np.array(Data) #converting data into array
x = np.linspace(38, 842, 238)
y = np.linspace(1e9, 6e9, 100)
X, Y = np.meshgrid(x, y)

pl.rc('text', usetex=True)
pl.rc('font', family='serif', size = 20)
#pl.xlim(0, 640)
#pl.ylim(0, 512)
pl.gca().invert_yaxis() #Inverting Y axis
#pl.imshow(I, cmap = 'jet')
pl.pcolormesh(X, Y, I) #Amplitude
#pl.pcolormesh(X, Y, I, cmap = 'hsv_r', vmax = 360.0, vmin = 0.0) #Phase
cb = pl.colorbar(label = 'dBm')
#cb = pl.colorbar(label = r'\phi (deg)')
tick_locator = ticker.MaxNLocator(nbins=6)
cb.locator = tick_locator
cb.update_ticks()
frame1 = pl.gca()
#frame1.axes.get_xaxis().set_visible(False)
#frame1.axes.get_yaxis().set_visible(False)

pl.savefig(path+ filename + r".png", format = 'png')
pl.show()

#print (I)
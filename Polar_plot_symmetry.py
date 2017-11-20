# -*- coding: utf-8 -*-
"""
Created on Fri Sep 23 14:10:27 2016

@author: dwi
"""

import numpy as np
import matplotlib.pyplot as pl
import pandas as pd

path = r'/home/dwi/ownCloud/Anlabshared/Dwi/researchdata/160830 Magnetic particle stray mapping/scan 25 x 25 um/Analysis/'
fi = r'/home/dwi/ownCloud/Anlabshared/Dwi/researchdata/160830 Magnetic particle stray mapping/scan 25 x 25 um/Analysis/Pos93_for_Symmetry_plotC2.csv'
Data = pd.read_csv (fi)
#angle = input ('Select angle to plot (Vertical or Horizontal): ')
phi = np.array(Data['T'] * np.pi/180.0) 
B = np.array(Data['B'] * 1e3)

pl.rc('text', usetex=False)
pl.rc('font', family='sans', size = 18)
#pl.title(r"$\phi$")
pl.figure()
ax = pl.axes(polar = True)
#ax.set_rmax(3)
#ax.set_rticks([1.0, 1.5, 2.0, 2.5])  # less radial ticks
ax.set_rlabel_position(0)  # get radial labels away from plotted line
#ax.grid(True)
pl.yticks([0.0, 0.5, 1.0, 1.5])
pl.ylim(0, 1.7)
pl.xlabel(r"$\theta$", fontsize = 22, fontname = "arial")
#pl.ylabel('B (mT)', fontsize = 22)
pl.plot(phi, B, 'b^', markersize = 10)
pl.savefig(path + 'Polar_Theta_P93.png', format = 'png', bbox_inches = 'tight', frameon = False)
pl.show()
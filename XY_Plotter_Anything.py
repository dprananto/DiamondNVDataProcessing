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

#fi = input ('File path to open: ')
Data = pd.read_csv ('/home/dwi/ownCloud/Anlabshared/Dwi/researchdata/170118_temperature_callibration/ZFS_temperature_nofield.csv')
ttl = input (r"Graph title: " )
pl.title (ttl)
#Data.columns = ['X', 'Y']
x = np.array(Data['Temperature'])
y = np.array(Data['ZFSfit'])
#y2 = np.array(Data['frequency2'])
#y3 = np.array(Data['I3'])
#y4 = np.array(Data['I4'])
#y5 = np.array(Data['I5'])
#y6 = np.array(Data['I6'])
#y7 = np.array(Data['I7'])
#y8 = np.array(Data['I8'])

#Y = 1.0 - y
#Y2 = 1.0 - y2
#Y3 = 1.0 - y3
#Y4 = 1.0 - y4

pl.rc('text', usetex=True)
pl.rc('font', family='serif', size = 20)
pl.plot (x, y, 'ro', markersize = 10.0, label = r'$\downarrow$')
#pl.plot (x, y2, 'bo', markersize = 10.0, label= r'$\uparrow$')
#pl.plot (x, Y3, 'go', markersize = 10.0, label= 'Dip 3')
#pl.plot (x, Y4, 'yo', markersize = 10.0, label= 'Dip 4')
#pl.plot (x, y5, 'yo', markersize = 10.0, label= 'Dip 5')
#pl.plot (x, y6, 'go', markersize = 10.0, label= 'Dip 6')
#pl.plot (x, y7, 'bo', markersize = 10.0, label= 'Dip 7')
#pl.plot (x, y8, 'ro', markersize = 10.0, label= 'Dip 8')

pl.xlabel(r"Temperature (K)", fontsize = 20)
pl.ylabel(r'D (GHz)', fontsize = 20)

pl.legend(loc=1, ncol=4)
pl.show()

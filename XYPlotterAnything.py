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
path = r'/home/dwi/ownCloud/Anlabshared/Dwi/researchdata/171101SSENVC/'
filename = r'171101SSENVC_B1.3A1.39V_1to4.5GHz_10dBm'
Data = pd.read_csv (path + filename + r'.csv')
#ttl = input (r"Graph title: " )
#pl.title ('Microwave Absorption Spectrum of YIG Disk H = 670 Oe')
#Data.columns = ['f', 'absorption']
x = np.array(Data['freq[Hz]'])
y = np.array(Data['Trc1_S11[dB]'])
#y2 = np.array(Data['f2'])
#y3 = np.array(Data['f3'])
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
pl.plot (x, y, lw = 3.0, label = 'pos 1', color = 'k')
#pl.plot (x, y, 'o', markersize = 15.0, markerfacecolor='None', markeredgewidth=2 ,markeredgecolor='k', label = r'$f_1$')
#pl.plot (x, y2, 'o', markersize = 15.0, markerfacecolor='None', markeredgewidth=2 ,markeredgecolor='g', label = r'$f_2$')
#pl.plot (x, y3, 'o', markersize = 15.0, markerfacecolor='None', markeredgewidth=2 ,markeredgecolor='b', label = r'$f_3$')
#pl.plot (x, Y4, 'yo', markersize = 10.0, label= 'Dip 4')
#pl.plot (x, y5, 'yo', markersize = 10.0, label= 'Dip 5')
#pl.plot (x, y6, 'go', markersize = 10.0, label= 'Dip 6')
#pl.plot (x, y7, 'bo', markersize = 10.0, label= 'Dip 7')
#pl.plot (x, y8, 'ro', markersize = 10.0, label= 'Dip 8')

pl.ylabel(r"Absorption, $S_{11}$ (dB)", fontsize = 20)
pl.xlabel(r'$f$ (GHz) ', fontsize = 20)
#pl.ylim(2.7, 3.6)
#pl.xlim(-16, 16)
###########################PRINTING FMR FREQUENCY#############################
#minFMR = y[0:1096].min()


#for index, item in enumerate(y):
 #   if item == minFMR:
  #      print(Data['f'].iloc[index])   
###############################################################################

#pl.legend(loc=0, ncol=3, fontsize=16)
pl.savefig(path + filename + '.png', format = 'png', bbox_inches = 'tight', frameon = False)
pl.show()

# -*- coding: utf-8 -*-
"""
Created on Wed Aug 24 09:06:21 2016
Modified on Mon Nov 6 12:00:00 2017
@author: dwi
"""

###############################################################################
#THIS PROGRAM PLOT FMR SPECTRUM ###############################################
###############################################################################
import pandas as pd
import numpy as np
import matplotlib.pyplot as pl

#fi = input ('File path to open: ')
path = r'/home/dwi/ownCloud/researchdata/171112SSE/VNA/'
filename = r'171110233V0.18ASSE'
filenameref = r'171110210V0.05ASSE' #Reference Filename
Data = pd.read_csv (path + filename + r'.csv')
Dataref = pd.read_csv (path + filenameref + r'.csv') #Reference Data
#ttl = input (r"Graph title: " )
#pl.title ('Microwave Absorption Spectrum of YIG Disk H = 670 Oe')
Data.columns = ['freq', 'S11']#, 'phase', 'SWR', 'Imaginary' ]
Dataref.columns = ['freq', 'S11']#, 'phase', 'SWR', 'Imaginary' ]
x = np.array(Data['freq'])
y = np.array(Data['S11'])
yref = np.array(Dataref['S11'])
ysub = y - yref #Subtraction with reference
#y2 = np.array(Data['f2'])
#y3 = np.array(Data['f3'])
#y4 = np.array(Data['I4'])
#y5 = np.array(Data['I5'])
#y6 = np.array(Data['I6'])
#y7 = np.array(Data['I7'])
#y8 = np.array(Data['I8'])

#####DEFINE MAGNETIC FIELD#####################################################
#Vdc = float(filename[10:])
#V = 5.311823977 * Vdc - 0.0039008842 #Voltage of Power Source Instrument
#H = 27.816023953826136 + 125.41195086834962 * V #Conversion to External Fields

#Y = 1.0 - y
#Y2 = 1.0 - y2
#Y3 = 1.0 - y3
#Y4 = 1.0 - y4

pl.rc('text', usetex=True)
pl.rc('font', family='serif', size = 20)
pl.plot (x, ysub, lw = 2.5, label = 'pos 1', color = 'k')
#pl.plot (x, y, 'o', markersize = 15.0, markerfacecolor='None', markeredgewidth=2 ,markeredgecolor='k', label = r'$f_1$')
#pl.plot (x, y2, 'o', markersize = 15.0, markerfacecolor='None', markeredgewidth=2 ,markeredgecolor='g', label = r'$f_2$')
#pl.plot (x, y3, 'o', markersize = 15.0, markerfacecolor='None', markeredgewidth=2 ,markeredgecolor='b', label = r'$f_3$')
#pl.plot (x, Y4, 'yo', markersize = 10.0, label= 'Dip 4')
#pl.plot (x, y5, 'yo', markersize = 10.0, label= 'Dip 5')
#pl.plot (x, y6, 'go', markersize = 10.0, label= 'Dip 6')
#pl.plot (x, y7, 'bo', markersize = 10.0, label= 'Dip 7')
#pl.plot (x, y8, 'ro', markersize = 10.0, label= 'Dip 8')

pl.ylabel(r"$S_{11}$ (dB)", fontsize = 20)
pl.xlabel(r'$f$ (GHz) ', fontsize = 20)
#pl.ylim(2.7, 3.6)
pl.xlim(1e9, 4e9)
###########################PRINTING FMR FREQUENCY#############################
#minFMR = y[0:1096].min()


#for index, item in enumerate(y):
 #   if item == minFMR:
  #      print(Data['f'].iloc[index])   
###############################################################################

#pl.legend(loc=0, ncol=3, fontsize=16)
#pl.savefig(path + filename[0:7] + r'H' + str(round(H, 2)) + r'Oe' + '.png', format = 'png', bbox_inches = 'tight', frameon = False)
pl.savefig(path + filename[0:7] + '.png', format = 'png', bbox_inches = 'tight', frameon = False)
#pl.title(r'$H$=' + str(round(H, 2)) + r'Oe')
pl.show()
#print(H)

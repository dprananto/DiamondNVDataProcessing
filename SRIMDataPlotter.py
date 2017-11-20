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
path = r'/home/dwi/ownCloud/Anlabshared/Dwi/Projects/CoherentManipulationWithSpinWave/SRIM/'
Data = pd.read_csv (r'/home/dwi/ownCloud/Anlabshared/Dwi/Projects/CoherentManipulationWithSpinWave/SRIM/N14IonstoDiamond30keV.csv')
#ttl = input (r"Graph title: " )
#pl.title ('Microwave Absorption Spectrum of YIG Disk H = 670 Oe')
density = 1.13e23 #atomic density of Carbon (atoms/cm^3)
dose = 2.0e11 #implantation dose (ions/cm^2)
Data.columns = ['Depth', 'N', 'C']
x = np.array(Data['Depth'])/10.0
y = np.array(Data['N']) * dose
y2 = np.array(Data['C']) * dose
#y3 = np.array(Data['T1MWPulse2ms'])
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

fig, ax1 = pl.subplots()

#Make Twin y-axis
ax2 = ax1.twinx()
ax1.plot(x, y, 'b-', lw = 3.0, label = r'N')
ax2.plot(x, y2, 'g-', lw = 3.0, label= r'C')

ax1.set_xlabel(r"Depth (nm)", fontsize = 20)
ax1.set_ylabel(r'$^{14}$N (cm$^{-3}$)', fontsize = 20, color = 'b')
#ax1.set_yscale('log')
ax2.set_ylabel(r'V (cm$^{-3}$)', fontsize = 20, color = 'g')
#ax2.set_yscale('log')

#pl.plot (x, y, 'bo', markersize = 10.0, label = r'N')
#pl.plot (x, y2, 'go', markersize = 10.0, label= r'C')
#pl.plot (x, y3, 'ro', markersize = 10.0, label= 'MW Pulse 2 ms')
#pl.plot (x, Y4, 'yo', markersize = 10.0, label= 'Dip 4')
#pl.plot (x, y5, 'yo', markersize = 10.0, label= 'Dip 5')
#pl.plot (x, y6, 'go', markersize = 10.0, label= 'Dip 6')
#pl.plot (x, y7, 'bo', markersize = 10.0, label= 'Dip 7')
#pl.plot (x, y8, 'ro', markersize = 10.0, label= 'Dip 8')

#pl.xlabel(r"Position z ($\AA$)", fontsize = 20)
#pl.ylabel(r'(Atoms/cm$^3$)/(Atoms/cm$^2$)', fontsize = 20)

#pl.ylim(2, 6)
pl.xlim(0, 80)
###########################PRINTING FMR FREQUENCY#############################
#minFMR = y[0:1096].min()


#for index, item in enumerate(y):
 #   if item == minFMR:
  #      print(Data['f'].iloc[index])   
###############################################################################

#ax1.legend(loc=2, ncol=2, fontsize=14)
#ax2.legend(loc=1, ncol=2, fontsize=14)
pl.savefig(path + 'SRIM_N14_30keV.eps', format = 'eps', bbox_inches = 'tight')
pl.show()

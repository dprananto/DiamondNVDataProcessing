# -*- coding: utf-8 -*-
"""
Created on Wed Aug 24 09:06:21 2016

@author: dwi
"""

###############################################################################
#THIS PROGRAM PLOT MULTIPLE ESR SPECTRUM ###############################################
###############################################################################
import pandas as pd
import numpy as np
import matplotlib.pyplot as pl

fi1 = input ('File1 path to open: ')
Data1 = pd.read_csv (fi1)
Label1 = input ('Data Label1:')

fi2 = input ('File2 path to open: ')
Data2 = pd.read_csv (fi2)
Label2 = input ('Data Label2:')

#fi3 = input ('File3 path to open: ')
#Data3 = pd.read_csv (fi3)
#Label3 = input ('Data Label3:')

#fi4 = input ('File4 path to open: ')
#Data4 = pd.read_csv (fi4)
#Label4 = input ('Data Label4:')

#fi5 = input ('File5 path to open: ')
#Data5 = pd.read_csv (fi5)
#Label5 = input ('Data Label5:')

#fi6 = input ('File6 path to open: ')
#Data6 = pd.read_csv (fi6)
#Label6 = input ('Data Label6:')

#fi7 = input ('File7 path to open: ')
#Data7 = pd.read_csv (fi7)
#Label7 = input ('Data Label7:')

#fi8 = input ('File8 path to open: ')
#Data8 = pd.read_csv (fi8)
#Label8 = input ('Data Label8:')

#fi9 = input ('File9 path to open: ')
#Data9 = pd.read_csv (fi9)
#Label9 = input ('Data Label9:')

#fi10 = input ('File10 path to open: ')
#Data10 = pd.read_csv (fi10)
#Label10 = input ('Data Label10:')

#fi11 = input ('File10 path to open: ')
#Data11 = pd.read_csv (fi11)
#Label11 = input ('Data Label11:')

ttl = input (r"Graph title: " )
pl.title (ttl)
Data1.columns = ['Microwave frequency', 'Intensity']
Data2.columns = ['Microwave frequency', 'Intensity']
#Data3.columns = ['Microwave frequency', 'Intensity']
#Data4.columns = ['Microwave frequency', 'Intensity']
#Data5.columns = ['Microwave frequency', 'Intensity']
#Data6.columns = ['Microwave frequency', 'Intensity']
#Data7.columns = ['Microwave frequency', 'Intensity']
#Data8.columns = ['Microwave frequency', 'Intensity']
#Data9.columns = ['Microwave frequency', 'Intensity']
#Data10.columns = ['Microwave frequency', 'Intensity']
#Data11.columns = ['Microwave frequency', 'Intensity']

x1 = np.array(Data1['Microwave frequency'])
y1 = np.array(Data1['Intensity']/Data1['Intensity'].max())

x2 = np.array(Data2['Microwave frequency'])
y2 = np.array(Data2['Intensity']/Data2['Intensity'].max()) #- 0.03

#x3 = np.array(Data3['Microwave frequency'])
#y3 = np.array(Data3['Intensity']/Data3['Intensity'].max()) #- 0.06

#x4 = np.array(Data4['Microwave frequency'])
#y4 = np.array(Data4['Intensity']/Data4['Intensity'].max()) #- 0.09

#x5 = np.array(Data5['Microwave frequency'])
#y5 = np.array(Data5['Intensity']/Data5['Intensity'].max()) - 0.12

#x6 = np.array(Data6['Microwave frequency'])
#y6 = np.array(Data6['Intensity']/Data6['Intensity'].max()) - 0.15
 
#x7 = np.array(Data7['Microwave frequency'])
#y7 = np.array(Data7['Intensity']/Data7['Intensity'].max()) - 0.18

#x8 = np.array(Data8['Microwave frequency'])
#y8 = np.array(Data8['Intensity']/Data8['Intensity'].max()) - 0.21

#x9 = np.array(Data9['Microwave frequency'])
#y9 = np.array(Data9['Intensity']/Data9['Intensity'].max()) - 0.24

#x10 = np.array(Data10['Microwave frequency'])
#y10 = np.array(Data10['Intensity']/Data10['Intensity'].max()) - 0.27

#x11 = np.array(Data11['Microwave frequency'])
#y11 = np.array(Data11['Intensity']/Data11['Intensity'].max()) - 0.30

#offset = 1.0 - np.mean (y[0:30])
#yo = y + offset


pl.rc('text', usetex=True)
pl.rc('font', family='serif', size = 14)
pl.plot (x1, y1, lw = 2, label = Label1)
pl.plot (x2, y2, lw = 2, label = Label2)
#pl.plot (x3, y3, lw = 2, label = Label3)
#pl.plot (x4, y4, lw = 2, label = Label4)
#pl.plot (x5, y5, lw = 2, label = Label5)
#pl.plot (x6, y6, lw = 2, label = Label6)
#pl.plot (x7, y7, lw = 2, label = Label7)
#pl.plot (x8, y8, lw = 2, label = Label8)
#pl.plot (x9, y9, lw = 2, label = Label9)
#pl.plot (x10, y10, lw = 2, label = Label10)
#pl.plot (x11, y11, lw = 2, label = Label11)

pl.xlabel('Microwave frequency (GHz)', fontsize = 22)
pl.ylabel('Intensity (arb.)', fontsize = 22)

pl.legend(loc=8, ncol=2)
pl.show()

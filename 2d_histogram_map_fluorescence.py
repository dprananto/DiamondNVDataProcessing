# -*- coding: utf-8 -*-
"""
Created on Mon Dec  5 16:39:15 2016
This program will plot spatial intensity. Can be used as spectogram for ESR spectra
@author: dwi
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as pl

path = '/home/dwi/ownCloud/Anlabshared/Dwi/researchdata/161229 LSEE (111) Magnet Field/Analysis/'
#fi = input ('File path to open: ')
Data1 = pd.read_csv (r'/home/dwi/ownCloud/Anlabshared/Dwi/researchdata/161229 LSEE (111) Magnet Field/2016_12_29_LSSE_(111)_FieldPerpdToTgrad_0K_Mf=0.0_Theta=0.0_Phi=0.0.csv', header = None)
Data1.columns = ['Microwave frequency', 'Intensity']
Data2 = pd.read_csv (r'/home/dwi/ownCloud/Anlabshared/Dwi/researchdata/161229 LSEE (111) Magnet Field/2016_12_29_LSSE_(111)_FieldPerpdToTgrad_1K_Mf=0.0_Theta=0.0_Phi=0.0.csv', header = None)
Data2.columns = ['Microwave frequency', 'Intensity']
Data3 = pd.read_csv (r'/home/dwi/ownCloud/Anlabshared/Dwi/researchdata/161229 LSEE (111) Magnet Field/2016_12_29_LSSE_(111)_FieldPerpdToTgrad_2K_Mf=0.0_Theta=0.0_Phi=0.0.csv', header = None)
Data3.columns = ['Microwave frequency', 'Intensity']
Data4 = pd.read_csv (r'/home/dwi/ownCloud/Anlabshared/Dwi/researchdata/161229 LSEE (111) Magnet Field/2016_12_29_LSSE_(111)_FieldPerpdToTgrad_3K_Mf=0.0_Theta=0.0_Phi=0.0.csv', header = None)
Data4.columns = ['Microwave frequency', 'Intensity']
Data5 = pd.read_csv (r'/home/dwi/ownCloud/Anlabshared/Dwi/researchdata/161229 LSEE (111) Magnet Field/2016_12_29_LSSE_(111)_FieldPerpdToTgrad_4K_Mf=0.0_Theta=0.0_Phi=0.0.csv', header = None)
Data5.columns = ['Microwave frequency', 'Intensity']
Data6 = pd.read_csv (r'/home/dwi/ownCloud/Anlabshared/Dwi/researchdata/161229 LSEE (111) Magnet Field/2016_12_29_LSSE_(111)_FieldPerpdToTgrad_5K_Mf=0.0_Theta=0.0_Phi=0.0.csv', header = None)
Data6.columns = ['Microwave frequency', 'Intensity']
Data7 = pd.read_csv (r'/home/dwi/ownCloud/Anlabshared/Dwi/researchdata/161229 LSEE (111) Magnet Field/2016_12_29_LSSE_(111)_FieldPerpdToTgrad_6K_Mf=0.0_Theta=0.0_Phi=0.0.csv', header = None)
Data7.columns = ['Microwave frequency', 'Intensity']
Data8 = pd.read_csv (r'/home/dwi/ownCloud/Anlabshared/Dwi/researchdata/161229 LSEE (111) Magnet Field/2016_12_29_LSSE_(111)_FieldPerpdToTgrad_7K_Mf=0.0_Theta=0.0_Phi=0.0.csv', header = None)
Data8.columns = ['Microwave frequency', 'Intensity']
Data9 = pd.read_csv (r'/home/dwi/ownCloud/Anlabshared/Dwi/researchdata/161229 LSEE (111) Magnet Field/2016_12_29_LSSE_(111)_FieldPerpdToTgrad_8K_Mf=0.0_Theta=0.0_Phi=0.0.csv', header = None)
Data9.columns = ['Microwave frequency', 'Intensity']
Data10 = pd.read_csv (r'/home/dwi/ownCloud/Anlabshared/Dwi/researchdata/161229 LSEE (111) Magnet Field/2016_12_29_LSSE_(111)_FieldPerpdToTgrad_9K_Mf=0.0_Theta=0.0_Phi=0.0.csv', header = None)
Data10.columns = ['Microwave frequency', 'Intensity']
Data11 = pd.read_csv (r'/home/dwi/ownCloud/Anlabshared/Dwi/researchdata/161229 LSEE (111) Magnet Field/2016_12_29_LSSE_(111)_FieldPerpdToTgrad_10K_Mf=0.0_Theta=0.0_Phi=0.0.csv', header = None)
Data11.columns = ['Microwave frequency', 'Intensity']
Data12 = pd.read_csv (r'/home/dwi/ownCloud/Anlabshared/Dwi/researchdata/161229 LSEE (111) Magnet Field/2016_12_29_LSSE_(111)_FieldPerpdToTgrad_11K_Mf=0.0_Theta=0.0_Phi=0.0.csv', header = None)
Data12.columns = ['Microwave frequency', 'Intensity']
Data13 = pd.read_csv (r'/home/dwi/ownCloud/Anlabshared/Dwi/researchdata/161229 LSEE (111) Magnet Field/2016_12_29_LSSE_(111)_FieldPerpdToTgrad_12K_Mf=0.0_Theta=0.0_Phi=0.0.csv', header = None)
Data13.columns = ['Microwave frequency', 'Intensity']
Data14 = pd.read_csv (r'/home/dwi/ownCloud/Anlabshared/Dwi/researchdata/161229 LSEE (111) Magnet Field/2016_12_29_LSSE_(111)_FieldPerpdToTgrad_13K_Mf=0.0_Theta=0.0_Phi=0.0.csv', header = None)
Data14.columns = ['Microwave frequency', 'Intensity']
Data15 = pd.read_csv (r'/home/dwi/ownCloud/Anlabshared/Dwi/researchdata/161229 LSEE (111) Magnet Field/2016_12_29_LSSE_(111)_FieldPerpdToTgrad_14K_Mf=0.0_Theta=0.0_Phi=0.0.csv', header = None)
Data15.columns = ['Microwave frequency', 'Intensity']
Data16 = pd.read_csv (r'/home/dwi/ownCloud/Anlabshared/Dwi/researchdata/161229 LSEE (111) Magnet Field/2016_12_29_LSSE_(111)_FieldPerpdToTgrad_15K_Mf=0.0_Theta=0.0_Phi=0.0.csv', header = None)
Data16.columns = ['Microwave frequency', 'Intensity']
Data17 = pd.read_csv (r'/home/dwi/ownCloud/Anlabshared/Dwi/researchdata/161229 LSEE (111) Magnet Field/2016_12_29_LSSE_(111)_FieldPerpdToTgrad_16K_Mf=0.0_Theta=0.0_Phi=0.0.csv', header = None)
Data17.columns = ['Microwave frequency', 'Intensity']
Data18 = pd.read_csv (r'/home/dwi/ownCloud/Anlabshared/Dwi/researchdata/161229 LSEE (111) Magnet Field/2016_12_29_LSSE_(111)_FieldPerpdToTgrad_17K_Mf=0.0_Theta=0.0_Phi=0.0.csv', header = None)
Data18.columns = ['Microwave frequency', 'Intensity']
Data19 = pd.read_csv (r'/home/dwi/ownCloud/Anlabshared/Dwi/researchdata/161229 LSEE (111) Magnet Field/2016_12_29_LSSE_(111)_FieldPerpdToTgrad_18K_Mf=0.0_Theta=0.0_Phi=0.0.csv', header = None)
Data19.columns = ['Microwave frequency', 'Intensity']
Data20 = pd.read_csv (r'/home/dwi/ownCloud/Anlabshared/Dwi/researchdata/161229 LSEE (111) Magnet Field/2016_12_29_LSSE_(111)_FieldPerpdToTgrad_19K_Mf=0.0_Theta=0.0_Phi=0.0.csv', header = None)
Data20.columns = ['Microwave frequency', 'Intensity']
Data21 = pd.read_csv (r'/home/dwi/ownCloud/Anlabshared/Dwi/researchdata/161229 LSEE (111) Magnet Field/2016_12_29_LSSE_(111)_FieldPerpdToTgrad_20K_Mf=0.0_Theta=0.0_Phi=0.0.csv', header = None)
Data21.columns = ['Microwave frequency', 'Intensity']
Data22 = pd.read_csv (r'/home/dwi/ownCloud/Anlabshared/Dwi/researchdata/161229 LSEE (111) Magnet Field/2016_12_29_LSSE_(111)_FieldPerpdToTgrad_21K_Mf=0.0_Theta=0.0_Phi=0.0.csv', header = None)
Data22.columns = ['Microwave frequency', 'Intensity']
Data23 = pd.read_csv (r'/home/dwi/ownCloud/Anlabshared/Dwi/researchdata/161229 LSEE (111) Magnet Field/2016_12_29_LSSE_(111)_FieldPerpdToTgrad_22K_Mf=0.0_Theta=0.0_Phi=0.0.csv', header = None)
Data23.columns = ['Microwave frequency', 'Intensity']
#Data23 = pd.read_csv ('/home/dwi/ownCloud/Anlabshared/Dwi/Data and Analysis/170103 LSSE (110) Diamond/2017_01_03_LSSE_(110)Diamond_FieldPerpdTgrad_-3K_Mf=0.0_Theta=0.0_Phi=90.0.csv')
#Data23.columns = ['Microwave frequency', 'Intensity']
#Data24 = pd.read_csv ('/home/dwi/ownCloud/Anlabshared/Dwi/Data and Analysis/170103 LSSE (110) Diamond/2017_01_03_LSSE_(110)Diamond_FieldPerpdTgrad_-4K_Mf=0.0_Theta=0.0_Phi=90.0.csv')
#Data24.columns = ['Microwave frequency', 'Intensity']
#Data25 = pd.read_csv ('/home/dwi/ownCloud/Anlabshared/Dwi/Data and Analysis/170103 LSSE (110) Diamond/2017_01_03_LSSE_(110)Diamond_FieldPerpdTgrad_-5K_Mf=0.0_Theta=0.0_Phi=90.0.csv')
#Data25.columns = ['Microwave frequency', 'Intensity']

x = np.array(Data2['Microwave frequency'])
y = np.arange(0.0, 23.0, 1.0)
y = np.array(y)
i1 = np.array(Data1['Intensity']/Data1['Intensity'].max())
i2 = np.array(Data2['Intensity']/Data2['Intensity'].max())
i3 = np.array(Data3['Intensity']/Data3['Intensity'].max())
i4 = np.array(Data4['Intensity']/Data4['Intensity'].max())
i5 = np.array(Data5['Intensity']/Data5['Intensity'].max())
i6 = np.array(Data6['Intensity']/Data6['Intensity'].max())
i7 = np.array(Data7['Intensity']/Data7['Intensity'].max())
i8 = np.array(Data8['Intensity']/Data8['Intensity'].max())
i9 = np.array(Data9['Intensity']/Data9['Intensity'].max())
i10 = np.array(Data10['Intensity']/Data10['Intensity'].max())
i11 = np.array(Data11['Intensity']/Data11['Intensity'].max())
i12 = np.array(Data12['Intensity']/Data12['Intensity'].max())
i13 = np.array(Data13['Intensity']/Data13['Intensity'].max())
i14 = np.array(Data14['Intensity']/Data14['Intensity'].max())
i15 = np.array(Data15['Intensity']/Data15['Intensity'].max())
i16 = np.array(Data16['Intensity']/Data16['Intensity'].max())
i17 = np.array(Data17['Intensity']/Data17['Intensity'].max())
i18 = np.array(Data18['Intensity']/Data18['Intensity'].max())
i19 = np.array(Data19['Intensity']/Data19['Intensity'].max())
i20 = np.array(Data20['Intensity']/Data20['Intensity'].max())
i21 = np.array(Data21['Intensity']/Data21['Intensity'].max())
i22 = np.array(Data22['Intensity']/Data22['Intensity'].max())
i23 = np.array(Data23['Intensity']/Data23['Intensity'].max())
#i24 = np.array(Data24['Intensity']/Data24['Intensity'].max())
#i25 = np.array(Data25['Intensity']/Data25['Intensity'].max())

I = np.vstack((i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11, i12, i14, i15, i16, i17, i18, i19, i20, i21, i22, i23))#, i24, i25))
#I = I.transpose()
Im = np.array(I)
X, Y = np.meshgrid(x, y)

pl.rc('text', usetex=True)
pl.rc('font', family='serif', size = 20)
pl.xlabel('$f$ (GHz)', fontsize=20)
pl.ylabel(r'$\Delta$T (K)', fontsize = 20)
pl.ylim(0,22)
pl.pcolormesh(X, Y, Im)
pl.colorbar()
#pl.title(r'Temperature Gradient Dependence of NV ESR on YIG with No External Field')
pl.savefig(path + 'SpectogramFieldPerpTgrad.ps', format = 'ps', bbox_inches = 'tight', frameon = False)
pl.show()
#print (y)
#print (I)
#print (np.shape(X))
#print (np.shape(Y))
#print (np.shape(Im))
#print(Im)
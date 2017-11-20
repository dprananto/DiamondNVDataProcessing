# -*- coding: utf-8 -*-
"""
Created on Mon Dec  5 16:39:15 2016
This program will plot spatial intensity. Can be used as spectogram for ESR spectra
@author: dwi
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as pl


path = '/home/dwi/ownCloud/Anlabshared/Dwi/researchdata/170711SpinSeebeck110Diamond/ESR/'
Data1 = pd.read_csv ('/home/dwi/ownCloud/Anlabshared/Dwi/researchdata/170711SpinSeebeck110Diamond/ESR/1to3.5GHzH0OeDT5K.csv')
Data1.columns = ['Intensity']
Data2 = pd.read_csv ('/home/dwi/ownCloud/Anlabshared/Dwi/researchdata/170711SpinSeebeck110Diamond/ESR/1to3.5GHzH15OeDT5K.csv')
Data2.columns = ['Intensity']
Data3 = pd.read_csv ('/home/dwi/ownCloud/Anlabshared/Dwi/researchdata/170711SpinSeebeck110Diamond/ESR/1to3.5GHzH30OeDT5K.csv')
Data3.columns = ['Intensity']
Data4 = pd.read_csv ('/home/dwi/ownCloud/Anlabshared/Dwi/researchdata/170711SpinSeebeck110Diamond/ESR/1to3.5GHzH45OeDT5K.csv')
Data4.columns = ['Intensity']
Data5 = pd.read_csv ('/home/dwi/ownCloud/Anlabshared/Dwi/researchdata/170711SpinSeebeck110Diamond/ESR/1to3.5GHzH60OeDT5K.csv')
Data5.columns = ['Intensity']
Data6 = pd.read_csv ('/home/dwi/ownCloud/Anlabshared/Dwi/researchdata/170711SpinSeebeck110Diamond/ESR/1to3.5GHzH75OeDT5K.csv')
Data6.columns = ['Intensity']
Data7 = pd.read_csv ('/home/dwi/ownCloud/Anlabshared/Dwi/researchdata/170711SpinSeebeck110Diamond/ESR/1to3.5GHzH90OeDT5K.csv')
Data7.columns = ['Intensity']
Data8 = pd.read_csv ('/home/dwi/ownCloud/Anlabshared/Dwi/researchdata/170711SpinSeebeck110Diamond/ESR/1to3.5GHzH105OeDT5K.csv')
Data8.columns = ['Intensity']
Data9 = pd.read_csv ('/home/dwi/ownCloud/Anlabshared/Dwi/researchdata/170711SpinSeebeck110Diamond/ESR/1to3.5GHzH120OeDT5K.csv')
Data9.columns = ['Intensity']
Data10 = pd.read_csv ('/home/dwi/ownCloud/Anlabshared/Dwi/researchdata/170711SpinSeebeck110Diamond/ESR/1to3.5GHzH135OeDT5K.csv')
Data10.columns = ['Intensity']
Data11 = pd.read_csv ('/home/dwi/ownCloud/Anlabshared/Dwi/researchdata/170711SpinSeebeck110Diamond/ESR/1to3.5GHzH150OeDT5K.csv')
Data11.columns = ['Intensity']
Data12 = pd.read_csv ('/home/dwi/ownCloud/Anlabshared/Dwi/researchdata/170711SpinSeebeck110Diamond/ESR/1to3.5GHzH165OeDT5K.csv')
Data12.columns = ['Intensity']
Data13 = pd.read_csv ('/home/dwi/ownCloud/Anlabshared/Dwi/researchdata/170711SpinSeebeck110Diamond/ESR/1to3.5GHzH180OeDT5K.csv')
Data13.columns = ['Intensity']
Data14 = pd.read_csv ('/home/dwi/ownCloud/Anlabshared/Dwi/researchdata/170711SpinSeebeck110Diamond/ESR/1to3.5GHzH195OeDT5K.csv')
Data14.columns = ['Intensity']
Data15 = pd.read_csv ('/home/dwi/ownCloud/Anlabshared/Dwi/researchdata/170711SpinSeebeck110Diamond/ESR/1to3.5GHzH210OeDT5K.csv')
Data15.columns = ['Intensity']
Data16 = pd.read_csv ('/home/dwi/ownCloud/Anlabshared/Dwi/researchdata/170711SpinSeebeck110Diamond/ESR/1to3.5GHzH225OeDT5K.csv')
Data16.columns = ['Intensity']
Data17 = pd.read_csv ('/home/dwi/ownCloud/Anlabshared/Dwi/researchdata/170711SpinSeebeck110Diamond/ESR/1to3.5GHzH240OeDT5K.csv')
Data17.columns = ['Intensity']
Data18 = pd.read_csv ('/home/dwi/ownCloud/Anlabshared/Dwi/researchdata/170711SpinSeebeck110Diamond/ESR/1to3.5GHzH255OeDT5K.csv')
Data18.columns = ['Intensity']
Data19 = pd.read_csv ('/home/dwi/ownCloud/Anlabshared/Dwi/researchdata/170711SpinSeebeck110Diamond/ESR/1to3.5GHzH270OeDT5K.csv')
Data19.columns = ['Intensity']
Data20 = pd.read_csv ('/home/dwi/ownCloud/Anlabshared/Dwi/researchdata/170711SpinSeebeck110Diamond/ESR/1to3.5GHzH285OeDT5K.csv')
Data20.columns = ['Intensity']
Data21 = pd.read_csv ('/home/dwi/ownCloud/Anlabshared/Dwi/researchdata/170711SpinSeebeck110Diamond/ESR/1to3.5GHzH300OeDT5K.csv')
Data21.columns = ['Intensity']
#Data22 = pd.read_csv ('/home/dwi/ownCloud/Anlabshared/Dwi/researchdata/161229 LSEE (111) Magnet Field/2016_12_29_LSSE_(111)_FieldPerpdToTgrad_21K_Mf=0.0_Theta=0.0_Phi=0.0.csv')
#Data22.columns = ['Microwave frequency', 'Intensity']
#Data23 = pd.read_csv ('/home/dwi/ownCloud/Anlabshared/Dwi/researchdata/161229 LSEE (111) Magnet Field/2016_12_29_LSSE_(111)_FieldPerpdToTgrad_22K_Mf=0.0_Theta=0.0_Phi=0.0.csv')
#Data23.columns = ['Microwave frequency', 'Intensity']
#Data24 = pd.read_csv ('/home/dwi/ownCloud/Anlabshared/Dwi/researchdata/170103 LSSE (110) Diamond/2017_01_03_LSSE_(110)Diamond_FieldPerpdTgrad_18K_Mf=0.0_Theta=0.0_Phi=90.0.csv')
#Data24.columns = ['Microwave frequency', 'Intensity']
#Data25 = pd.read_csv ('/home/dwi/ownCloud/Anlabshared/Dwi/researchdata/170103 LSSE (110) Diamond/2017_01_03_LSSE_(110)Diamond_FieldPerpdTgrad_19K_Mf=0.0_Theta=0.0_Phi=90.0.csv')
#Data25.columns = ['Microwave frequency', 'Intensity']
#Data26 = pd.read_csv ('/home/dwi/ownCloud/Anlabshared/Dwi/researchdata/170103 LSSE (110) Diamond/2017_01_03_LSSE_(110)Diamond_FieldPerpdTgrad_20K_Mf=0.0_Theta=0.0_Phi=90.0.csv')
#Data26.columns = ['Microwave frequency', 'Intensity']

f = np.linspace(1e9, 3.5e9, 303)
y = np.array([f])#Data2['Microwave frequency'])
H = np.arange(0, 302, 15)
x = np.array([H])
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
#i22 = np.array(Data22['Intensity']/Data22['Intensity'].max())
#i23 = np.array(Data23['Intensity']/Data23['Intensity'].max())
#i24 = np.array(Data24['Intensity']/Data24['Intensity'].max())
#i25 = np.array(Data25['Intensity']/Data25['Intensity'].max())
#i26 = np.array(Data26['Intensity']/Data26['Intensity'].max())
offset1 = 1.0 - np.mean (i1[0:20])
io1 = i1 + offset1
offset2 = 1.0 - np.mean (i2[0:20])
io2 = i2 + offset2
offset3 = 1.0 - np.mean (i3[0:20])
io3 = i3 + offset3
offset4 = 1.0 - np.mean (i4[0:20])
io4 = i4 + offset4
offset5 = 1.0 - np.mean (i5[0:20])
io5 = i5 + offset5
offset6 = 1.0 - np.mean (i6[0:20])
io6 = i6 + offset6
offset7 = 1.0 - np.mean (i7[0:20])
io7 = i7 + offset7
offset8 = 1.0 - np.mean (i8[0:20])
io8 = i8 + offset8
offset9 = 1.0 - np.mean (i9[0:20])
io9 = i9 + offset9
offset10 = 1.0 - np.mean (i10[0:20])
io10 = i10 + offset10
offset11 = 1.0 - np.mean (i11[0:20])
io11 = i11 + offset11
offset12 = 1.0 - np.mean (i12[0:20])
io12 = i12 + offset12
offset13 = 1.0 - np.mean (i13[0:20])
io13 = i13 + offset13
offset14 = 1.0 - np.mean (i14[0:20])
io14 = i14 + offset14
offset15 = 1.0 - np.mean (i15[0:20])
io15 = i15 + offset15
offset16 = 1.0 - np.mean (i16[0:20])
io16 = i16 + offset16
offset17 = 1.0 - np.mean (i17[0:20])
io17 = i17 + offset17
offset18 = 1.0 - np.mean (i18[0:20])
io18 = i18 + offset18
offset19 = 1.0 - np.mean (i19[0:20])
io19 = i19 + offset19
offset20 = 1.0 - np.mean (i20[0:20])
io20 = i20 + offset20
offset21 = 1.0 - np.mean (i21[0:20])
io21 = i21 + offset21
I = np.vstack((io1, io2, io3, io4, io5, io6, io7, io8, io9, io10, io11, io12, io13, io14, io15, io16, io17, io18, io19, io20, io21))#, i19, i20, i21, i22, i23))
#I = I.transpose()
Im = np.array(I)
Y, X = np.meshgrid(y, x)

pl.rc('text', usetex=True)
pl.rc('font', family='serif', size = 22)
pl.xlabel('H (Oe)', fontsize=22)
pl.ylabel(r'f (GHz)', fontsize = 22)
#pl.xlim(2.6e9, 3.15e9)
pl.pcolormesh(X, Y, Im)
pl.colorbar()
#pl.title(r'Temperature Gradient Dependence of NV ESR on YIG with No External Field')
pl.savefig(path + r"ESR_110_DT5K.png", format = 'png', bbox_inches = 'tight', frameon = False)
pl.show()
#print (f)
print (I)
####################PRINTING INTENSITY OF DIPS##########################
#print (np.mean(i1[0:100])-i1[148:247].min())
#print (np.mean(i2[0:100])-i2[148:247].min())
#print (np.mean(i3[0:100])-i3[148:247].min())
#print (np.mean(i4[0:100])-i4[148:247].min())
#print (np.mean(i5[0:100])-i5[148:247].min())
#print (np.mean(i6[0:100])-i6[148:247].min())
#print (np.mean(i7[0:100])-i7[148:247].min())
#print (np.mean(i8[0:100])-i8[148:247].min())
#print (np.mean(i9[0:100])-i9[148:247].min())
#print (np.mean(i10[0:100])-i10[148:247].min())
#print (np.mean(i11[0:100])-i11[148:247].min())
#print (np.mean(i12[0:100])-i12[148:247].min())
#print (np.mean(i13[0:100])-i13[148:247].min())
#print (np.mean(i14[0:100])-i14[148:247].min())
#print (np.mean(i15[0:100])-i15[148:247].min())
#print (np.mean(i16[0:100])-i16[148:247].min())
#print (np.mean(i17[0:100])-i17[148:247].min())
#print (np.mean(i18[0:100])-i18[148:247].min())
#print (np.mean(i19[0:100])-i19[148:247].min())
#print (np.mean(i20[0:100])-i20[148:247].min())
#print (np.mean(i21[0:100])-i21[148:247].min())
#print (np.mean(i22[0:100])-i22[148:247].min())
#print (np.mean(i23[0:100])-i23[148:247].min())
#print (np.mean(i24[0:100])-i24[226:264].min())
#print (np.mean(i25[0:100])-i25[226:264].min())
#print (np.mean(i26[0:100])-i26[226:264].min())

#d1 = i1[148:247].min()
#d2 = i2[148:247].min()
#d3 = i3[148:247].min()
#d4 = i4[148:247].min()
#d5 = i5[148:247].min()
#d6 = i6[148:247].min()
#d7 = i7[148:247].min()
#d8 = i8[148:247].min()
#d9 = i9[148:247].min()
#d10 = i10[148:247].min()
#d11 = i11[148:247].min()
#d12 = i12[148:247].min()
#d13 = i13[148:247].min()
#d14 = i14[247:499].min()
#d15 = i15[247:499].min()
#d16 = i16[247:499].min()
#d17 = i17[247:499].min()
#d18 = i18[247:499].min()
#d19 = i19[247:499].min()
#d20 = i20[247:499].min()
#d21 = i21[247:499].min()
#d22 = i22[247:499].min()
#d23 = i23[247:499].min()

#D = np.array([ d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, d12])#, d13])#, d14, d15, d16, d17, d18, d19, d20, d21, d22, d23])
#R = Data1['Intensity'].index
#print(R)
#print(Data1.where(Data1 = Q))

#############PRINTING FREQUENCIES OF DIPS###############################
#for index, item in enumerate(i1):
    #if item == d1:
       # print(Data1['Microwave frequency'].iloc[index])

#for index, item in enumerate(i2):
#    if item == d2:
#        print(Data2['Microwave frequency'].iloc[index])   

#for index, item in enumerate(i3):
#    if item == d3:
#        print(Data3['Microwave frequency'].iloc[index])

#for index, item in enumerate(i4):
#    if item == d4:
#        print(Data4['Microwave frequency'].iloc[index])

#for index, item in enumerate(i5):
#    if item == d5:
#        print(Data5['Microwave frequency'].iloc[index])

#for index, item in enumerate(i6):
#    if item == d6:
#        print(Data6['Microwave frequency'].iloc[index])

#for index, item in enumerate(i7):
#    if item == d7:
#       print(Data7['Microwave frequency'].iloc[index])

#for index, item in enumerate(i8):
#    if item == d8:
#        print(Data8['Microwave frequency'].iloc[index])

#for index, item in enumerate(i9):
#    if item == d9:
#        print(Data9['Microwave frequency'].iloc[index])

#for index, item in enumerate(i10):
#    if item == d10:
#        print(Data10['Microwave frequency'].iloc[index])

#for index, item in enumerate(i11):
#    if item == d11:
#        print(Data11['Microwave frequency'].iloc[index])

#for index, item in enumerate(i12):
#    if item == d12:
#        print(Data12['Microwave frequency'].iloc[index])

#for index, item in enumerate(i13):
 #   if item == d13:
  #      print(Data13['Microwave frequency'].iloc[index])

#for index, item in enumerate(i14):
#    if item == d14:
#        print(Data14['Microwave frequency'].iloc[index])

#for index, item in enumerate(i15):
#    if item == d15:
#        print(Data15['Microwave frequency'].iloc[index])

#for index, item in enumerate(i16):
#    if item == d16:
#        print(Data16['Microwave frequency'].iloc[index])

#for index, item in enumerate(i17):
#    if item == d17:
#        print(Data17['Microwave frequency'].iloc[index])

#for index, item in enumerate(i18):
 #   if item == d18:
  #      print(Data18['Microwave frequency'].iloc[index])

#for index, item in enumerate(i19):
 #   if item == d19:
  #      print(Data19['Microwave frequency'].iloc[index])

#for index, item in enumerate(i20):
 #   if item == d20:
  #      print(Data20['Microwave frequency'].iloc[index])

#for index, item in enumerate(i21):
 #   if item == d21:
  #      print(Data21['Microwave frequency'].iloc[index])

#for index, item in enumerate(i22):
 #   if item == d22:
  #      print(Data22['Microwave frequency'].iloc[index])

#for index, item in enumerate(i23):
 #   if item == d23:
  #      print(Data23['Microwave frequency'].iloc[index])

#print (d1[1])
#print(D)
#print(I)
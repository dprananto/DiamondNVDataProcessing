import numpy as np
import matplotlib.pyplot as pl
import pandas as pd
from pandas import DataFrame, read_csv
Data1 = pd.read_csv ('/home/dwi/ownCloud/Anlabshared/Dwi/Data/160506/2016_05_06_Yellow_Diamond_10.csv')
Data1.columns = ['Microwave Frequency [Hz]', 'Counts']
Data2 = pd.read_csv ('/home/dwi/ownCloud/Anlabshared/Dwi/Data/160506/2016_05_06_Yellow_Diamond_11.csv')
Data2.columns = ['Microwave Frequency [Hz]', 'Counts']
Data3 = pd.read_csv ('/home/dwi/ownCloud/Anlabshared/Dwi/Data/160506/2016_05_06_Yellow_Diamond_12.csv')
Data3.columns = ['Microwave Frequency [Hz]', 'Counts']
Data4 = pd.read_csv ('/home/dwi/ownCloud/Anlabshared/Dwi/Data/160506/2016_05_06_Yellow_Diamond_13.csv')
Data4.columns = ['Microwave Frequency [Hz]', 'Counts']
Data5 = pd.read_csv ('/home/dwi/ownCloud/Anlabshared/Dwi/Data/160506/2016_05_06_Yellow_Diamond_14.csv')
Data5.columns = ['Microwave Frequency [Hz]', 'Counts']
Data6 = pd.read_csv ('/home/dwi/ownCloud/Anlabshared/Dwi/Data/160506/2016_05_06_Yellow_Diamond_15.csv')
Data6.columns = ['Microwave Frequency [Hz]', 'Counts']
Data7 = pd.read_csv ('/home/dwi/ownCloud/Anlabshared/Dwi/Data/160506/2016_05_06_Yellow_Diamond_16.csv')
Data7.columns = ['Microwave Frequency [Hz]', 'Counts']
Data8 = pd.read_csv ('/home/dwi/ownCloud/Anlabshared/Dwi/Data/160506/2016_05_06_Yellow_Diamond_17.csv')
Data8.columns = ['Microwave Frequency [Hz]', 'Counts']
Data9 = pd.read_csv ('/home/dwi/ownCloud/Anlabshared/Dwi/Data/160506/2016_05_06_Yellow_Diamond_18.csv')
Data9.columns = ['Microwave Frequency [Hz]', 'Counts']
Data10 = pd.read_csv ('/home/dwi/ownCloud/Anlabshared/Dwi/Data/160506/2016_05_06_Yellow_Diamond_19.csv')
Data10.columns = ['Microwave Frequency [Hz]', 'Counts']
Data11 = pd.read_csv ('/home/dwi/ownCloud/Anlabshared/Dwi/Data/160506/2016_05_06_Yellow_Diamond_20.csv')
Data11.columns = ['Microwave Frequency [Hz]', 'Counts']
Concat = pd.concat ([Data1, Data2['Counts']], axis=1)
A, B, C, D, E, F, G, H, I, J, K = Data1['Counts'], Data2['Counts'], Data3['Counts'], Data4['Counts'], Data5['Counts'], Data6['Counts'], Data7['Counts'], Data8['Counts'], Data9['Counts'], Data10['Counts'], Data11['Counts']
X1 = Data2['Microwave Frequency [Hz]']
X2 = Data3['Microwave Frequency [Hz]']
#print Data1[0:100], Data2
#pl.plot (X2, F, label='0.620 mW')
#pl.plot (X2, E, label='0.546 mW')
#pl.plot (X2, D, label='0.374 mW')
#pl.plot (X2, C, label='0.297 mW')
#pl.plot (X1, B, label='0.228 mW')
#pl.plot (X1, A, label='0.079 mW')
#pl.plot (X2, G, label='0.049 mW')
#pl.plot (X2, H, label='0.044 mW')
#pl.plot (X2, I, label='0.030 mW')
#pl.plot (X2, J, label='0.024 mW')
#pl.plot (X2, K, label='0.018 mW')
#pl.xlabel("Microwave Frequency [GHz]")
#pl.ylabel("Photon Counts per Second")
pl.legend()

#Data1.plot('Microwave Frequency [Hz]', 'Counts', color='Red')
#Data2.plot(x = 'Microwave Frequency [Hz]', y= 'Counts')
#Data3.plot(x = 'Microwave Frequency [Hz]', y= 'Counts')
#Data4.plot(x = 'Microwave Frequency [Hz]', y= 'Counts')
#Data5.plot(x = 'Microwave Frequency [Hz]', y= 'Counts')
#Data6.plot(x = 'Microwave Frequency [Hz]', y= 'Counts')
#Data7.plot(x = 'Microwave Frequency [Hz]', y= 'Counts')
#Data8.plot(x = 'Microwave Frequency [Hz]', y= 'Counts')
#Data9.plot(x = 'Microwave Frequency [Hz]', y= 'Counts')
#Data10.plot(x = 'Microwave Frequency [Hz]', y= 'Counts')
#Data11.plot(x = 'Microwave Frequency [Hz]', y= 'Counts')


#Signal to noise ratio calculation
N1 = Data1[0:99]
F1=N1['Counts'].mean()
ESR1=Data1['Counts'].min()
Noise1= N1['Counts'].max() - N1['Counts'].min()
S1= F1 - ESR1
SNR1=S1/Noise1
print N1
print 'Properties of Data1:'
print 'Fuorescence signal=' + str(F1)
print 'ESR minimum signal=' + str(ESR1)
print 'Noise signal=' + str(Noise1)
print 'SNR=' + str(SNR1)
#N1.plot('Microwave Frequency [Hz]', 'Counts')

N2 = Data2[0:99]
F2=N2['Counts'].mean()
ESR2=Data2['Counts'].min()
Noise2= N2['Counts'].max() - N2['Counts'].min()
S2= F2 - ESR2
SNR2=S2/Noise2
#print N1
print 'Properties of Data2:'
print 'Fuorescence signal=' + str(F2)
print 'ESR minimum signal=' + str(ESR2)
print 'Noise signal=' + str(Noise2)
print 'SNR=' + str(SNR2)
#N1.plot('Microwave Frequency [Hz]', 'Counts')

N3 = Data3[0:99]
F3=N3['Counts'].mean()
ESR3=Data3['Counts'].min()
Noise3= N3['Counts'].max() - N3['Counts'].min()
S3= F3 - ESR3
SNR3=S3/Noise3
#print N3
print 'Properties of Data3:'
print 'Fuorescence signal=' + str(F3)
print 'ESR minimum signal=' + str(ESR3)
print 'Noise signal=' + str(Noise3)
print 'SNR=' + str(SNR3)
#N1.plot('Microwave Frequency [Hz]', 'Counts')

N4 = Data4[0:99]
F4=N4['Counts'].mean()
ESR4=Data4['Counts'].min()
Noise4= N4['Counts'].max() - N4['Counts'].min()
S4= F4 - ESR4
SNR4=S4/Noise4
#print N3
print 'Properties of Data4:'
print 'Fuorescence signal=' + str(F4)
print 'ESR minimum signal=' + str(ESR4)
print 'Noise signal=' + str(Noise4)
print 'SNR=' + str(SNR4)
#N1.plot('Microwave Frequency [Hz]', 'Counts')

N5 = Data5[0:99]
F5=N5['Counts'].mean()
ESR5=Data5['Counts'].min()
Noise5= N5['Counts'].max() - N5['Counts'].min()
S5= F5 - ESR5
SNR5=S5/Noise5
#print N3
print 'Properties of Data5:'
print 'Fuorescence signal=' + str(F5)
print 'ESR minimum signal=' + str(ESR5)
print 'Noise signal=' + str(Noise5)
print 'SNR=' + str(SNR5)
#N1.plot('Microwave Frequency [Hz]', 'Counts')

N6 = Data6[0:99]
F6=N6['Counts'].mean()
ESR6=Data6['Counts'].min()
Noise6= N6['Counts'].max() - N6['Counts'].min()
S6= F6 - ESR6
SNR6=S6/Noise6
#print N3
print 'Properties of Data6:'
print 'Fuorescence signal=' + str(F6)
print 'ESR minimum signal=' + str(ESR6)
print 'Noise signal=' + str(Noise6)
print 'SNR=' + str(SNR6)
#N1.plot('Microwave Frequency [Hz]', 'Counts')

N7 = Data7[0:99]
F7=N7['Counts'].mean()
ESR7=Data7['Counts'].min()
Noise7= N7['Counts'].max() - N7['Counts'].min()
S7= F7 - ESR7
SNR7=S7/Noise7
#print N3
print 'Properties of Data7:'
print 'Fuorescence signal=' + str(F7)
print 'ESR minimum signal=' + str(ESR7)
print 'Noise signal=' + str(Noise7)
print 'SNR=' + str(SNR7)
#N1.plot('Microwave Frequency [Hz]', 'Counts')

N8 = Data8[0:99]
F8=N8['Counts'].mean()
ESR8=Data8['Counts'].min()
Noise8= N8['Counts'].max() - N8['Counts'].min()
S8= F8 - ESR8
SNR8=S8/Noise8
#print N3
print 'Properties of Data8:'
print 'Fuorescence signal=' + str(F8)
print 'ESR minimum signal=' + str(ESR8)
print 'Noise signal=' + str(Noise8)
print 'SNR=' + str(SNR8)
#N1.plot('Microwave Frequency [Hz]', 'Counts')

N9 = Data9[0:99]
F9=N9['Counts'].mean()
ESR9=Data9['Counts'].min()
Noise9= N9['Counts'].max() - N9['Counts'].min()
S9= F9 - ESR9
SNR9=S9/Noise9
#print N3
print 'Properties of Data9:'
print 'Fuorescence signal=' + str(F9)
print 'ESR minimum signal=' + str(ESR9)
print 'Noise signal=' + str(Noise9)
print 'SNR=' + str(SNR9)
#N1.plot('Microwave Frequency [Hz]', 'Counts')

N10 = Data10[0:99]
F10=N10['Counts'].mean()
ESR10=Data10['Counts'].min()
Noise10= N10['Counts'].max() - N10['Counts'].min()
S10= F10 - ESR10
SNR10=S10/Noise10
#print N3
print 'Properties of Data10:'
print 'Fuorescence signal=' + str(F10)
print 'ESR minimum signal=' + str(ESR10)
print 'Noise signal=' + str(Noise10)
print 'SNR=' + str(SNR10)
#N1.plot('Microwave Frequency [Hz]', 'Counts')

N11 = Data11[0:99]
F11=N11['Counts'].mean()
ESR11=Data11['Counts'].min()
Noise11= N11['Counts'].max() - N11['Counts'].min()
S11= F11 - ESR11
SNR11=S11/Noise11
#print N3
print 'Properties of Data11:'
print 'Fuorescence signal=' + str(F11)
print 'ESR minimum signal=' + str(ESR11)
print 'Noise signal=' + str(Noise11)
print 'SNR=' + str(SNR11)
#N1.plot('Microwave Frequency [Hz]', 'Counts')

LSNR = pd.DataFrame({'Laser Power [mW]': [0.079, 0.228, 0.297, 0.374, 0.546, .62, .049, .044, .030, .024, .018], 'SNR' : [SNR1, SNR2, SNR3, SNR4, SNR5, SNR6, SNR7, SNR8, SNR9, SNR10, SNR11]})
LS = pd.DataFrame({'Laser Power [mW]': [0.079, 0.228, 0.297, 0.374, 0.546, .62, .049, .044, .030, .024, .018], 'ESR Signal' : [S1, S2, S3, S4, S5, S6, S7, S8, S9, S10, S11]})
fit = np.polyfit(LSNR['Laser Power [mW]'], LSNR['SNR'], 2) 
p = np.poly1d(fit)
x= np.linspace(0.018, 0.62, 100)
print LSNR
print LS
pl.scatter(LSNR['Laser Power [mW]'], LSNR['SNR'])
pl.plot(x, p(x), color = 'Red')
pl.xlabel('Laser Power [mW]')
pl.ylabel('SNR')
pl.title('Laser Power Dependent of SNR')

#pl.scatter(LS['Laser Power [mW]'], LS['ESR Signal'], color = 'Red')
#pl.xlabel('Laser Power [mW]')
#pl.ylabel('ESR Signal [cps]')
#pl.title('Laser Power Dependent of ESR Signal')

pl.show()
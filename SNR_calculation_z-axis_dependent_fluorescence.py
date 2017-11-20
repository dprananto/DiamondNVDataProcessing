import matplotlib.pyplot as pl
import pandas as pd
Data1 = pd.read_csv('/home/dwi/ownCloud/Anlabshared/Dwi/Data/160509/2016_05_09_YD_zDependent_390.09um.csv')
Data1.columns  = ['Microwave Frequency [Hz]', 'Counts']
Data2 = pd.read_csv('/home/dwi/ownCloud/Anlabshared/Dwi/Data/160509/2016_05_09_YD_zDependent_397.29um.csv')
Data2.columns  = ['Microwave Frequency [Hz]', 'Counts']
Data3 = pd.read_csv('/home/dwi/ownCloud/Anlabshared/Dwi/Data/160509/2016_05_09_YD_zDependent_432.37um.csv')
Data3.columns  = ['Microwave Frequency [Hz]', 'Counts']

Data1.plot('Microwave Frequency [Hz]', 'Counts')
Data2.plot('Microwave Frequency [Hz]', 'Counts')
#pl.show()
#Signal to Noise Ratio Calculation
N1 = Data1[0:99]
F1=N1['Counts'].mean()
ESR1=Data1['Counts'].min()
Noise1= N1['Counts'].max() - N1['Counts'].min()
S1= F1 - ESR1
SNR1=S1/Noise1
#print N1
print 'Properties of Data 1 z = 390.09um:'
print 'Fuorescence signal=' + str(F1)
print 'ESR minimum signal=' + str(ESR1)
print 'ESR signal=' + str(S1)
print 'Noise signal=' + str(Noise1)
print 'SNR=' + str(SNR1)

N2 = Data2[0:99]
F2=N2['Counts'].mean()
ESR2=Data2['Counts'].min()
Noise2= N2['Counts'].max() - N2['Counts'].min()
S2= F2 - ESR2
SNR2=S2/Noise2
#print N3
print 'Properties of Data 2 z = 397.29um:'
print 'Fuorescence signal=' + str(F2)
print 'ESR minimum signal=' + str(ESR2)
print 'ESR signal=' + str(S2)
print 'Noise signal=' + str(Noise2)
print 'SNR=' + str(SNR2)

N3 = Data3[0:99]
F3=N3['Counts'].mean()
ESR3=Data3['Counts'].min()
Noise3= N3['Counts'].max() - N3['Counts'].min()
SNR3=(F3 - ESR3)/Noise2
#print N3
print 'Properties of Data 3 z = 432.37um:'
print 'Fuorescence signal=' + str(F3)
print 'ESR minimum signal=' + str(ESR3)
print 'Noise signal=' + str(Noise3)
print 'SNR=' + str(SNR3)
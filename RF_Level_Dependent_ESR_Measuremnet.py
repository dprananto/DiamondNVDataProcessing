import matplotlib.pyplot as pl
import pandas as pd
Data1 = pd.read_csv('/home/dwi/ownCloud/Anlabshared/Dwi/Data/160509/2016_05_09_YD_RFLevel_dependent_16dbm.csv')
Data1.columns  = ['Microwave Frequency [Hz]', 'Counts']
Data2 = pd.read_csv('/home/dwi/ownCloud/Anlabshared/Dwi/Data/160509/2016_05_09_YD_RFLevel_dependent_16.5dbm.csv')
Data2.columns  = ['Microwave Frequency [Hz]', 'Counts']
Data3 = pd.read_csv('/home/dwi/ownCloud/Anlabshared/Dwi/Data/160509/2016_05_09_YD_RFLevel_dependent_17dbm.csv')
Data3.columns  = ['Microwave Frequency [Hz]', 'Counts']
Data4 = pd.read_csv('/home/dwi/ownCloud/Anlabshared/Dwi/Data/160509/2016_05_09_YD_RFLevel_dependent_17.5dbm.csv')
Data4.columns  = ['Microwave Frequency [Hz]', 'Counts']
Data5 = pd.read_csv('/home/dwi/ownCloud/Anlabshared/Dwi/Data/160509/2016_05_09_YD_RFLevel_dependent_18dbm.csv')
Data5.columns  = ['Microwave Frequency [Hz]', 'Counts']

pl.plot(Data1['Microwave Frequency [Hz]'], Data1['Counts'], label='16 dBm')
pl.plot(Data2['Microwave Frequency [Hz]'], Data2['Counts'], label='16.5 dBm')
pl.plot(Data3['Microwave Frequency [Hz]'], Data3['Counts'], label='17 dBm')
pl.plot(Data4['Microwave Frequency [Hz]'], Data4['Counts'], label='17.5 dBm')
pl.plot(Data5['Microwave Frequency [Hz]'], Data5['Counts'], label='18 dBm')
pl.legend()
pl.title('RF Level Dependent Measurement')
pl.xlabel('Microwave Frequency [GHz]')
pl.ylabel('Photon Counts per Second')
pl.show()

#Signal to Noise Ratio Calculation
N1 = Data1[0:99]
F1=N1['Counts'].mean()
ESR1=Data1['Counts'].min()
Noise1= N1['Counts'].max() - N1['Counts'].min()
S1= F1 - ESR1
SNR1=S1/Noise1
#print N1
print 'Properties of Data 1 16dBm:'
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
#print N1
print 'Properties of Data 2 16.5dBm:'
print 'Fuorescence signal=' + str(F2)
print 'ESR minimum signal=' + str(ESR2)
print 'ESR signal=' + str(S2)
print 'Noise signal=' + str(Noise2)
print 'SNR=' + str(SNR2)

N3 = Data3[0:99]
F3=N3['Counts'].mean()
ESR3=Data3['Counts'].min()
Noise3= N3['Counts'].max() - N3['Counts'].min()
S3= F3 - ESR3
SNR3=S3/Noise3
#print N1
print 'Properties of Data 3 17dBm:'
print 'Fuorescence signal=' + str(F3)
print 'ESR minimum signal=' + str(ESR3)
print 'ESR signal=' + str(S3)
print 'Noise signal=' + str(Noise3)
print 'SNR=' + str(SNR3)

N4 = Data4[0:99]
F4=N4['Counts'].mean()
ESR4=Data4['Counts'].min()
Noise4= N4['Counts'].max() - N4['Counts'].min()
S4= F4 - ESR4
SNR4=S4/Noise4
#print N1
print 'Properties of Data 4 17.5dBm:'
print 'Fuorescence signal=' + str(F4)
print 'ESR minimum signal=' + str(ESR4)
print 'ESR signal=' + str(S4)
print 'Noise signal=' + str(Noise4)
print 'SNR=' + str(SNR4)

N5 = Data5[0:99]
F5=N5['Counts'].mean()
ESR5=Data5['Counts'].min()
Noise5= N5['Counts'].max() - N5['Counts'].min()
S5= F5 - ESR5
SNR5=S5/Noise5
#print N1
print 'Properties of Data 5 18dBm:'
print 'Fuorescence signal=' + str(F5)
print 'ESR minimum signal=' + str(ESR5)
print 'ESR signal=' + str(S5)
print 'Noise signal=' + str(Noise5)
print 'SNR=' + str(SNR5)
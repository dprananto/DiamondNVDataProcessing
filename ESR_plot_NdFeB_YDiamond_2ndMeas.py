import matplotlib.pyplot as pl
import numpy as np
import pandas as pd

Data1 = pd.read_csv ('/home/dwi/ownCloud/Anlabshared/Dwi/Data and Analysis/160525/2016_05_25_magparticle_.post1.csv')
Data1.columns = ['Microwave Frequency [Hz]', 'Counts']
Data2 = pd.read_csv ('/home/dwi/ownCloud/Anlabshared/Dwi/Data and Analysis/160525/2016_05_25_magparticle_.post2.csv')
Data2.columns = ['Microwave Frequency [Hz]', 'Counts']
Data3 = pd.read_csv ('/home/dwi/ownCloud/Anlabshared/Dwi/Data and Analysis/160525/2016_05_25_magparticle_.post3.csv')
Data3.columns = ['Microwave Frequency [Hz]', 'Counts']
Data4 = pd.read_csv ('/home/dwi/ownCloud/Anlabshared/Dwi/Data and Analysis/160525/2016_05_25_magparticle_.post4.csv')
Data4.columns = ['Microwave Frequency [Hz]', 'Counts']
Data5 = pd.read_csv ('/home/dwi/ownCloud/Anlabshared/Dwi/Data and Analysis/160525/2016_05_25_magparticle_.post5.csv')
Data5.columns = ['Microwave Frequency [Hz]', 'Counts']


Data1_norm = Data1['Counts'] / Data1['Counts'].max()
Data2_norm = Data2['Counts'] / Data2['Counts'].max()
Data3_norm = Data3['Counts'] / Data3['Counts'].max()
Data4_norm = Data4['Counts'] / Data4['Counts'].max()
Data5_norm = Data5['Counts'] / Data5['Counts'].max()

#print Data1_norm
norm_Data1 = pd.concat ([Data1['Microwave Frequency [Hz]'], Data1_norm], axis = 1)
norm_Data1.columns = ['Microwave Frequency [Hz]', 'Counts']
norm_Data2 = pd.concat ([Data2['Microwave Frequency [Hz]'], Data2_norm], axis = 1)
norm_Data2.columns = ['Microwave Frequency [Hz]', 'Counts']
norm_Data3 = pd.concat ([Data3['Microwave Frequency [Hz]'], Data3_norm], axis = 1)
norm_Data3.columns = ['Microwave Frequency [Hz]', 'Counts']
norm_Data4 = pd.concat ([Data4['Microwave Frequency [Hz]'], Data4_norm], axis = 1)
norm_Data4.columns = ['Microwave Frequency [Hz]', 'Counts']
norm_Data5 = pd.concat ([Data5['Microwave Frequency [Hz]'], Data5_norm], axis = 1)
norm_Data5.columns = ['Microwave Frequency [Hz]', 'Counts']

#print norm_Data1

pl.plot (norm_Data1['Microwave Frequency [Hz]'], norm_Data1['Counts'], linewidth = 2, label ='Position 1')
#pl.plot (norm_Data2['Microwave Frequency [Hz]'], norm_Data2['Counts'], linewidth = 2, label ='Position 2')
#pl.plot (norm_Data3['Microwave Frequency [Hz]'], norm_Data3['Counts'], linewidth = 2, label ='Position 3')
#pl.plot (norm_Data4['Microwave Frequency [Hz]'], norm_Data4['Counts'], linewidth = 2, label ='Position 4')
#pl.plot (norm_Data5['Microwave Frequency [Hz]'], norm_Data5['Counts'], linewidth = 2,  label ='Position 5')

pl.legend()
pl.title('ESR Spectra of Diamond NV Ensemble Around Magnetic Particle of NdFeB')
pl.xlabel ('Microwave Frequency [GHz]')
pl.ylabel ('Photon Counts per Second (norm.)')
pl.show()
#print Data1
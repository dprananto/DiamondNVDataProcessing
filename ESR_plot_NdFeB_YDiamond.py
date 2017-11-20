import matplotlib.pyplot as pl
import numpy as np
import pandas as pd

Data1 = pd.read_csv ('/home/dwi/ownCloud/Anlabshared/Dwi/Data/160518/2016_05_18_MagParticleNdFeB_YDiamond_pos10.csv')
Data1.columns = ['Microwave Frequency [Hz]', 'Counts']
Data2 = pd.read_csv ('/home/dwi/ownCloud/Anlabshared/Dwi/Data/160518/2016_05_18_MagParticleNdFeB_YDiamond_pos11.csv')
Data2.columns = ['Microwave Frequency [Hz]', 'Counts']
Data3 = pd.read_csv ('/home/dwi/ownCloud/Anlabshared/Dwi/Data/160518/2016_05_18_MagParticleNdFeB_YDiamond_pos12.csv')
Data3.columns = ['Microwave Frequency [Hz]', 'Counts']
Data4 = pd.read_csv ('/home/dwi/ownCloud/Anlabshared/Dwi/Data/160518/2016_05_18_MagParticleNdFeB_YDiamond_pos12#1.csv')
Data4.columns = ['Microwave Frequency [Hz]', 'Counts']
Data5 = pd.read_csv ('/home/dwi/ownCloud/Anlabshared/Dwi/Data/160518/2016_05_18_MagParticleNdFeB_YDiamond_pos14.csv')
Data5.columns = ['Microwave Frequency [Hz]', 'Counts']
Data6 = pd.read_csv ('/home/dwi/ownCloud/Anlabshared/Dwi/Data/160518/2016_05_18_MagParticleNdFeB_YDiamond_pos15.csv')
Data6.columns = ['Microwave Frequency [Hz]', 'Counts']
Data7 = pd.read_csv ('/home/dwi/ownCloud/Anlabshared/Dwi/Data/160518/2016_05_18_MagParticleNdFeB_YDiamond_pos16.csv')
Data7.columns = ['Microwave Frequency [Hz]', 'Counts']
Data8 = pd.read_csv ('/home/dwi/ownCloud/Anlabshared/Dwi/Data/160518/2016_05_18_MagParticleNdFeB_YDiamond_pos17.csv')
Data8.columns = ['Microwave Frequency [Hz]', 'Counts']
Data9 = pd.read_csv ('/home/dwi/ownCloud/Anlabshared/Dwi/Data/160518/2016_05_18_MagParticleNdFeB_YDiamond_pos18.csv')
Data9.columns = ['Microwave Frequency [Hz]', 'Counts']
Data10 = pd.read_csv ('/home/dwi/ownCloud/Anlabshared/Dwi/Data/160518/2016_05_18_MagParticleNdFeB_YDiamond_pos19.csv')
Data10.columns = ['Microwave Frequency [Hz]', 'Counts']
Data11 = pd.read_csv ('/home/dwi/ownCloud/Anlabshared/Dwi/Data/160518/2016_05_18_MagParticleNdFeB_YDiamond_pos20.csv')
Data11.columns = ['Microwave Frequency [Hz]', 'Counts']
Data12 = pd.read_csv ('/home/dwi/ownCloud/Anlabshared/Dwi/Data/160518/2016_05_18_MagParticleNdFeB_YDiamond_pos21.csv')
Data12.columns = ['Microwave Frequency [Hz]', 'Counts']
Data13 = pd.read_csv ('/home/dwi/ownCloud/Anlabshared/Dwi/Data/160518/2016_05_18_MagParticleNdFeB_YDiamond_pos22.csv')
Data13.columns = ['Microwave Frequency [Hz]', 'Counts']
Data14 = pd.read_csv ('/home/dwi/ownCloud/Anlabshared/Dwi/Data/160518/2016_05_18_MagParticleNdFeB_YDiamond_pos23.csv')
Data14.columns = ['Microwave Frequency [Hz]', 'Counts']
Data15 = pd.read_csv ('/home/dwi/ownCloud/Anlabshared/Dwi/Data/160518/2016_05_18_MagParticleNdFeB_YDiamond_pos24.csv')
Data15.columns = ['Microwave Frequency [Hz]', 'Counts']

Data1_norm = Data1['Counts'] / Data1['Counts'].max()
Data2_norm = Data2['Counts'] / Data2['Counts'].max()
Data3_norm = Data3['Counts'] / Data3['Counts'].max()
Data4_norm = Data4['Counts'] / Data4['Counts'].max()
Data5_norm = Data5['Counts'] / Data5['Counts'].max()
Data6_norm = Data6['Counts'] / Data6['Counts'].max()
Data7_norm = Data7['Counts'] / Data7['Counts'].max()
Data8_norm = Data8['Counts'] / Data8['Counts'].max()
Data9_norm = Data9['Counts'] / Data9['Counts'].max()
Data10_norm = Data10['Counts'] / Data10['Counts'].max()
Data11_norm = Data11['Counts'] / Data11['Counts'].max()
Data12_norm = Data12['Counts'] / Data12['Counts'].max()
Data13_norm = Data13['Counts'] / Data13['Counts'].max()
Data14_norm = Data14['Counts'] / Data14['Counts'].max()
Data15_norm = Data15['Counts'] / Data15['Counts'].max()
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
norm_Data6 = pd.concat ([Data6['Microwave Frequency [Hz]'], Data6_norm], axis = 1)
norm_Data6.columns = ['Microwave Frequency [Hz]', 'Counts']
norm_Data7 = pd.concat ([Data7['Microwave Frequency [Hz]'], Data7_norm], axis = 1)
norm_Data7.columns = ['Microwave Frequency [Hz]', 'Counts']
norm_Data8 = pd.concat ([Data8['Microwave Frequency [Hz]'], Data8_norm], axis = 1)
norm_Data8.columns = ['Microwave Frequency [Hz]', 'Counts']
norm_Data9 = pd.concat ([Data9['Microwave Frequency [Hz]'], Data9_norm], axis = 1)
norm_Data9.columns = ['Microwave Frequency [Hz]', 'Counts']
norm_Data10 = pd.concat ([Data10['Microwave Frequency [Hz]'], Data10_norm], axis = 1)
norm_Data10.columns = ['Microwave Frequency [Hz]', 'Counts']
norm_Data11 = pd.concat ([Data11['Microwave Frequency [Hz]'], Data11_norm], axis = 1)
norm_Data11.columns = ['Microwave Frequency [Hz]', 'Counts']
norm_Data12 = pd.concat ([Data12['Microwave Frequency [Hz]'], Data12_norm], axis = 1)
norm_Data12.columns = ['Microwave Frequency [Hz]', 'Counts']
norm_Data13 = pd.concat ([Data13['Microwave Frequency [Hz]'], Data13_norm], axis = 1)
norm_Data13.columns = ['Microwave Frequency [Hz]', 'Counts']
norm_Data14 = pd.concat ([Data14['Microwave Frequency [Hz]'], Data14_norm], axis = 1)
norm_Data14.columns = ['Microwave Frequency [Hz]', 'Counts']
norm_Data15 = pd.concat ([Data15['Microwave Frequency [Hz]'], Data15_norm], axis = 1)
norm_Data15.columns = ['Microwave Frequency [Hz]', 'Counts']
#print norm_Data1

#pl.plot (norm_Data1['Microwave Frequency [Hz]'], norm_Data1['Counts'], label ='Position 1')
pl.plot (norm_Data2['Microwave Frequency [Hz]'], norm_Data2['Counts'], linewidth = 2, label ='Position 2')
#pl.plot (norm_Data3['Microwave Frequency [Hz]'], norm_Data3['Counts'], label ='Position 3')
pl.plot (norm_Data4['Microwave Frequency [Hz]'], norm_Data4['Counts'], linewidth = 2, label ='Position 4')
pl.plot (norm_Data5['Microwave Frequency [Hz]'], norm_Data5['Counts'], linewidth = 2,  label ='Position 5')
pl.plot (norm_Data6['Microwave Frequency [Hz]'], norm_Data6['Counts'], linewidth = 2, label ='Position 6')
#pl.plot (norm_Data7['Microwave Frequency [Hz]'], norm_Data7['Counts'], label ='Position 7')
pl.plot (norm_Data8['Microwave Frequency [Hz]'], norm_Data8['Counts'], linewidth = 2, label ='Position 8')
pl.plot (norm_Data9['Microwave Frequency [Hz]'], norm_Data9['Counts'], linewidth = 2, label ='Position 9')
pl.plot (norm_Data10['Microwave Frequency [Hz]'], norm_Data10['Counts'], linewidth = 2, label ='Position 10')
#pl.plot (norm_Data11['Microwave Frequency [Hz]'], norm_Data11['Counts'], linewidth = 2, label ='Position 11')
#pl.plot (norm_Data12['Microwave Frequency [Hz]'], norm_Data12['Counts'], linewidth = 2, label ='Position 12')
#pl.plot (norm_Data13['Microwave Frequency [Hz]'], norm_Data13['Counts'], linewidth = 2, label ='Position 13')
#pl.plot (norm_Data14['Microwave Frequency [Hz]'], norm_Data14['Counts'], linewidth = 2, label ='Position 14')
#pl.plot (norm_Data15['Microwave Frequency [Hz]'], norm_Data15['Counts'], linewidth = 2, label ='Position 15')
pl.legend()
pl.title('ESR Spectra of Diamond NV Ensemble Around Magnetic Particle of NdFeB')
pl.xlabel ('Microwave Frequency [GHz]')
pl.ylabel ('Photon Counts per Second (norm.)')
pl.show()
#print Data1
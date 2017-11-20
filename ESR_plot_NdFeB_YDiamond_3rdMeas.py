import pandas as pd
import numpy as np
import matplotlib.pyplot as pl


Data1 = pd.read_csv ('/home/dwi/ownCloud/Anlabshared/Dwi/Data and Analysis/160627/2016_06_27_MagPar_pos14.csv')
Data1.columns = ['Microwave freq. (GHz)', 'Intensity (cps)']
x = Data1['Microwave freq. (GHz)']
y = Data1['Intensity (cps)']/Data1['Intensity (cps)'].max()


Data2 = pd.read_csv ('/home/dwi/ownCloud/Anlabshared/Dwi/Data and Analysis/160627/2016_06_27_MagPar_pos15.csv')
Data2.columns = ['Microwave freq. (GHz)', 'Intensity (cps)']
x2 = Data2['Microwave freq. (GHz)']
y2 = Data2['Intensity (cps)']/Data2['Intensity (cps)'].max()

Data3 = pd.read_csv ('/home/dwi/ownCloud/Anlabshared/Dwi/Data and Analysis/160627/2016_06_27_MagPar_pos16.csv')
Data3.columns = ['Microwave freq. (GHz)', 'Intensity (cps)']
x3 = Data3['Microwave freq. (GHz)']
y3 = Data3['Intensity (cps)']/Data3['Intensity (cps)'].max()


Data4 = pd.read_csv ('/home/dwi/ownCloud/Anlabshared/Dwi/Data and Analysis/160627/2016_06_27_MagPar_pos8.csv')
Data4.columns = ['Microwave freq. (GHz)', 'Intensity (cps)']
x4 = Data4['Microwave freq. (GHz)']
y4 = Data4['Intensity (cps)']/Data4['Intensity (cps)'].max()

Data5 = pd.read_csv ('/home/dwi/ownCloud/Anlabshared/Dwi/Data and Analysis/160627/2016_06_27_MagPar_pos4.csv')
Data5.columns = ['Microwave freq. (GHz)', 'Intensity (cps)']
x5 = Data5['Microwave freq. (GHz)']
y5 = Data5['Intensity (cps)']/Data5['Intensity (cps)'].max()


Data6 = pd.read_csv ('/home/dwi/ownCloud/Anlabshared/Dwi/Data and Analysis/160627/2016_06_27_MagPar_pos20.csv')
Data6.columns = ['Microwave freq. (GHz)', 'Intensity (cps)']
x6 = Data6['Microwave freq. (GHz)']
y6 = Data6['Intensity (cps)']/Data6['Intensity (cps)'].max()

Data7 = pd.read_csv ('/home/dwi/ownCloud/Anlabshared/Dwi/Data and Analysis/160627/2016_06_27_MagPar_pos9.csv')
Data7.columns = ['Microwave freq. (GHz)', 'Intensity (cps)']
x7 = Data7['Microwave freq. (GHz)']
y7 = Data7['Intensity (cps)']/Data7['Intensity (cps)'].max()


Data8 = pd.read_csv ('/home/dwi/ownCloud/Anlabshared/Dwi/Data and Analysis/160627/2016_06_27_MagPar_pos21.csv')
Data8.columns = ['Microwave freq. (GHz)', 'Intensity (cps)']
x8 = Data8['Microwave freq. (GHz)']
y8 = Data8['Intensity (cps)']/Data8['Intensity (cps)'].max()

pl.rc('text', usetex=True)
pl.rc('font', family='serif')

pl.plot (x, y, lw = 2,label='pos 14')
pl.plot (x2, y2, lw = 2, label='pos 15')
pl.plot (x3, y3, lw = 2, label='pos 16')
#pl.plot (x4, y4, lw = 2, label='pos 8')
#pl.plot (x5, y5, lw = 2, label='pos 4')
#pl.plot (x6, y6, lw = 2, label='pos 20')
#pl.plot (x7, y7, lw = 2,  label='pos 9')
#pl.plot (x8, y8, lw = 2, label='pos 21')
pl.xlabel('Microwave frequency (GHz)', fontsize = 18)
pl.ylabel('Intensity (norm.)', fontsize = 18)
pl.legend(loc=4)
pl.show()
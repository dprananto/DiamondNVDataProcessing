import pandas as pd
import matplotlib.pyplot as pl

Data1 = pd.read_csv('/home/dwi/ownCloud/Anlabshared/Dwi/Data and Analysis/160718/2016_07_18_5.csv')
Data1.columns = ['Microwave Frequency', 'Intensity']
x1 = Data1['Microwave Frequency']
y1 = Data1['Intensity']/Data1['Intensity'].max()
Data2 = pd.read_csv('/home/dwi/ownCloud/Anlabshared/Dwi/Data and Analysis/160718/2016_07_18_6.csv')
Data2.columns = ['Microwave Frequency', 'Intensity']
x2 = Data2['Microwave Frequency']
y2 = (Data2['Intensity']/Data2['Intensity'].max()) - 0.02
Data3 = pd.read_csv('/home/dwi/ownCloud/Anlabshared/Dwi/Data and Analysis/160715/2016_07_15_8-5.csv')
Data3.columns = ['Microwave Frequency', 'Intensity']
x3 = Data3['Microwave Frequency']
y3 = (Data3['Intensity']/Data3['Intensity'].max()) - 0.04
Data4 = pd.read_csv('/home/dwi/ownCloud/Anlabshared/Dwi/Data and Analysis/160718/2016_07_18_4.csv')
Data4.columns = ['Microwave Frequency', 'Intensity']
x4 = Data4['Microwave Frequency']
y4 = Data4['Intensity']/Data4['Intensity'].max() - 0.06
Data5 = pd.read_csv('/home/dwi/ownCloud/Anlabshared/Dwi/Data and Analysis/160715/2016_07_15_5-5.csv')
Data5.columns = ['Microwave Frequency', 'Intensity']
x5 = Data5['Microwave Frequency']
y5 = Data5['Intensity']/Data5['Intensity'].max() - 0.08
Data6 = pd.read_csv('/home/dwi/ownCloud/Anlabshared/Dwi/Data and Analysis/160715/2016_07_15_5-6.csv')
Data6.columns = ['Microwave Frequency', 'Intensity']
x6 = Data6['Microwave Frequency']
y6 = Data6['Intensity']/Data6['Intensity'].max() - 0.1
Data7 = pd.read_csv('/home/dwi/ownCloud/Anlabshared/Dwi/Data and Analysis/160715/2016_07_15_5-7.csv')
Data7.columns = ['Microwave Frequency', 'Intensity']
x7 = Data7['Microwave Frequency']
y7 = Data7['Intensity']/Data7['Intensity'].max() - 0.02
Data8 = pd.read_csv('/home/dwi/ownCloud/Anlabshared/Dwi/Data and Analysis/160715/2016_07_15_5-8.csv')
Data8.columns = ['Microwave Frequency', 'Intensity'] 
x8 = Data8['Microwave Frequency']
y8 = Data8['Intensity']/Data8['Intensity'].max() - 0.14
Data9 = pd.read_csv('/home/dwi/ownCloud/Anlabshared/Dwi/Data and Analysis/160715/2016_07_15_5-9.csv')
Data9.columns = ['Microwave Frequency', 'Intensity']
x9 = Data9['Microwave Frequency']
y9 = Data9['Intensity']/Data9['Intensity'].max() - 0.16
Data10 = pd.read_csv('/home/dwi/ownCloud/Anlabshared/Dwi/Data and Analysis/160705/2016_07_05_pos10.csv')
Data10.columns = ['Microwave Frequency', 'Intensity']
x10 = Data10['Microwave Frequency']
y10 = Data10['Intensity']/Data10['Intensity'].max() - 0.08

pl.rc('text', usetex=True)
pl.rc('font', family='serif')
pl.xlabel(r"Microwave Frequency (GHz)", fontsize = 18)
pl.ylabel(r"Intensity (arb.)", fontsize = 18)
pl.plot(x1, y1, lw =1.5, label = r"$[111]$")
pl.plot(x2, y2, lw =1.5, label = r"$[1\bar{1}\bar{1} ]$")
pl.plot(x3, y3, lw =1.5, label = r"$[\bar{1}\bar{1}1]$")
pl.plot(x4, y4, lw =1.5, label = r"$[\bar{1}1\bar{1}]$")
#pl.plot(x5, y5, lw =1.5, label = r"5")
#pl.plot(x6, y6, lw =1.5, label = r"6")
#pl.plot(x7, y7, lw =1.5, label = r"7")
#pl.plot(x8, y8, lw =1.5, label = r"8")
#pl.plot(x9, y9, lw =1.5, label = r"9")
#pl.plot(x10, y10, lw =1.5, label = r"Position 10")
pl.legend(loc = 8)
pl.show()

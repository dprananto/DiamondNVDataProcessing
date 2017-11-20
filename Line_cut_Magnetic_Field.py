import pandas as pd
import matplotlib.pyplot as pl

Data = pd.read_csv('/home/dwi/ownCloud/Anlabshared/Dwi/Data and Analysis/160705/linecut_Magnetic_Field_Magnitude.csv')
Data.columns = ['Distance','Magnetic Field']
x = Data['Distance']
y = Data['Magnetic Field']

pl.rc('text', usetex=True)
pl.rc('font', family='serif')
pl.xlabel(r"Distance ($\mu \text{m}$)", fontsize=18)
pl.ylabel(r"Magnetic Field (mT)", fontsize=18)
pl.plot(x, y, 'ro')
pl.show()

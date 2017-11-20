import numpy as np
import matplotlib.pyplot as pl
import pandas as pd
from pandas import read_csv
Data = pd.read_csv('/home/dwi/ownCloud/Anlabshared/Dwi/Data and Analysis/160525/Graphs/LineProfile_NdFeB_on_YDNVC_2,894GHz.csv')
Data.columns = ['X', 'Y']

pl.rc('text', usetex=True)
pl.rc('font', family='serif')

pl.plot(Data['X'], Data['Y'])
pl.xlabel(r'x ($\mu$m)', fontsize=18)
pl.ylabel(r'Intensity diff.(arb.)', fontsize=18)
pl.show()
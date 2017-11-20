import pandas as pd
import matplotlib.pyplot as pl
import numpy as np
Data = pd.read_csv('/home/dwi/ownCloud/Anlabshared/Dwi/Data and Analysis/160511/2016_05_11_nanodiamond_OD16_100times_4_best.csv')
Data.columns = ['Microwave Frequency [GHz]', 'Intensity']
fit = np.polyfit(Data['Microwave Frequency [GHz]'], Data['Intensity'], 9)
p=np.poly1d(fit)
x=np.linspace(2801500000, 2950000000, 200)
pl.plot (Data['Microwave Frequency [GHz]'], Data['Intensity'])
#pl.scatter (Data['Microwave Frequency [GHz]'], Data['Photons Count per Second'])
pl.plot(x, p(x), linewidth = 1.5, color="Red")
pl.title ('ESR Spectrum of Nanodiamond NV Center')
pl.xlabel('Microwave Frequency [Ghz]')
pl.ylabel ('Intensity [cps.]')
pl.show()
print Data
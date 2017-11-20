import numpy as np
import matplotlib.pyplot as pl
import pandas as pd
from pandas import DataFrame, read_csv
Z = pd.read_csv('/home/dwi/ownCloud/Anlabshared/Dwi/Data and Analysis/160610/z-dependent_upsidedown_YellowDiamond.csv')
Z.columns = ['z-axis', 'counts']
#fit = np.polyfit(Z['z-axis'], Z['counts'], 30)
#p=np.poly1d(fit)
#x=np.linspace(1175, 1215, 100)
#print fit 
#Z.plot(x="z-axis", y='counts')
pl.scatter(Z['z-axis'], Z['counts'])
#pl.plot(x, p(x))
pl.title('z-axis Dependent of Diamond NV Center Fluorescence (upsidedown-side)')
pl.xlabel('Microscope z-axis [um]')
pl.ylabel('Fluorescence Counts per Second')
pl.show()
print Z
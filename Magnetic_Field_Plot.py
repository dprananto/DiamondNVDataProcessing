# -*- coding: utf-8 -*-
"""
Created on Wed Aug 24 18:22:16 2016

@author: dwi
"""
import pandas as pd
import numpy as np
#import matplotlib.pyplot as pl

fi = input ("File path to open: ")
Data = pd.read_csv (fi)
ttl = input (r"Graph title: ")
#pl.title (ttl)
Data.columns = ['Microwave frequency', 'Intensity']
x = np.array(Data['PosX']) #* (-1) - 10
y = np.array(Data['B'])
#X = np.linspace(1, 16, 100)
#Y = (1/(X**2)) * (1/1000)

#pl.rc('text', usetex=True)
#pl.rc('font', family='serif', size = 20)
#pl.plot (x, y, 'ro', lw = 2, label = 'Experiment', color = 'b')
#pl.plot (X, Y, color = 'r')
#pl.xlabel(r"Position ($\mu$m)", fontsize = 18)
#pl.ylabel('Magnetic Field (T)', fontsize = 18)

#pl.legend(loc=0, ncol=3)
#pl.show()
print (x)
print (y)
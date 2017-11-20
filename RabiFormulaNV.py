# -*- coding: utf-8 -*-
"""
Created on Wed Aug  2 19:42:58 2017

@author: dwi
"""

import numpy as np
from numpy import cos
import matplotlib.pyplot as pl

x0 = np.array([2.65e6, 2.65e6]) #Initial coefficient

t = np.linspace (0, 10e-6, 300)
def func (t, O_bar, O):
    return ((O_bar/O)**2 * (1 - cos(O*t))/2)

pl.rc('text', usetex=True)
pl.rc('font', family='serif', size = 20)    
pl.xlabel (r'$\tau$ ($\mu$s)', fontsize = 20.0)
pl.ylabel (r'$P_{|-1\rangle}$', fontsize = 20.0)
print(func(t, *x0))
pl.plot(t*1e6, func(t, *x0))
pl.show()

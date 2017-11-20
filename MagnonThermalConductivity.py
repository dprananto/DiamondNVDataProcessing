# -*- coding: utf-8 -*-
"""
Created on Tue Jan 31 17:51:15 2017

@author: dwi
"""

#YIG Magnon and Phonon thermal Characteristics based on Sanders and Walton Paper PRB 15, 1489
import matplotlib.pyplot as pl
import numpy as np

T = np.linspace (300, 350, 100) #Temperature
Km = 0.0183 * T ** 2.0 #Magnon thermal conductivity
Kp = 0.0104 * T **3.0 #Phonon thermal conductivity
Cm = 3.34e-6 * T ** (3/2) #Magnon specific heat
Cp = 1.52e-6 * T ** 3.0


pl.plot (T, Cm, label = 'Cm')
pl.plot (T, Cp, label = 'Cp')
#pl.plot (T, Km, label = 'Km')
#pl.plot (T, Kp, label = 'Kp')
pl.legend()
pl.show()
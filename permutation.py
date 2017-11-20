# -*- coding: utf-8 -*-
"""
Created on Fri Jul 29 11:47:38 2016

@author: dwi
"""

import numpy as np
import itertools

for c in itertools.permutations (np.linspace (0.0, 360, (36*2 + 1)), 2):
    print (c)

#print (c[0])
#com = np.linspace (0.0, 360, 37)
#print (com)
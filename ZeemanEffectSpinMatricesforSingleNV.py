# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 22:58:48 2017

@author: Dwi Prananto
"""

from physicsconstants import h_bar, G_e, k_B, u0, g_e, u_B, h
import numpy as np
from numpy import exp, cos, sin, pi, sqrt, mat, dot, linalg
import matplotlib.pyplot as pl

#B = np.linspace(0.0, 20.0e-3, 100)
D = 2.873e9 #Zero-field splitting
E = 1e6 #strain
t = 0.0 * pi /180.0 #\theta
p = 0.0 * pi / 180.0 #\varphi
B_0 = np.linspace(0.0, 30e-3, 100)

#Pauli matrices
s_x = mat([[0, 1.0, 0], [1.0, 0, 1.0], [0, 1.0, 0]]) / sqrt(2)
s_y = mat([[0, 1.0, 0], [-1.0, 0, 1.0], [0, -1.0, 0]]) / sqrt(2) * 1j
s_z = mat([[1.0, 0, 0], [0, 0, 0], [0, 0, -1.0]])


fplist = []
fmlist = []
for B in B_0:
    B_x = sin(t) * cos(p) * B
    B_y = sin(t) * sin(p) * B
    B_z = cos(t) * B
    Hp = h_bar * D * (s_z**2) + g_e * u_B * B_z * s_z + g_e * u_B * (B_x * s_x + B_y * s_y) + h_bar * E * (s_x ** 2 - s_y ** 2)
    Hm = h_bar * D * (s_z**2) - g_e * u_B * B_z * s_z + g_e * u_B * (B_x * s_x + B_y * s_y) + h_bar * E * (s_x ** 2 - s_y ** 2)
    EigP = linalg.eigvalsh(Hp / h_bar)#Eigenvalues of Hamiltonian matrix (Hermitian)
    EigM = linalg.eigvalsh(Hm / h_bar)#Eigenvalues of Hamiltonian matrix (Hermitian)
    f_p = EigP[2]
    f_m = EigM[1]
    fplist.append(f_p)
    fmlist.append(f_m)

print(Hp)    

#print(fplist)
pl.plot(B_0, fplist)
pl.plot(B_0, fmlist)
#pl.ylim(2.7e9, 3.05e9)
pl.show()
print(Hp.conjugate())
#print(f_m)
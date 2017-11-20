# -*- coding: utf-8 -*-
"""
Created on Wed Feb 15 09:20:04 2017

@author: dwi
"""

from physicsconstants import h_bar, G_e, k_B, u0, g_e, u_B, h
import numpy as np
from numpy import exp, cos, sin, pi, sqrt, mat, dot, linalg
import matplotlib.pyplot as pl

#B = np.linspace(0.0, 20.0e-3, 100)
Diamond = 111
D = 2.873e9 #Zero-field splitting
E = 1e6 #strain
t = 109.5 * pi /180.0 #\theta
#t2 = 45.0 * pi /180.0 #second \theta
p = 0.0 * pi / 180.0 #\varphi
B_0 = np.linspace(0.0, 60e-3, 100)

#Pauli spin matrices
s_x = mat([[0, 1.0, 0], [1.0, 0, 1.0], [0, 1.0, 0]]) / sqrt(2)
s_y = mat([[0, 1.0, 0], [-1.0, 0, 1.0], [0, -1.0, 0]]) / sqrt(2) * 1j
s_z = mat([[1.0, 0, 0], [0, 0, 0], [0, 0, -1.0]])
s_p = mat([[0, 1, 0], [0, 0, 1], [0, 0, 0]]) * sqrt(2)
s_m = mat([[0, 0, 0], [1, 0, 0], [0, 1, 0]])


if Diamond == 110:
    NV_1 = mat([sin(90.0 * pi / 180.0) * cos (35.25 * pi / 180.0), sin(90.0 * pi /180.0) * sin(35.25 * pi /180.0), cos(90.0 * pi /180.0)]) / sqrt(3) #NV1 vector
    NV_1x = mat([sin(0.0 * pi / 180.0) * cos (0.0 * pi / 180.0), sin(0.0 * pi /180.0) * sin(0.0 * pi /180.0), cos(0.0 * pi /180.0)]) / sqrt(3)
    NV_1y = mat([sin(90.0 * pi / 180.0) * cos (305.25 * pi / 180.0), sin(90.0 * pi /180.0) * sin(305.25 * pi /180.0), cos(90.0 * pi /180.0)]) / sqrt(3)
    NV_2 = mat([sin(90.0 * pi / 180.0) * cos(144.75 * pi /180.0), sin(90.0 * pi /180.0) * sin(144.75 * pi /180.0), cos(90.0 * pi /180.0)]) / sqrt(3) #NV2 vector
    NV_3 = mat([sin(35.25 * pi /180.0) * cos(270.0 * pi /180.0), sin(35.25 * pi /180.0) * sin(270.0 * pi /180.0), cos(35.25 * pi /180.0)]) / sqrt(3) #NV3 vector
    NV_3x = mat([sin(125.25 * pi /180.0) * cos(270.0 * pi /180.0), sin(125.25 * pi /180.0) * sin(270.0 * pi /180.0), cos(125.25 * pi /180.0)]) / sqrt(3)
    NV_3y = mat([sin(90.0 * pi /180.0) * cos(360.0 * pi /180.0), sin(90.0 * pi /180.0) * sin(360.0 * pi /180.0), cos(90.0 * pi /180.0)]) / sqrt(3)
    NV_4 = mat([sin(144.75 * pi /180.0) * cos(270.0 * pi /180.0), sin(144.75 * pi /180.0) * sin(270.0 * pi /180.0), cos(144.75 * pi /180.0)]) /sqrt(3) #NV4 vector
    
if Diamond == 111:
    NV_1 = mat([sin(0.0 * pi / 180.0) * cos (0.0 * pi / 180.0), sin(0.0 * pi /180.0) * sin(0.0 * pi /180.0), cos(0.0 * pi /180.0)]) / sqrt(3) #NV1 vector
    NV_1x = mat([sin(90.0 * pi / 180.0) * cos (0.0 * pi / 180.0), sin(90.0 * pi /180.0) * sin(0.0 * pi /180.0), cos(90.0 * pi /180.0)]) / sqrt(3)
    NV_1y = mat([sin(90.0 * pi / 180.0) * cos (90.0 * pi / 180.0), sin(90.0 * pi /180.0) * sin(90.0 * pi /180.0), cos(90.0 * pi /180.0)]) / sqrt(3)
    NV_2 = mat([sin(109.5 * pi / 180.0) * cos(0.0 * pi /180.0), sin(109.5 * pi /180.0) * sin(0.0 * pi /180.0), cos(109.5 * pi /180.0)]) / sqrt(3) #NV2 vector
    NV_3 = mat([sin(109.5 * pi /180.0) * cos(120.0 * pi /180.0), sin(109.5 * pi /180.0) * sin(120.0 * pi /180.0), cos(109.5 * pi /180.0)]) / sqrt(3) #NV3 vector
    NV_4 = mat([sin(109.5 * pi /180.0) * cos(240.0 * pi /180.0), sin(109.5 * pi /180.0) * sin(240.0 * pi /180.0), cos(109.5 * pi /180.0)]) / sqrt(3) #NV4 vector
if Diamond == 100:
    #S_1 = mat([1, 1, 1])/ sqrt(3)
    #S_1 = mat([-1, 1, -1])/ sqrt(3)
    #S_1 = mat([-1, -1, 1])/ sqrt(3)
    #S_1 = mat([1, -1, -1])/ sqrt(3)
    NV_1 = mat([sin(54.75 * pi / 180.0) * cos (0.0 * pi / 180.0), sin(54.75 * pi /180.0) * sin(0.0 * pi /180.0), cos(54.75 * pi /180.0)]) #/ sqrt(3) #NV1 vector
    NV_2 = mat([sin(125.25 * pi / 180.0) * cos(90.0 * pi /180.0), sin(125.25 * pi /180.0) * sin(90.0 * pi /180.0), cos(125.25 * pi /180.0)]) #/ sqrt(3) #NV2 vector
    NV_3 = mat([sin(54.75 * pi /180.0) * cos(180.0 * pi /180.0), sin(54.75 * pi /180.0) * sin(180.0 * pi /180.0), cos(54.75 * pi /180.0)]) #/ sqrt(3) #NV3 vector
    NV_4 = mat([sin(125.25 * pi /180.0) * cos(270.0 * pi /180.0), sin(125.25 * pi /180.0) * sin(270.0 * pi /180.0), cos(125.25 * pi /180.0)]) #/ sqrt(3) #NV4 vector
#S = mat([[s_x], [s_y], [s_z]])


fplist1 = []
fmlist1 = []
#fplist1_2 = []
#fmlist1_2 = []
fplist2 = []
fmlist2 = []
fplist3 = []
fmlist3 = []
fplist4 = []
fmlist4 = []
for B in B_0:
    B_x = sin(t) * cos(p) * B
    B_y = sin(t) * sin(p) * B
    B_z = cos(t) * B
    #B_x2 = sin(t2) * cos(p) * B
    #B_y2 = sin(t2) * sin(p) * B
    #B_z2 = cos(t2) * B
    B_vec = mat([[B_x], [B_y], [B_z]])
    #B_vec2 = mat([[B_x2], [B_y2], [B_z2]])
    #Hp1 = (h * D * (s_z**2) + g_e * u_B * B_z * s_z + g_e * u_B * (B_x * s_x + B_y * s_y) + h * E * (s_x ** 2 - s_y ** 2)) 
    Hp1 = h * D * (s_z**2) + g_e * u_B * linalg.norm(dot(NV_1, B_vec)) * s_z + g_e * u_B * (B_x * s_x + B_y * s_y) + h * E * (s_x ** 2 - s_y ** 2)
    #Hm1 = (h * D * (s_z**2) - g_e * u_B * B_z * s_z + g_e * u_B * (B_x * s_x + B_y * s_y) + h * E * (s_x ** 2 - s_y ** 2)) 
    Hm1 = h * D * (s_z**2) - g_e * u_B * linalg.norm(dot(NV_1, B_vec)) * s_z - g_e * u_B * (linalg.norm(dot(NV_1x, B_vec)) * s_x + linalg.norm(dot(NV_1y, B_vec)) * s_y) + h * E * (s_x ** 2 - s_y ** 2)
    Hp2 = h * D * (s_z**2) + g_e * u_B * linalg.norm(dot(NV_2, B_vec)) * s_z + g_e * u_B * (linalg.norm(dot(NV_1x, B_vec)) * s_x + linalg.norm(dot(NV_1y, B_vec)) * s_y) + h * E * (s_x ** 2 - s_y ** 2)
    Hm2 = h * D * (s_z**2) - g_e * u_B * linalg.norm(dot(NV_2, B_vec)) * s_z - g_e * u_B * (B_x * s_x + B_y * s_y) + h * E * (s_x ** 2 - s_y ** 2)
    Hp3 = h * D * (s_z**2) + g_e * u_B * linalg.norm(dot(NV_3, B_vec)) * s_z + g_e * u_B * (B_x * s_x + B_y * s_y) + h * E * (s_x ** 2 - s_y ** 2)
    Hm3 = h * D * (s_z**2) - g_e * u_B * linalg.norm(dot(NV_3, B_vec)) * s_z - g_e * u_B * (B_x * s_x + B_y * s_y) + h * E * (s_x ** 2 - s_y ** 2)
    Hp4 = h * D * (s_z**2) + g_e * u_B * linalg.norm(dot(NV_4, B_vec)) * s_z + g_e * u_B * (B_x * s_x + B_y * s_y) + h * E * (s_x ** 2 - s_y ** 2)
    Hm4 = h * D * (s_z**2) - g_e * u_B * linalg.norm(dot(NV_4, B_vec)) * s_z - g_e * u_B * (B_x * s_x + B_y * s_y) + h * E * (s_x ** 2 - s_y ** 2)
    EigP1 = linalg.eigvalsh(Hp1 /h)#Eigenvalues of Hamiltonian matrix (Hermitian)
    EigM1 = linalg.eigvalsh(Hm1 /h)#Eigenvalues of Hamiltonian matrix (Hermitian)
    #EigP1_2 = linalg.eigvalsh(Hp1_2 / h)#Eigenvalues of Hamiltonian matrix (Hermitian)
    #EigM1_2 = linalg.eigvalsh(Hm1_2 / h)#Eigenvalues of Hamiltonian matrix (Hermitian)
    EigP2 = linalg.eigvalsh(Hp2 / h)#Eigenvalues of Hamiltonian matrix (Hermitian)
    EigM2 = linalg.eigvalsh(Hm2 / h)#Eigenvalues of Hamiltonian matrix (Hermitian)    
    EigP3 = linalg.eigvalsh(Hp3 / h)#Eigenvalues of Hamiltonian matrix (Hermitian)
    EigM3 = linalg.eigvalsh(Hm3 / h)#Eigenvalues of Hamiltonian matrix (Hermitian)
    EigP4 = linalg.eigvalsh(Hp4 / h)#Eigenvalues of Hamiltonian matrix (Hermitian)
    EigM4 = linalg.eigvalsh(Hm4 / h)#Eigenvalues of Hamiltonian matrix (Hermitian)    
    f_p1 = EigP1[2]
    f_m1 = EigM1[1]
    #f_p1_2 = EigP1_2[2]
    #f_m1_2 = EigM1_2[1]
    f_p2 = EigP2[2]
    f_m2 = EigM2[1]
    f_p3 = EigP3[2]
    f_m3 = EigM3[1]
    f_p4 = EigP4[2]
    f_m4 = EigM4[1]
    fplist1.append(f_p1)
    fmlist1.append(f_m1)
    #fplist1_2.append(f_p1_2)
    #fmlist1_2.append(f_m1_2)
    fplist2.append(f_p2)
    fmlist2.append(f_m2)
    fplist3.append(f_p3)
    fmlist3.append(f_m3)   
    fplist4.append(f_p4)
    fmlist4.append(f_m4)
    

#print(B_vec)
#print(NV_1)
#print(dot (NV_4, B_vec))    

#print(fplist)
pl.rc('text', usetex=True)
pl.rc('font', family='serif', size = 22) 
pl.xlabel(r'Magnetic field (T)', fontsize = 28)
pl.ylabel(r'$f$ (GHz)', fontsize = 28)
pl.plot(B_0, fplist1, lw = 2, color = 'r')
pl.plot(B_0, fmlist1, lw = 2, color = 'r', label = 'NV1')
#pl.plot(B_0, fplist1_2, lw = 2, color = 'b')
#pl.plot(B_0, fmlist1_2, lw = 2, color = 'b', label = r'$\theta$ = 60^o')
pl.plot(B_0, fplist2, lw = 2, color = 'green', label = 'NV2')
pl.plot(B_0, fmlist2, lw = 2, color = 'g',)
pl.plot(B_0, fplist3, lw = 2, color = 'blue', label = 'NV3')
pl.plot(B_0, fmlist3, lw = 2, color = 'b',)
#pl.plot(B_0, fplist4, lw = 2, color = 'black', label = 'NV4')
#pl.plot(B_0, fmlist4, lw = 2, color = 'k',)
pl.ylim(0.0e9, 5e9)
pl.legend(loc = 2, ncol= 2, fontsize = 14)
pl.show()

#print(Hp.conjugate())
#print(s_m)
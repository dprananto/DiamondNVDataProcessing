# -*- coding: utf-8 -*-
"""
Created on Thu Feb  9 16:51:47 2017

@author: dwi
"""
#Zeeman Effects on (111) diamond nv

from physicsconstants import h_bar, G_e, k_B, u0, g_e, u_B
import numpy as np
from numpy import exp, cos, sin, pi, sqrt, mat, dot, linalg
import matplotlib.pyplot as pl

Diamond = 111
t = 90.0 * pi /180.0 #\theta
p = 0.0 * pi / 180.0 #\varphi
D = 2.873e9 #Zero-field splitting
B_0 = np.linspace(0, 40.0e-3, 100) #External field
E = 100.0e3
#B_0 = 10.0e-3 #External field statics
#B_y = 5.0e-3

#w = G_e * B_0 #Magnon frequency

#T_m = np.linspace(370.0, 300.0, 100) #Magnon temperature
#n = 1 / (exp(h_bar * w / k_B * T_m) - 1.0)

#D = 2.25849832e9 + 5.66491395e6 * T_m - 1.71685005e4 * T_m**2.0 + 16.9451600 * T_m**3.0

B = mat([[sin(t) * cos(p)], [sin(t) * sin(p)], [cos(t)]]) / sqrt(3) #External field vector

#for i in B_0:
    #B_xyz = i * mat([[sin(t) * cos(p)], [sin(t) * cos(p)], [cos(t)]])
    
#B_m = linalg.norm(B_xyz)
    
#NV axes vector
if Diamond == 110:
    S_1 = mat([sin(90.0 * pi / 180.0) * cos (0.0 * pi / 180.0), sin(90.0 * pi /180.0) * sin(90.0 * pi /180.0), cos(0.0 * pi /180.0)]) / sqrt(3) #NV1 vector
    S_2 = mat([sin(90.0 * pi / 180.0) * cos(109.5 * pi /180.0), sin(90.0 * pi /180.0) * sin(109.5 * pi /180.0), cos(90.0 * pi /180.0)]) / sqrt(3) #NV2 vector
    S_3 = mat([sin(19.5 * pi /180.0) * cos(234.75 * pi /180.0), sin(19.5 * pi /180.0) * sin(234.75 * pi /180.0), cos(19.5 * pi /180.0)]) / sqrt(3) #NV3 vector
    S_4 = mat([sin(144.75 * pi /180.0) * cos(234.75 * pi /180.0), sin(144.75 * pi /180.0) * sin(234.75 * pi /180.0), cos(144.75 * pi /180.0)]) /sqrt(3) #NV4 vector
if Diamond == 111:
    S_1 = mat([sin(0.0 * pi / 180.0) * cos (0.0 * pi / 180.0), sin(0.0 * pi /180.0) * sin(0.0 * pi /180.0), cos(0.0 * pi /180.0)]) / sqrt(3) #NV1 vector
    S_2 = mat([sin(109.5 * pi / 180.0) * cos(0.0 * pi /180.0), sin(109.5 * pi /180.0) * sin(0.0 * pi /180.0), cos(109.5 * pi /180.0)]) / sqrt(3) #NV2 vector
    S_3 = mat([sin(109.5 * pi /180.0) * cos(120.0 * pi /180.0), sin(109.5 * pi /180.0) * sin(120.0 * pi /180.0), cos(109.5 * pi /180.0)]) / sqrt(3) #NV3 vector
    S_4 = mat([sin(109.5 * pi /180.0) * cos(240.0 * pi /180.0), sin(109.5 * pi /180.0) * sin(240.0 * pi /180.0), cos(109.5 * pi /180.0)]) / sqrt(3) #NV4 vector
if Diamond == 100:
    #S_1 = mat([1, 1, 1])/ sqrt(3)
    #S_1 = mat([-1, 1, -1])/ sqrt(3)
    #S_1 = mat([-1, -1, 1])/ sqrt(3)
    #S_1 = mat([1, -1, -1])/ sqrt(3)
    S_1 = mat([sin(54.75 * pi / 180.0) * cos (0.0 * pi / 180.0), sin(54.75 * pi /180.0) * sin(0.0 * pi /180.0), cos(54.75 * pi /180.0)]) #/ sqrt(3) #NV1 vector
    S_2 = mat([sin(125.25 * pi / 180.0) * cos(90.0 * pi /180.0), sin(125.25 * pi /180.0) * sin(90.0 * pi /180.0), cos(125.25 * pi /180.0)]) #/ sqrt(3) #NV2 vector
    S_3 = mat([sin(54.75 * pi /180.0) * cos(180.0 * pi /180.0), sin(54.75 * pi /180.0) * sin(180.0 * pi /180.0), cos(54.75 * pi /180.0)]) #/ sqrt(3) #NV3 vector
    S_4 = mat([sin(125.25 * pi /180.0) * cos(270.0 * pi /180.0), sin(125.25 * pi /180.0) * sin(270.0 * pi /180.0), cos(125.25 * pi /180.0)]) #/ sqrt(3) #NV4 vector

BS1 = dot(S_1, B) #Dot product B dot S1
BS2 = dot(S_2, B) #Dot product B dot S2
BS3 = dot(S_3, B) #Dot product B dot S3
BS4 = dot(S_4, B) #Dot product B dot S4
BS11 = BS1[0,0] #Take value of dot product
BS21 = BS2[0,0] #Take value of dot product
BS31 = BS3[0,0] #Take value of dot product
BS41 = BS4[0,0] #Take value of dot product

#Frequency as a function of magnon temperature
#f_p1 = D + G_e * (B_0 + B_y * n) * BS11
#f_p2 = D + G_e * (B_0 + B_y * n) * BS21
#f_p3 = D + G_e * (B_0 + B_y * n) * BS31
#f_p4 = D + G_e * (B_0 + B_y * n) * BS41
#f_m1 = D - G_e * (B_0 + B_y * n) * BS11
#f_m2 = D - G_e * (B_0 + B_y * n) * BS21
#f_m3 = D - G_e * (B_0 + B_y * n) * BS31
#f_m4 = D - G_e * (B_0 + B_y * n) * BS41 

#Frequency as a function of external field
f_p1 = D + G_e * B_0 * BS11 + E
f_p2 = D + G_e * B_0 * BS21 + E
f_p3 = D + G_e * B_0 * BS31 + E
f_p4 = D + G_e * B_0 * BS41 + E
f_m1 = D - G_e * B_0 * BS11 + E
f_m2 = D - G_e * B_0 * BS21 + E
f_m3 = D - G_e * B_0 * BS31 + E
f_m4 = D - G_e * B_0 * BS41 + E

#print (BS31)
pl.plot(B_0, f_p1, B_0, f_m1, color = 'red')
pl.plot(B_0, f_p2, B_0, f_m2, color = 'green') 
pl.plot(B_0, f_p3, B_0, f_m3, color = 'blue') 
pl.plot(B_0, f_p4, B_0, f_m4, color = 'black')

#pl.plot(T_m, f_p1, T_m, f_m1, color = 'red')
#pl.plot(T_m, f_p2, T_m, f_m2, color = 'green') 
#pl.plot(T_m, f_p3, T_m, f_m3, color = 'blue') 
#pl.plot(T_m, f_p4, T_m, f_m4, color = 'black')

pl.xlabel(r'Magnetic field (T)', fontsize = 22)
pl.ylabel(r'$\omega$(Hz)', fontsize = 22)
pl.ylim(1.5e9, 4.0e9)
#pl.xlim(0.0, 0.03)
pl.rc('text', usetex=True)
pl.rc('font', family='serif', size = 22)  
pl.show()
#print(B_xyz)
print(S_1)
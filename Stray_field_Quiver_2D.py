# -*- coding: utf-8 -*-
"""
Created on Mon Nov 14 09:58:54 2016
Ploting vector Field in 2D quiver plot
@author: dwi
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as pl
#import matplotlib.image as impng
#from mpl_toolkits.mplot3d import Axes3D
#from matplotlib.cbook import get_sample_data
#from matplotlib._png import read_png
from matplotlib.pyplot import cm
#from matplotlib.ticker import LinearLocator, FormatStrFormatter

#fig = pl.figure()
#ax = fig.add_subplot (111, projection='3d')

Data = pd.read_csv('/home/dwi/ownCloud/Anlabshared/Dwi/researchdata/160830 Magnetic particle stray mapping/scan 25 x 25 um/Analysis/analysis_selected.csv')
X = np.array([Data['PosX']])
Y = np.array([Data['PosY']])
#Z = np.array([Data['PosZ']])
#X, Y , Z = np.meshgrid(x, y, z)
#X = Data[0:10]['x']
#Y = Data[0:10]['y']
#Z = Data[0:10]['z']
B = Data['Bz'] #MAGNETIC FIELD
T = Data['T'] #THETA
P = Data['P'] #PHI

U = np.array(Data['Bx'])
V = np.array(Data['By'])
W = np.array(Data['Bz'])

length = np.sqrt(U**2 + V**2) * 10**3
UN = U/length
VN = V /length

XX, YY = np.meshgrid (X, Y)
#U = B * np.sin(T * np.pi/180.0) * np.cos(P * np.pi/180.0) #VECTOR COMPONENT-X
#V = B * np.sin(T * np.pi/180.0) * np.sin(P * np.pi/180.0) #VECTOR COMPONENT-Y
#W = B * np.cos(T * np.pi/180.0) #VECTOR COMPONENT-Z
#U, V, W = np.meshgrid (u, v, w)
 
#U = Data.iloc[0:10]['Bx']
#V = Data.iloc[0:10]['By']
#W = Data.iloc[0:10]['Bz']

#################### DRAW SPHERE ##############################################
#a = np.linspace (0.0, np.pi, 100)
#b = np.linspace (0.0, 2 * np.pi, 100)

#A, B = np.meshgrid (a, b)

#I = 12 * np.sin (A) * np.cos (B)
#J = 12 * np.sin (A) * np.sin (B)
#K = 12 * np.cos (A)

#fn = get_sample_data('/home/dwi/ownCloud/Anlabshared/Dwi/Paper writing/Stray Field Measurement/Figures/Fluorescence_image2.png', asfileobj=False)
#img = read_png (fn)

pl.gca().invert_yaxis()
pl.gca().invert_xaxis()
pl.rc('text', usetex=True)
pl.rc('font', family='serif', size = 22)

#ax.set_zlim (-25.0, 25.0)
#ax.set_xlim (-25.0, 25.0)
#ax.set_ylim (-25.0, 25.0)
#ax.plot_wireframe (I, J, K, rstride=9, cstride=9, color = 'b')
#ax.set_zlabel('z ($\mu$um)')
#ax.zaxis.set_major_locator(LinearLocator(10))
#ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
###############################################################################
#print (A)
#print (B)
#print (Z)
#print (U)
#print (V)
#print (W)

#stepX, stepY = 50./img.shape[0], 50./img.shape[1]
#C = np.arange(-25., 25., stepX)
#D = np.arange(-25., 25., stepY)
#C, D = np.meshgrid(C, D)

pl.figure()
pl.xlabel (r"x ($\mu$m)")
pl.ylabel (r"y ($\mu$m)")
pl.quiver(X, Y, UN, VN, length, cmap = cm.winter, headlength = 5)
clb = pl.colorbar()
clb.ax.set_title(r'$\textbf{B} (mT)$')
#pl.plot_surface(D, C, -15.0, rstride=6, cstride=6, facecolors=img)
pl.show()

#print (C)
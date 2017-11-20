import numpy as np
import pandas as pd
import matplotlib.pyplot as pl
#import matplotlib.image as impng
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.cbook import get_sample_data
from matplotlib._png import read_png
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from matplotlib.pyplot import cm

fig = pl.figure()
ax = fig.add_subplot (111, projection='3d')

savpath = r'/home/dwi/ownCloud/Anlabshared/Dwi/Projects/StrayFieldMeasurement/Figures/'
Data = pd.read_csv(r'/home/dwi/ownCloud/Anlabshared/Dwi/researchdata/160830 Magnetic particle stray mapping/scan 25 x 25 um/Analysis/analysis_selected.csv')
X = np.array([Data['PosX']])
Y = np.array([Data['PosY']])
Z = np.array([Data['PosZ']])
#X, Y , Z = np.meshgrid(x, y, z)
#X = Data[0:10]['x']
#Y = Data[0:10]['y']
#Z = Data[0:10]['z']
B = Data['B'] #MAGNETIC FIELD
T = Data['T'] #THETA
P = Data['P'] #PHI

U = np.array(Data['Bx'])
V = np.array(Data['By'])
W = np.array(Data['Bz'])
#U = B * np.sin(T * np.pi/180.0) * np.cos(P * np.pi/180.0) #VECTOR COMPONENT-X
#V = B * np.sin(T * np.pi/180.0) * np.sin(P * np.pi/180.0) #VECTOR COMPONENT-Y
#W = B * np.cos(T * np.pi/180.0) #VECTOR COMPONENT-Z
#U, V, W = np.meshgrid (u, v, w)
 
#U = Data.iloc[0:10]['Bx']
#V = Data.iloc[0:10]['By']
#W = Data.iloc[0:10]['Bz']

#################### DRAW SPHERE ##############################################
a = np.linspace (0.0, np.pi, 100)
b = np.linspace (0.0, 2 * np.pi, 100)

A, B = np.meshgrid (a, b)

I = 12 * np.sin (A) * np.cos (B) 
J = 12 * np.sin (A) * np.sin (B)
K = 12 * np.cos (A) + 12

fn = get_sample_data(r'/home/dwi/ownCloud/Anlabshared/Dwi/Projects/StrayFieldMeasurement/Figures/Fluorescence_image2.png', asfileobj=False)
img = read_png (fn)


pl.rc('text', usetex=True)
pl.rc('font', family='serif', size = 20)

ax.set_zlim (-10.0, 20.0)
ax.set_xlim (-26.0, 26.0)
ax.set_ylim (-26.0, 26.0)
ax.plot_wireframe (I, J, K, rstride=9, cstride=9, color = 'r')
#ax.zaxis.set_major_locator(LinearLocator(10))
#ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
###############################################################################
#print (A)
#print (B)
#print (Z)
#print (U)
#print (V)
#print (W)

stepX, stepY = 50./img.shape[0], 50./img.shape[1]
C = np.arange(-25., 25., stepX)
D = np.arange(-25., 25., stepY)
C, D = np.meshgrid(C, D)

#pl.xlabel (r"x ($\mu$m)")
#pl.ylabel (r"y ($\mu$m)")
#ax.set_zlabel('z ($\mu$m)')
ax.quiver(X, Y, Z, U, V, W, length = 5, color = 'c', arrow_length_ratio = .6, pivot = 'tail')
ax.plot_surface(D, C, -2, rstride=3, cstride=3, facecolors=img)
ax.view_init(40, 45)
pl.tight_layout(pad=0.0001, w_pad=10, h_pad=10)
pl.savefig(savpath + 'quiver3D.png', format = 'png', bbox_inches = 'tight', frameon = True)
pl.show()

#print (C)
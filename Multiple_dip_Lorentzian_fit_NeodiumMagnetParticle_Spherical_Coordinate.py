import pandas as pd
import numpy as np
import matplotlib.pyplot as pl
from scipy.optimize import leastsq
from scipy.signal import argrelmin

Data = pd.read_csv ('/home/dwi/ownCloud/Anlabshared/Dwi/Data and Analysis/160705/2016_07_05_pos7.csv')
Data.columns = ['Microwave frequency', 'Intensity']
x = np.array(Data['Microwave frequency'])
y = np.array(Data['Intensity']/Data['Intensity'].max())


############################## PEAK SEARCH #####################################
dipind = argrelmin(y, order=19, mode='clip')
ind = np.array(dipind[0])

xdip = x[dipind]
ydip = y[dipind]

############################## INITIAL PARAMETERS ##############################
p = [5.0*10**6, 5.0*10**6, 5.0*10**6, 5.0*10**6, 0.003, 0.003, 0.003, 0.003, 0.5*10**(-3), 90.0, 10.0]  
# [HWHM1-4, dip 1-4, B0, Theta, Phi]
############################# LORENTZIAN FUNCTION ##############################
def Lorentzian(x, p):
    num = p[0]**2
    num2 = p[1]**2
    num3= p[2]**2
    num4 = p[3]**2
    num5 = p[3]**2
    num6 = p[2]**2
    num7 = p[1]**2
    num8 = p[0]**2
    denum = ( x - (2.87*10**9 - (28*10**9 * p[8] * (1/np.sqrt(3)) * (np.sin (p[9]*np.pi/180.0) * np.cos(p[10]*np.pi/180.0) + np.sin(p[9]*np.pi/180.0) * np.sin(p[10]*np.pi/180.0) + np.cos(p[9]*np.pi/180.0)))))**2 + p[0]**2
    denum2 = ( x - (2.87*10**9 - (28*10**9 * p[8] * (1/np.sqrt(3)) * (- np.sin (p[9]*np.pi/180.0) * np.cos(p[10]*np.pi/180.0) - np.sin(p[9]*np.pi/180.0) * np.sin(p[10]*np.pi/180.0) + np.cos(p[9]*np.pi/180.0)))))**2 + p[1]**2
    denum3 = ( x - (2.87*10**9 - (28*10**9 * p[8] * (1/np.sqrt(3)) * (np.sin (p[9]*np.pi/180.0) * np.cos(p[10]*np.pi/180.0) - np.sin(p[9]*np.pi/180.0) * np.sin(p[10]*np.pi/180.0) - np.cos(p[9]*np.pi/180.0)))))**2 + p[2]**2
    denum4 = ( x - (2.87*10**9 - (28*10**9 * p[8] * (1/np.sqrt(3)) * (- np.sin (p[9]*np.pi/180.0) * np.cos(p[10]*np.pi/180.0) + np.sin(p[9]*np.pi/180.0) * np.sin(p[10]*np.pi/180.0) - np.cos(p[9]*np.pi/180.0)))))**2 + p[3]**2
    denum5 = ( x - (2.87*10**9 + (28*10**9 * p[8] * (1/np.sqrt(3)) * (- np.sin (p[9]*np.pi/180.0) * np.cos(p[10]*np.pi/180.0) + np.sin(p[9]*np.pi/180.0) * np.sin(p[10]*np.pi/180.0) - np.cos(p[9]*np.pi/180.0)))))**2 + p[3]**2
    denum6 = ( x - (2.87*10**9 + (28*10**9 * p[8] * (1/np.sqrt(3)) * (np.sin (p[9]*np.pi/180.0) * np.cos(p[10]*np.pi/180.0) - np.sin(p[9]*np.pi/180.0) * np.sin(p[10]*np.pi/180.0) - np.cos(p[9]*np.pi/180.0)))))**2 + p[2]**2
    denum7 = ( x - (2.87*10**9 + (28*10**9 * p[8] * (1/np.sqrt(3)) * (- np.sin (p[9]*np.pi/180.0) * np.cos(p[10]*np.pi/180.0) - np.sin(p[9]*np.pi/180.0) * np.sin(p[10]*np.pi/180.0) + np.cos(p[9]*np.pi/180.0)))))**2 + p[1]**2
    denum8 = ( x - (2.87*10**9 + (28*10**9 * p[8] * (1/np.sqrt(3)) * (np.sin (p[9]*np.pi/180.0) * np.cos(p[10]*np.pi/180.0) + np.sin(p[9]*np.pi/180.0) * np.sin(p[10]*np.pi/180.0) + np.cos(p[9]*np.pi/180.0)))))**2 + p[0]**2
    Lor = p[4]*(num/denum)/np.pi
    Lor2 = p[5]*(num2/denum2)/np.pi
    Lor3 = p[6]*(num3/denum3)/np.pi
    Lor4 = p[7]*(num4/denum4)/np.pi
    Lor5 = p[7]*(num5/denum5)/np.pi
    Lor6 = p[6]*(num6/denum6)/np.pi
    Lor7 = p[5]*(num7/denum7)/np.pi
    Lor8 = p[4]*(num8/denum8)/np.pi
    return (1 - (Lor + Lor2 + Lor3 + Lor4 + Lor5 + Lor6 + Lor7 + Lor8))

############################### LEAST-SQUARE FIT ###############################

def residuals(p, y, x):
    return y - Lorentzian(x, p)

pbest = leastsq(residuals, p ,args=(y,x), full_output=1)
best_parameters1 = pbest[0]

fit = Lorentzian(x, best_parameters1)

################# INDIVIDUAL ESR SPECTRA ########################################
num = best_parameters1[0]**2
num2 = best_parameters1[1]**2
num3 = best_parameters1[2]**2
num4 = best_parameters1[3]**2
num5 = best_parameters1[3]**2
num6 = best_parameters1[2]**2
num7 = best_parameters1[1]**2
num8 = best_parameters1[0]**2
denum = ( x - (2.87*10**9 - (28*10**9 * best_parameters1[8] * (1/np.sqrt(3)) * (np.sin (best_parameters1[9]*np.pi/180.0) * np.cos(best_parameters1[10]*np.pi/180.0) + np.sin(best_parameters1[9]*np.pi/180.0) * np.sin(best_parameters1[10]*np.pi/180.0) + np.cos(best_parameters1[9]*np.pi/180.0)))))**2 + best_parameters1[0]**2
denum2 = ( x - (2.87*10**9 - (28*10**9 * best_parameters1[8] * (1/np.sqrt(3)) * (- np.sin (best_parameters1[9]*np.pi/180.0) * np.cos(best_parameters1[10]*np.pi/180.0) - np.sin(best_parameters1[9]*np.pi/180.0) * np.sin(best_parameters1[10]*np.pi/180.0) + np.cos(best_parameters1[9]*np.pi/180.0)))))**2 + best_parameters1[1]**2
denum3 = ( x - (2.87*10**9 - (28*10**9 * best_parameters1[8] * (1/np.sqrt(3)) * (np.sin (best_parameters1[9]*np.pi/180.0) * np.cos(best_parameters1[10]*np.pi/180.0) - np.sin(best_parameters1[9]*np.pi/180.0) * np.sin(best_parameters1[10]*np.pi/180.0) - np.cos(best_parameters1[9]*np.pi/180.0)))))**2 + best_parameters1[2]**2
denum4 = ( x - (2.87*10**9 - (28*10**9 * best_parameters1[8] * (1/np.sqrt(3)) * (- np.sin (best_parameters1[9]*np.pi/180.0) * np.cos(best_parameters1[10]*np.pi/180.0) + np.sin(best_parameters1[9]*np.pi/180.0) * np.sin(best_parameters1[10]*np.pi/180.0) - np.cos(best_parameters1[9]*np.pi/180.0)))))**2 + best_parameters1[3]**2
denum5 = ( x - (2.87*10**9 + (28*10**9 * best_parameters1[8] * (1/np.sqrt(3)) * (- np.sin (best_parameters1[9]*np.pi/180.0) * np.cos(best_parameters1[10]*np.pi/180.0) + np.sin(best_parameters1[9]*np.pi/180.0) * np.sin(best_parameters1[10]*np.pi/180.0) - np.cos(best_parameters1[9]*np.pi/180.0)))))**2 + best_parameters1[3]**2
denum6 = ( x - (2.87*10**9 + (28*10**9 * best_parameters1[8] * (1/np.sqrt(3)) * (np.sin (best_parameters1[9]*np.pi/180.0) * np.cos(best_parameters1[10]*np.pi/180.0) - np.sin(best_parameters1[9]*np.pi/180.0) * np.sin(best_parameters1[10]*np.pi/180.0) - np.cos(best_parameters1[9]*np.pi/180.0)))))**2 + best_parameters1[2]**2
denum7 = ( x - (2.87*10**9 + (28*10**9 * best_parameters1[8] * (1/np.sqrt(3)) * (- np.sin (best_parameters1[9]*np.pi/180.0) * np.cos(best_parameters1[10]*np.pi/180.0) - np.sin(best_parameters1[9]*np.pi/180.0) * np.sin(best_parameters1[10]*np.pi/180.0) + np.cos(best_parameters1[9]*np.pi/180.0)))))**2 + best_parameters1[1]**2
denum8 = ( x - (2.87*10**9 + (28*10**9 * best_parameters1[8] * (1/np.sqrt(3)) * (np.sin (best_parameters1[9]*np.pi/180.0) * np.cos(best_parameters1[10]*np.pi/180.0) + np.sin(best_parameters1[9]*np.pi/180.0) * np.sin(best_parameters1[10]*np.pi/180.0) + np.cos(best_parameters1[9]*np.pi/180.0)))))**2 + best_parameters1[0]**2
Lor = 1- (best_parameters1[4]*(num/denum)/np.pi)
Lor2 = 1- (best_parameters1[5]*(num2/denum2)/np.pi)
Lor3 = 1- (best_parameters1[6]*(num3/denum3)/np.pi)
Lor4 = 1- (best_parameters1[7]*(num4/denum4)/np.pi)
Lor5 = 1- (best_parameters1[7]*(num5/denum5)/np.pi)
Lor6 = 1- (best_parameters1[6]*(num6/denum6)/np.pi)
Lor7 = 1- (best_parameters1[5]*(num7/denum7)/np.pi)
Lor8 = 1- (best_parameters1[4]*(num8/denum8)/np.pi)

#Lorentzian = 1 - (Lor + Lor2 + Lor3 + Lor4 + Lor5 + Lor6 + Lor7 + Lor8)

################ VECTOR COMPONENTS ############################################
X = best_parameters1[4] * np.sin (best_parameters1[5] * np.pi /180.0) * np.cos (best_parameters1[6] * np.pi / 180.0) 
Y = best_parameters1[4] * np.sin (best_parameters1[5] * np.pi /180.0) * np.sin (best_parameters1[6] * np.pi / 180.0)
Z = best_parameters1[4] * np.cos (best_parameters1[5] * np.pi /180.0)
###############################################################################

pl.rc('text', usetex=True)
pl.rc('font', family='serif')
pl.plot (x, y, lw = 1.5, label = 'Experiment')
#pl.plot (x[dipind],y[dipind], 'ro')
pl.plot(x, fit, lw=1.5, color = 'r', label = r"Least-square fit") 
pl.plot (x, Lor, color = 'g', label = r"NV1 $\parallel [111]$")
pl.plot (x, Lor8, color = 'g')
pl.plot (x, Lor2, color = 'y', label = r"NV2 $\parallel [\bar{1} \bar{1} 1]$")
pl.plot (x, Lor7, color = 'y')
pl.plot (x, Lor3, color = 'c', label = r"NV3 $\parallel [1 \bar{1} \bar{1}]$")
pl.plot (x, Lor6, color = 'c')
pl.plot (x, Lor4, color = 'm', label = r"NV4 $\parallel [\bar{1} 1 \bar{1}]$")
pl.plot (x, Lor5, color = 'm')
#pl.plot (x, Lorentzian)
pl.xlabel('Microwave frequency (GHz)', fontsize = 18)
pl.ylabel('Intensity (norm.)', fontsize = 18)
pl.title (r"Position 5")
pl.legend(loc = 4)
pl.show()
#print dipind
#print ind[1]
#print x[ind[1]]
#print p[1]
print ("NV")
print ("B0 =" + str(best_parameters1[8])+r"T")
print ("Theta =" + str(best_parameters1[9]))
print ("Phi =" + str(best_parameters1[10]))
print ("x =" + str(X))
print ("y =" + str(Y))
print ("z =" + str(Z)) 
#print residuals(best_parameters1, y, x)

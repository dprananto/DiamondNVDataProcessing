# -*- coding: utf-8 -*-
"""
Created on Sat Jul 30 14:59:43 2016

@author: dwi
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as pl
from scipy.optimize import leastsq
#from scipy.signal import argrelmin
from scipy import stats
import itertools

Data = pd.read_csv ('/home/dwi/ownCloud/Anlabshared/Dwi/Data and Analysis/160705/2016_07_05_pos1-3.csv')
pl.title (r"Position 1")
Data.columns = ['Microwave frequency', 'Intensity']
x = np.array(Data['Microwave frequency'])
y = np.array(Data['Intensity']/Data['Intensity'].max())
offset = 1.0 - np.mean(y[0:30])
yo = y + offset
posX = 10.0
posY = 0.0
posZ = -11.0
H = 4 * 10 ** 6.0 #HWHM
############################## PEAK SEARCH #####################################
#dipind = argrelmin(y, order=19, mode='clip')
#ind = np.array(dipind[0])

#xdip = x[dipind]
#ydip = y[dipind]
###############################################################################
for c in itertools.permutations (np.linspace(0.0, 360.0, (36*2 + 1)), 2):
############################## INITIAL PARAMETERS ##############################
    p = [0.02, 0.02, 0.02, 0.02, 2*10**(-3), c[0], c[1]]  
# [dip 1-4, B0, Theta, Phi]
 
############################# LORENTZIAN FUNCTION ##############################
    def Lorentzian(x, p):
        num = H**2
        num2 = H**2
        num3= H**2
        num4 = H**2
        num5 = H**2
        num6 = H**2
        num7 = H**2
        num8 = H**2
        denum = ( x - (2.87*10**9 - (28*10**9 * p[4] * (1/np.sqrt(3)) * (np.sin (p[5]*np.pi/180.0) * np.cos(p[6]*np.pi/180.0) + np.sin(p[5]*np.pi/180.0) * np.sin(p[6]*np.pi/180.0) + np.cos(p[5]*np.pi/180.0)))))**2 + H**2
        denum2 = ( x - (2.87*10**9 - (28*10**9 * p[4] * (1/np.sqrt(3)) * (- np.sin (p[5]*np.pi/180.0) * np.cos(p[6]*np.pi/180.0) - np.sin(p[5]*np.pi/180.0) * np.sin(p[6]*np.pi/180.0) + np.cos(p[5]*np.pi/180.0)))))**2 + H**2
        denum3 = ( x - (2.87*10**9 - (28*10**9 * p[4] * (1/np.sqrt(3)) * (np.sin (p[5]*np.pi/180.0) * np.cos(p[6]*np.pi/180.0) - np.sin(p[5]*np.pi/180.0) * np.sin(p[6]*np.pi/180.0) - np.cos(p[5]*np.pi/180.0)))))**2 + H**2
        denum4 = ( x - (2.87*10**9 - (28*10**9 * p[4] * (1/np.sqrt(3)) * (- np.sin (p[5]*np.pi/180.0) * np.cos(p[6]*np.pi/180.0) + np.sin(p[5]*np.pi/180.0) * np.sin(p[6]*np.pi/180.0) - np.cos(p[5]*np.pi/180.0)))))**2 + H**2
        denum5 = ( x - (2.87*10**9 + (28*10**9 * p[4] * (1/np.sqrt(3)) * (- np.sin (p[5]*np.pi/180.0) * np.cos(p[6]*np.pi/180.0) + np.sin(p[5]*np.pi/180.0) * np.sin(p[6]*np.pi/180.0) - np.cos(p[5]*np.pi/180.0)))))**2 + H**2
        denum6 = ( x - (2.87*10**9 + (28*10**9 * p[4] * (1/np.sqrt(3)) * (np.sin (p[5]*np.pi/180.0) * np.cos(p[6]*np.pi/180.0) - np.sin(p[5]*np.pi/180.0) * np.sin(p[6]*np.pi/180.0) - np.cos(p[5]*np.pi/180.0)))))**2 + H**2
        denum7 = ( x - (2.87*10**9 + (28*10**9 * p[4] * (1/np.sqrt(3)) * (- np.sin (p[5]*np.pi/180.0) * np.cos(p[6]*np.pi/180.0) - np.sin(p[5]*np.pi/180.0) * np.sin(p[6]*np.pi/180.0) + np.cos(p[5]*np.pi/180.0)))))**2 + H**2
        denum8 = ( x - (2.87*10**9 + (28*10**9 * p[4] * (1/np.sqrt(3)) * (np.sin (p[5]*np.pi/180.0) * np.cos(p[6]*np.pi/180.0) + np.sin(p[5]*np.pi/180.0) * np.sin(p[6]*np.pi/180.0) + np.cos(p[5]*np.pi/180.0)))))**2 + H**2
        Lor = p[0]*(num/denum)/np.pi
        Lor2 = p[1]*(num2/denum2)/np.pi
        Lor3 = p[2]*(num3/denum3)/np.pi
        Lor4 = p[3]*(num4/denum4)/np.pi
        Lor5 = p[3]*(num5/denum5)/np.pi
        Lor6 = p[2]*(num6/denum6)/np.pi
        Lor7 = p[1]*(num7/denum7)/np.pi
        Lor8 = p[0]*(num8/denum8)/np.pi
        return (1 - (Lor + Lor2 + Lor3 + Lor4 + Lor5 + Lor6 + Lor7 + Lor8))

############################### LEAST-SQUARE FIT ###############################

    def residuals(p, yo, x):
        return yo - Lorentzian(x, p)

    pbest = leastsq(residuals, p ,args=(yo, x), full_output=1)
    best_parameters1 = pbest[0]

    fit = Lorentzian(x, best_parameters1)

################# INDIVIDUAL ESR SPECTRA #######################################
    #num = H**2
    #num2 = H**2
    #num3 = H**2
    #num4 = H**2
    #num5 = H**2
    #num6 = H**2
    #num7 = H**2
    #num8 = H**2
    #denum = ( x - (2.87*10**9 - (28*10**9 * best_parameters1[4] * (1/np.sqrt(3)) * (np.sin (best_parameters1[5]*np.pi/180.0) * np.cos(best_parameters1[6]*np.pi/180.0) + np.sin(best_parameters1[5]*np.pi/180.0) * np.sin(best_parameters1[6]*np.pi/180.0) + np.cos(best_parameters1[5]*np.pi/180.0)))))**2 + H**2
    #denum2 = ( x - (2.87*10**9 - (28*10**9 * best_parameters1[4] * (1/np.sqrt(3)) * (- np.sin (best_parameters1[5]*np.pi/180.0) * np.cos(best_parameters1[6]*np.pi/180.0) - np.sin(best_parameters1[5]*np.pi/180.0) * np.sin(best_parameters1[6]*np.pi/180.0) + np.cos(best_parameters1[5]*np.pi/180.0)))))**2 + H**2
    #denum3 = ( x - (2.87*10**9 - (28*10**9 * best_parameters1[4] * (1/np.sqrt(3)) * (np.sin (best_parameters1[5]*np.pi/180.0) * np.cos(best_parameters1[6]*np.pi/180.0) - np.sin(best_parameters1[5]*np.pi/180.0) * np.sin(best_parameters1[6]*np.pi/180.0) - np.cos(best_parameters1[5]*np.pi/180.0)))))**2 + H**2
    #denum4 = ( x - (2.87*10**9 - (28*10**9 * best_parameters1[4] * (1/np.sqrt(3)) * (- np.sin (best_parameters1[5]*np.pi/180.0) * np.cos(best_parameters1[6]*np.pi/180.0) + np.sin(best_parameters1[5]*np.pi/180.0) * np.sin(best_parameters1[6]*np.pi/180.0) - np.cos(best_parameters1[5]*np.pi/180.0)))))**2 + H**2
    #denum5 = ( x - (2.87*10**9 + (28*10**9 * best_parameters1[4] * (1/np.sqrt(3)) * (- np.sin (best_parameters1[5]*np.pi/180.0) * np.cos(best_parameters1[6]*np.pi/180.0) + np.sin(best_parameters1[5]*np.pi/180.0) * np.sin(best_parameters1[6]*np.pi/180.0) - np.cos(best_parameters1[5]*np.pi/180.0)))))**2 + H**2
    #denum6 = ( x - (2.87*10**9 + (28*10**9 * best_parameters1[4] * (1/np.sqrt(3)) * (np.sin (best_parameters1[5]*np.pi/180.0) * np.cos(best_parameters1[6]*np.pi/180.0) - np.sin(best_parameters1[5]*np.pi/180.0) * np.sin(best_parameters1[6]*np.pi/180.0) - np.cos(best_parameters1[5]*np.pi/180.0)))))**2 + H**2
    #denum7 = ( x - (2.87*10**9 + (28*10**9 * best_parameters1[4] * (1/np.sqrt(3)) * (- np.sin (best_parameters1[5]*np.pi/180.0) * np.cos(best_parameters1[6]*np.pi/180.0) - np.sin(best_parameters1[5]*np.pi/180.0) * np.sin(best_parameters1[6]*np.pi/180.0) + np.cos(best_parameters1[5]*np.pi/180.0)))))**2 + H**2
    #denum8 = ( x - (2.87*10**9 + (28*10**9 * best_parameters1[4] * (1/np.sqrt(3)) * (np.sin (best_parameters1[5]*np.pi/180.0) * np.cos(best_parameters1[6]*np.pi/180.0) + np.sin(best_parameters1[5]*np.pi/180.0) * np.sin(best_parameters1[6]*np.pi/180.0) + np.cos(best_parameters1[5]*np.pi/180.0)))))**2 + H**2
    #Lor = 1- (best_parameters1[0]*(num/denum)/np.pi)
    #Lor2 = 1- (best_parameters1[1]*(num2/denum2)/np.pi)
    #Lor3 = 1- (best_parameters1[2]*(num3/denum3)/np.pi)
    #Lor4 = 1- (best_parameters1[3]*(num4/denum4)/np.pi)
    #Lor5 = 1- (best_parameters1[3]*(num5/denum5)/np.pi)
    #Lor6 = 1- (best_parameters1[2]*(num6/denum6)/np.pi)
    #Lor7 = 1- (best_parameters1[1]*(num7/denum7)/np.pi)
    #Lor8 = 1- (best_parameters1[0]*(num8/denum8)/np.pi)

    #Lorentzian = 1 - (Lor + Lor2 + Lor3 + Lor4 + Lor5 + Lor6 + Lor7 + Lor8)
########################## VECTOR COMPONENTS ############################################
    Bx = best_parameters1[4] * np.sin (best_parameters1[5] * np.pi /180.0) * np.cos (best_parameters1[6] * np.pi / 180.0) 
    By = best_parameters1[4] * np.sin (best_parameters1[5] * np.pi /180.0) * np.sin (best_parameters1[6] * np.pi / 180.0)
    Bz = best_parameters1[4] * np.cos (best_parameters1[5] * np.pi /180.0)
###############################################################################
#print dipind
#print ind[1]
#print x[ind[1]]
#print p[1]
########################## KOLMOGOROV-SMIRNOV TEST ############################
    KS = stats.chisquare (yo, f_exp = fit)
####################### WRITE DATAFRAME TO CSV ################################

    data_frame = {'PosX': [posX], 'PosY': [posY], 'PosZ': [posZ], 'B': [best_parameters1[4]], 'T': [best_parameters1[5]], 'P': [best_parameters1[6]], 'Bx': [Bx], 'By': [By], 'Bz': [Bz], 'Bi': [p[4]], 'Ti': [p[5]], 'Pi': [p[6]], 'Chi': [KS[0]], 'p-value': [KS[1]], 'HWHM': [H]}
    df = pd.DataFrame(data_frame, columns = ['PosX', 'PosY', 'PosZ', 'B', 'T', 'P', 'Bx', 'By', 'Bz', 'Bi', 'Ti', 'Pi', 'Chi', 'p-value', 'HWHM'])
    df.to_csv ('/home/dwi/ownCloud/Anlabshared/Dwi/Data and Analysis/160705/analysis_position1-1chi.csv', mode = 'a', header = False)

###############################################################################

print ("Local Magnetic Field Vector at Certain Position")
#print ("B =" + str(best_parameters1[4])+r"T")
#print ("Theta =" + str(best_parameters1[5]))
#print ("Phi =" + str(best_parameters1[6]))
#print ("x =" + str(X))
#print ("y =" + str(Y))
#print ("z =" + str(Z)) 
#print (df)
#print (KS)
#print residuals(best_parameters1, y, x)
############################### PLOT ##########################################
#pl.rc('text', usetex=True)
#pl.rc('font', family='serif')
#pl.plot (x, y, lw = 1.5, label = 'Experiment')
#pl.plot (x[dipind],y[dipind], 'ro')
#pl.plot(x, fit, lw=1.5, color = 'r', label = r"Least-square fit") 
#pl.plot (x, Lor, color = 'g', label = r"NV1 $\parallel [111]$")
#pl.plot (x, Lor8, color = 'g')
#pl.plot (x, Lor2, color = 'y', label = r"NV2 $\parallel [\bar{1} \bar{1} 1]$")
#pl.plot (x, Lor7, color = 'y')
#pl.plot (x, Lor3, color = 'c', label = r"NV3 $\parallel [1 \bar{1} \bar{1}]$")
#pl.plot (x, Lor6, color = 'c')
#pl.plot (x, Lor4, color = 'm', label = r"NV4 $\parallel [\bar{1} 1 \bar{1}]$")
#pl.plot (x, Lor5, color = 'm')
#pl.plot (x, Lorentzian)
#pl.xlabel('Microwave frequency (GHz)', fontsize = 18)
#pl.ylabel('Intensity (norm.)', fontsize = 18)

#pl.legend(loc = 4)
#pl.show()
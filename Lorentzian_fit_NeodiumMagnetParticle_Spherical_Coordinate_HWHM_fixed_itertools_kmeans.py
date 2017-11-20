###############################################################################
#THIS PROGRAM WILL GENERATE ALL POSSIBLE COMBINATION OF THETA AND PHI AND 
#PERFORM FITTINGS FOR EACH COMBINATION
###############################################################################
import pandas as pd
import numpy as np
#import matplotlib.pyplot as pl
from scipy.optimize import leastsq
#from scipy.signal import argrelmin
from scipy import stats
import itertools

from sklearn.cluster import KMeans

fi = input ('File path to open: ')
fo = input ('File path to save: ')
fc = input ('File path to save clustered data: ')
Data = pd.read_csv (fi)
#pl.title (r"Position 1")
Data.columns = ['Microwave frequency', 'Intensity']
x = np.array(Data['Microwave frequency'])
y = np.array(Data['Intensity']/Data['Intensity'].max())
offset = 1.0 - np.mean(y[0:30])
yo = y + offset
posX = input ('Coordinate x:')
posY = input ('Coordinate y:')
posZ = input ('Coordinate z:')

h = input ('HWHM (*10^6) = ')
H = float (h) * 10 ** 6.0 #Half Width Half Maximum
p4 = input ('B (*10^-3)= ')
P4 = float (p4) *10**(-3)

############################## PEAK SEARCH #####################################
#dipind = argrelmin(y, order=19, mode='clip')
#ind = np.array(dipind[0])

#xdip = x[dipind]
#ydip = y[dipind]
###############################################################################
for c in itertools.permutations (np.linspace(0.0, 360.0, 37), 2):
############################## INITIAL PARAMETERS ##############################
    p = [0.02, 0.02, 0.02, 0.02, P4, c[0], c[1]]  # [dip 1-4, B0, Theta, Phi]
 
############################# LORENTZIAN FUNCTION ##############################
    def Lorentzian(x, p):
        num = H**2
        num2 = H**2
        num3 = H**2
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
    B = np.sqrt (Bx**2 + By**2 + Bz**2)
################################# ANGLE BETWEEN NV AND B #######################################
    BS1 = (best_parameters1[4] * (1/np.sqrt(3)) * (np.sin (best_parameters1[5]*np.pi/180.0) * np.cos(best_parameters1[6]*np.pi/180.0) + np.sin(best_parameters1[5]*np.pi/180.0) * np.sin(best_parameters1[6]*np.pi/180.0) + np.cos(best_parameters1[5]*np.pi/180.0)))
    BS2 = (best_parameters1[4] * (1/np.sqrt(3)) * (- np.sin (best_parameters1[5]*np.pi/180.0) * np.cos(best_parameters1[6]*np.pi/180.0) - np.sin(best_parameters1[5]*np.pi/180.0) * np.sin(best_parameters1[6]*np.pi/180.0) + np.cos(best_parameters1[5]*np.pi/180.0)))
    BS3 = (best_parameters1[4] * (1/np.sqrt(3)) * (np.sin (best_parameters1[5]*np.pi/180.0) * np.cos(best_parameters1[6]*np.pi/180.0) - np.sin(best_parameters1[5]*np.pi/180.0) * np.sin(best_parameters1[6]*np.pi/180.0) - np.cos(best_parameters1[5]*np.pi/180.0)))
    BS4 = (best_parameters1[4] * (1/np.sqrt(3)) * (- np.sin (best_parameters1[5]*np.pi/180.0) * np.cos(best_parameters1[6]*np.pi/180.0) + np.sin(best_parameters1[5]*np.pi/180.0) * np.sin(best_parameters1[6]*np.pi/180.0) - np.cos(best_parameters1[5]*np.pi/180.0)))

    t1 = np.arccos (BS1/(B*1)) * 180.0 / np.pi
    t2 = np.arccos (BS2/(B*1)) * 180.0 / np.pi
    t3 = np.arccos (BS3/(B*1)) * 180.0 / np.pi
    t4 = np.arccos (BS4/(B*1)) * 180.0 / np.pi
#print dipind
#print ind[1]
#print x[ind[1]]
#print p[1]
########################## KOLMOGOROV-SMIRNOV TEST ############################
    KS = stats.ks_2samp (yo, fit)
####################### WRITE DATAFRAME TO CSV ################################
    #if (best_parameters1[5] <= 360.0 and best_parameters1[5] >= 0.0 and best_parameters1[6] <= 360.0 and best_parameters1[6] >= 0.0):
    data_frame = {'PosX': [posX], 'PosY': [posY], 'PosZ': [posZ], 'B': [best_parameters1[4]], 'T': [best_parameters1[5]], 'P': [best_parameters1[6]], 'Bx': [Bx], 'By': [By], 'Bz': [Bz], 'Bi': [p[4]], 'Ti': [p[5]], 'Pi': [p[6]], 'D': [KS[0]], 'p-value': [KS[1]], 'HWHM': [H], 'theta1': [t1], 'theta2': [t2], 'theta3': [t3], 'theta4': [t4]}
    df = pd.DataFrame(data_frame, columns = ['PosX', 'PosY', 'PosZ', 'B', 'T', 'P', 'Bx', 'By', 'Bz', 'Bi', 'Ti', 'Pi', 'D', 'p-value', 'HWHM', 'theta1', 'theta2', 'theta3', 'theta4'])
    df.to_csv (fo, mode = 'a', header = False)

######################### K-means CLUSTERING############################################

ro = pd.read_csv(fo) #Read to be clustered data
ro.columns = ['Index', 'PosX', 'PosY', 'PosZ', 'B', 'T', 'P', 'Bx', 'By', 'Bz', 'Bi', 'Ti', 'Pi', 'D', 'p-value', 'HWHM', 'theta1', 'theta2', 'theta3', 'theta4']

ros = ro.iloc[:, 16:20] #Take theta 1-4 Data
rov = ros.values #Dataframe to array
Cluster = KMeans (4).fit(rov) #KMeans Clustering

CL = Cluster.labels_

df_processed = ro.copy()
df_processed['Cluster Class'] = pd.Series(CL, index=df_processed.index)
df_processed.to_csv(fc, header = True)

#print (ros)

#print ("Local Magnetic Field Vector at Certain Position")
#print (CL)
###############################################################################
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
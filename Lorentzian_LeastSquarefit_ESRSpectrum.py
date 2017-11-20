###############################################################################
#THIS PROGRAM PLOT ESR SPECTRUM WITH FITTING###################################
#THIS PROGRAM CAN BE USED TO PLOT ESR SPECTRUM WITH DTA FROM ITERATIVE FITTING#
###############################################################################
import pandas as pd
import numpy as np
import matplotlib.pyplot as pl
from scipy.optimize import leastsq
#from scipy.signal import argrelmin
from scipy import stats

#fi = input ('File path to open:')
path = r'/home/dwi/ownCloud/Anlabshared/Dwi/researchdata/160830 Magnetic particle stray mapping/scan 25 x 25 um/'
Data = pd.read_csv (r'/home/dwi/ownCloud/Anlabshared/Dwi/researchdata/160830 Magnetic particle stray mapping/scan 25 x 25 um/2016_08_30_S1_P93.csv')
#ttl = input (r"Graph title: " )
#pl.title (ttl)
Data.columns = ['Microwave frequency', 'Intensity']
x = np.array(Data['Microwave frequency'])
y = np.array(Data['Intensity']/Data['Intensity'].max())
offset = 1.0 - np.mean (y[0:20])
yo = y + offset
posX = 10.0
posY = 0.0
posZ = -11.0

h = float(input ('HWHM (*10^6) = '))
H = h * 10 ** 6.0 #Half Width Half Maximum
p4 = float(input ('B (*10^-3)= '))
P4 = p4 *10**(-3) #Magnetic Field
p5 = input ('Theta (deg.) = ')
P5 = float (p5)
p6 = input ('Phi (deg.) = ')
P6 = float (p6)
############################## PEAK SEARCH #####################################
#dipind = argrelmin(y, order=19, mode='clip')
#ind = np.array(dipind[0])

#xdip = x[dipind]
#ydip = y[dipind]

############################## INITIAL PARAMETERS ##############################
p = [0.02, 0.02, 0.02, 0.02, P4, P5, P6, H]  
# [dip 1-4, B0, Theta, Phi]
############################# LORENTZIAN FUNCTION ##############################
def Lorentzian(x, p):
    num = p[7]**2
    num2 = p[7]**2
    num3 = p[7]**2
    num4 = p[7]**2
    num5 = p[7]**2
    num6 = p[7]**2
    num7 = p[7]**2
    num8 = p[7]**2
    denum = ( x - (2.87*10**9 - (28*10**9 * p[4] * (1/np.sqrt(3)) * (np.sin (p[5]*np.pi/180.0) * np.cos(p[6]*np.pi/180.0) + np.sin(p[5]*np.pi/180.0) * np.sin(p[6]*np.pi/180.0) + np.cos(p[5]*np.pi/180.0)))))**2 + p[7]**2
    denum2 = ( x - (2.87*10**9 - (28*10**9 * p[4] * (1/np.sqrt(3)) * (- np.sin (p[5]*np.pi/180.0) * np.cos(p[6]*np.pi/180.0) - np.sin(p[5]*np.pi/180.0) * np.sin(p[6]*np.pi/180.0) + np.cos(p[5]*np.pi/180.0)))))**2 + p[7]**2
    denum3 = ( x - (2.87*10**9 - (28*10**9 * p[4] * (1/np.sqrt(3)) * (np.sin (p[5]*np.pi/180.0) * np.cos(p[6]*np.pi/180.0) - np.sin(p[5]*np.pi/180.0) * np.sin(p[6]*np.pi/180.0) - np.cos(p[5]*np.pi/180.0)))))**2 + p[7]**2
    denum4 = ( x - (2.87*10**9 - (28*10**9 * p[4] * (1/np.sqrt(3)) * (- np.sin (p[5]*np.pi/180.0) * np.cos(p[6]*np.pi/180.0) + np.sin(p[5]*np.pi/180.0) * np.sin(p[6]*np.pi/180.0) - np.cos(p[5]*np.pi/180.0)))))**2 + p[7]**2
    denum5 = ( x - (2.87*10**9 + (28*10**9 * p[4] * (1/np.sqrt(3)) * (- np.sin (p[5]*np.pi/180.0) * np.cos(p[6]*np.pi/180.0) + np.sin(p[5]*np.pi/180.0) * np.sin(p[6]*np.pi/180.0) - np.cos(p[5]*np.pi/180.0)))))**2 + p[7]**2
    denum6 = ( x - (2.87*10**9 + (28*10**9 * p[4] * (1/np.sqrt(3)) * (np.sin (p[5]*np.pi/180.0) * np.cos(p[6]*np.pi/180.0) - np.sin(p[5]*np.pi/180.0) * np.sin(p[6]*np.pi/180.0) - np.cos(p[5]*np.pi/180.0)))))**2 + p[7]**2
    denum7 = ( x - (2.87*10**9 + (28*10**9 * p[4] * (1/np.sqrt(3)) * (- np.sin (p[5]*np.pi/180.0) * np.cos(p[6]*np.pi/180.0) - np.sin(p[5]*np.pi/180.0) * np.sin(p[6]*np.pi/180.0) + np.cos(p[5]*np.pi/180.0)))))**2 + p[7]**2
    denum8 = ( x - (2.87*10**9 + (28*10**9 * p[4] * (1/np.sqrt(3)) * (np.sin (p[5]*np.pi/180.0) * np.cos(p[6]*np.pi/180.0) + np.sin(p[5]*np.pi/180.0) * np.sin(p[6]*np.pi/180.0) + np.cos(p[5]*np.pi/180.0)))))**2 + p[7]**2
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

pbest = leastsq(residuals, p ,args=(yo,x), full_output=1)
best_parameters1 = pbest[0]

fit = Lorentzian(x, best_parameters1)

################# INDIVIDUAL ESR SPECTRA #######################################
num = best_parameters1[7]**2
num2 = best_parameters1[7]**2
num3 = best_parameters1[7]**2
num4 = best_parameters1[7]**2
num5 = best_parameters1[7]**2
num6 = best_parameters1[7]**2
num7 = best_parameters1[7]**2
num8 = best_parameters1[7]**2
denum = ( x - (2.87*10**9 - (28*10**9 * best_parameters1[4] * (1/np.sqrt(3)) * (np.sin (best_parameters1[5]*np.pi/180.0) * np.cos(best_parameters1[6]*np.pi/180.0) + np.sin(best_parameters1[5]*np.pi/180.0) * np.sin(best_parameters1[6]*np.pi/180.0) + np.cos(best_parameters1[5]*np.pi/180.0)))))**2 + best_parameters1[7]**2
denum2 = ( x - (2.87*10**9 - (28*10**9 * best_parameters1[4] * (1/np.sqrt(3)) * (- np.sin (best_parameters1[5]*np.pi/180.0) * np.cos(best_parameters1[6]*np.pi/180.0) - np.sin(best_parameters1[5]*np.pi/180.0) * np.sin(best_parameters1[6]*np.pi/180.0) + np.cos(best_parameters1[5]*np.pi/180.0)))))**2 + best_parameters1[7]**2
denum3 = ( x - (2.87*10**9 - (28*10**9 * best_parameters1[4] * (1/np.sqrt(3)) * (np.sin (best_parameters1[5]*np.pi/180.0) * np.cos(best_parameters1[6]*np.pi/180.0) - np.sin(best_parameters1[5]*np.pi/180.0) * np.sin(best_parameters1[6]*np.pi/180.0) - np.cos(best_parameters1[5]*np.pi/180.0)))))**2 + best_parameters1[7]**2
denum4 = ( x - (2.87*10**9 - (28*10**9 * best_parameters1[4] * (1/np.sqrt(3)) * (- np.sin (best_parameters1[5]*np.pi/180.0) * np.cos(best_parameters1[6]*np.pi/180.0) + np.sin(best_parameters1[5]*np.pi/180.0) * np.sin(best_parameters1[6]*np.pi/180.0) - np.cos(best_parameters1[5]*np.pi/180.0)))))**2 + best_parameters1[7]**2
denum5 = ( x - (2.87*10**9 + (28*10**9 * best_parameters1[4] * (1/np.sqrt(3)) * (- np.sin (best_parameters1[5]*np.pi/180.0) * np.cos(best_parameters1[6]*np.pi/180.0) + np.sin(best_parameters1[5]*np.pi/180.0) * np.sin(best_parameters1[6]*np.pi/180.0) - np.cos(best_parameters1[5]*np.pi/180.0)))))**2 + best_parameters1[7]**2
denum6 = ( x - (2.87*10**9 + (28*10**9 * best_parameters1[4] * (1/np.sqrt(3)) * (np.sin (best_parameters1[5]*np.pi/180.0) * np.cos(best_parameters1[6]*np.pi/180.0) - np.sin(best_parameters1[5]*np.pi/180.0) * np.sin(best_parameters1[6]*np.pi/180.0) - np.cos(best_parameters1[5]*np.pi/180.0)))))**2 + best_parameters1[7]**2
denum7 = ( x - (2.87*10**9 + (28*10**9 * best_parameters1[4] * (1/np.sqrt(3)) * (- np.sin (best_parameters1[5]*np.pi/180.0) * np.cos(best_parameters1[6]*np.pi/180.0) - np.sin(best_parameters1[5]*np.pi/180.0) * np.sin(best_parameters1[6]*np.pi/180.0) + np.cos(best_parameters1[5]*np.pi/180.0)))))**2 + best_parameters1[7]**2
denum8 = ( x - (2.87*10**9 + (28*10**9 * best_parameters1[4] * (1/np.sqrt(3)) * (np.sin (best_parameters1[5]*np.pi/180.0) * np.cos(best_parameters1[6]*np.pi/180.0) + np.sin(best_parameters1[5]*np.pi/180.0) * np.sin(best_parameters1[6]*np.pi/180.0) + np.cos(best_parameters1[5]*np.pi/180.0)))))**2 + best_parameters1[7]**2
Lor = 1- (best_parameters1[0]*(num/denum)/np.pi)
Lor2 = 1- (best_parameters1[1]*(num2/denum2)/np.pi)
Lor3 = 1- (best_parameters1[2]*(num3/denum3)/np.pi)
Lor4 = 1- (best_parameters1[3]*(num4/denum4)/np.pi)
Lor5 = 1- (best_parameters1[3]*(num5/denum5)/np.pi)
Lor6 = 1- (best_parameters1[2]*(num6/denum6)/np.pi)
Lor7 = 1- (best_parameters1[1]*(num7/denum7)/np.pi)
Lor8 = 1- (best_parameters1[0]*(num8/denum8)/np.pi)

#Lorentzian = 1 - (Lor + Lor2 + Lor3 + Lor4 + Lor5 + Lor6 + Lor7 + Lor8)
###############################################################################
pl.rc('text', usetex=True)
pl.rc('font', family='serif', size = 24)
pl.plot (x, yo, lw = 1.5, color = 'r', label = 'Experiment')
#pl.plot (x[dipind],y[dipind], 'ro')
pl.plot(x, fit, lw=2.5, color = 'b', label = r"Least-square fit") 
pl.plot (x, Lor, lw = 2, color = 'g', label = r"NV1") #$\parallel [111]$")
pl.plot (x, Lor8, lw = 2,  color = 'g')
pl.plot (x, Lor2, lw = 2, color = 'y', label = r"NV2") #$\parallel [\bar{1} \bar{1} 1]$")
pl.plot (x, Lor7, lw = 2, color = 'y')
pl.plot (x, Lor3, lw = 2, color = 'c', label = r"NV3") #$\parallel [1 \bar{1} \bar{1}]$")
pl.plot (x, Lor6, lw = 2, color = 'c')
pl.plot (x, Lor4, lw = 2, color = 'm', label = r"NV4") #$\parallel [\bar{1} 1 \bar{1}]$")
pl.plot (x, Lor5, lw = 2, color = 'm')
#pl.plot (x, Lorentzian)
pl.xlabel(r"$f$ (GHz)", fontsize = 24)
pl.ylabel('Intensity (norm.)', fontsize = 24)

pl.legend(loc=0, ncol=3, fontsize = 16)
pl.savefig(path + 'ESR_P93.svg', format = 'svg', bbox_inches = 'tight', frameon = False)
pl.show()


################ VECTOR COMPONENTS ############################################
Bx = best_parameters1[4] * np.sin (best_parameters1[5] * np.pi /180.0) * np.cos (best_parameters1[6] * np.pi / 180.0) 
By = best_parameters1[4] * np.sin (best_parameters1[5] * np.pi /180.0) * np.sin (best_parameters1[6] * np.pi / 180.0)
Bz = best_parameters1[4] * np.cos (best_parameters1[5] * np.pi /180.0)

B = np.sqrt (Bx**2 + By**2 + Bz**2)
################################# ANGLE BETWEEN NV AND B #######################################
BS1 = (best_parameters1[4] * (1/np.sqrt(3)) * (np.sin (best_parameters1[5]*np.pi/180.0) * np.cos(best_parameters1[6]*np.pi/180.0) + np.sin(best_parameters1[5]*np.pi/180.0) * np.sin(best_parameters1[6]*np.pi/180.0) + np.cos(best_parameters1[5]*np.pi/180.0)))
BS2 = (best_parameters1[4] * (1/np.sqrt(3)) * (- np.sin (best_parameters1[5]*np.pi/180.0) * np.cos(best_parameters1[6]*np.pi/180.0) - np.sin(best_parameters1[5]*np.pi/180.0) * np.sin(best_parameters1[6]*np.pi/180.0) + np.cos(best_parameters1[5]*np.pi/180.0)))
BS3 = (best_parameters1[4] * (1/np.sqrt(3)) * (np.sin (best_parameters1[5]*np.pi/180.0) * np.cos(best_parameters1[6]*np.pi/180.0) - np.sin(best_parameters1[5]*np.pi/180.0) * np.sin(best_parameters1[6]*np.pi/180.0) - np.cos(best_parameters1[5]*np.pi/180.0)))
BS4 = (best_parameters1[4] * (1/np.sqrt(3)) * (- np.sin (best_parameters1[5]*np.pi/180.0) * np.cos(best_parameters1[6]*np.pi/180.0) + np.sin(best_parameters1[5]*np.pi/180.0) * np.sin(best_parameters1[6]*np.pi/180.0) - np.cos(best_parameters1[5]*np.pi/180.0)))

t1 = np.arccos (BS1/(B*1)) * 180 / np.pi
t2 = np.arccos (BS2/(B*1)) * 180 / np.pi
t3 = np.arccos (BS3/(B*1)) * 180 / np.pi
t4 = np.arccos (BS4/(B*1)) * 180 / np.pi

###############################################################################
#print dipind
#print ind[1]
#print x[ind[1]]
#print p[1]
########################## KOLMOGOROV-SMIRNOV TEST ############################
KS = stats.ks_2samp (yo, fit)
####################### WRITE DATAFRAME TO CSV ################################

data_frame = {'PosX': [posX], 'PosY': [posY], 'PosZ': [posZ], 'B': [best_parameters1[4]], 'T': [best_parameters1[5]], 'P': [best_parameters1[6]], 'Bx': [Bx], 'By': [By], 'Bz': [Bz], 'Bi': [p[4]], 'Ti': [p[5]], 'Pi': [p[6]], 'D': [KS[0]], 'p-value': [KS[1]]}
df = pd.DataFrame(data_frame, columns = ['PosX', 'PosY', 'PosZ', 'B', 'T', 'P', 'Bx', 'By', 'Bz', 'Bi', 'Ti', 'Pi', 'D', 'p-value'])
#df.to_csv ('/home/dwi/ownCloud/Anlabshared/Dwi/Data and Analysis/160705/analysis2.csv', mode = 'a', header = False)

###############################################################################

print ("Local Magnetic Field Vector at Certain Position")
print ("B =" + str(best_parameters1[4])+r"T")
print ("Theta =" + str(best_parameters1[5]))
print ("Phi =" + str(best_parameters1[6]))
print ("HWHM = " + str(best_parameters1[7]))
#print ("x =" + str(X))
#print ("y =" + str(Y))
#print ("z =" + str(Z)) 
#print (df)
#print (KS)
#print residuals(best_parameters1, y, x)
print ("theta1 = " + str(t1))
print ("theta2 = " + str(t2))
print ("theta3 = " + str(t3))
print ("theta4 = " + str(t4))
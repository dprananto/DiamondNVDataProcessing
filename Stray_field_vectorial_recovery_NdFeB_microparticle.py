###############################################################################
#THIS PROGRAM WILL GENERATE ALL POSSIBLE COMBINATION OF THETA AND PHI AND 
#PERFORM FITTINGS FOR EACH COMBINATION
###############################################################################
import pandas as pd
import numpy as np
from scipy.optimize import leastsq
from scipy import stats
import itertools

from sklearn.cluster import KMeans

fi = input ('File path to open: ')
fo = input ('File path to save: ')
fc = input ('File path to save clustered data: ')
Data = pd.read_csv (fi)
Data.columns = ['Microwave frequency', 'Intensity']
x = np.array(Data['Microwave frequency'])
y = np.array(Data['Intensity']/Data['Intensity'].max())
offset = 1.0 - np.mean(y[0:30])
yo = y + offset
posX = float(input ('Coordinate x:'))
posY = float(input ('Coordinate y:'))
posZ = float(input ('Coordinate z:'))
regionT = int(input ('Select region of vertical angle:'))
regionP = int(input ('Select region of horizontal angle:'))
h = float(input ('HWHM (*10^6) = '))
H = h * 10 ** 6.0 #Half Width Half Maximum
p4 = float (input ('B (*10^-3)= ')) #Magnetic Field
P4 = p4 * 10 ** -3.0
  
deg1 = 0.0 
deg2 = 90.0
deg3 = 90.0
deg4 = 180.0
deg5 = 180.0
deg6 = 270.0
deg7 = 270.0
deg8 = 360.0  
  
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

########################## KOLMOGOROV-SMIRNOV TEST ############################
    KS = stats.ks_2samp (yo, fit)
####################### WRITE DATAFRAME TO CSV ################################
    if (regionT == 1 and regionP == 1):
        if (best_parameters1[5] <= deg2 and best_parameters1[5] >= deg1 and best_parameters1[6] <= deg2 and best_parameters1[6] >= deg1):
            data_frame = {'PosX': [posX], 'PosY': [posY], 'PosZ': [posZ], 'B': [best_parameters1[4]], 'T': [best_parameters1[5]], 'P': [best_parameters1[6]], 'Bx': [Bx], 'By': [By], 'Bz': [Bz], 'Bi': [p[4]], 'Ti': [p[5]], 'Pi': [p[6]], 'D': [KS[0]], 'p-value': [KS[1]], 'HWHM': [H], 'theta1': [t1], 'theta2': [t2], 'theta3': [t3], 'theta4': [t4]}
            df = pd.DataFrame(data_frame, columns = ['PosX', 'PosY', 'PosZ', 'B', 'T', 'P', 'Bx', 'By', 'Bz', 'Bi', 'Ti', 'Pi', 'D', 'p-value', 'HWHM', 'theta1', 'theta2', 'theta3', 'theta4'])
            df.to_csv (fo, mode = 'a', header = False)
    elif (regionT == 1 and regionP == 2):
        if (best_parameters1[5] <= deg2 and best_parameters1[5] >= deg1 and best_parameters1[6] <= deg4 and best_parameters1[6] >= deg3):
            data_frame = {'PosX': [posX], 'PosY': [posY], 'PosZ': [posZ], 'B': [best_parameters1[4]], 'T': [best_parameters1[5]], 'P': [best_parameters1[6]], 'Bx': [Bx], 'By': [By], 'Bz': [Bz], 'Bi': [p[4]], 'Ti': [p[5]], 'Pi': [p[6]], 'D': [KS[0]], 'p-value': [KS[1]], 'HWHM': [H], 'theta1': [t1], 'theta2': [t2], 'theta3': [t3], 'theta4': [t4]}
            df = pd.DataFrame(data_frame, columns = ['PosX', 'PosY', 'PosZ', 'B', 'T', 'P', 'Bx', 'By', 'Bz', 'Bi', 'Ti', 'Pi', 'D', 'p-value', 'HWHM', 'theta1', 'theta2', 'theta3', 'theta4'])
            df.to_csv (fo, mode = 'a', header = False)
    elif (regionT == 1 and regionP == 3):
        if (best_parameters1[5] <= deg2 and best_parameters1[5] >= deg1 and best_parameters1[6] <= deg6 and best_parameters1[6] >= deg5):
            data_frame = {'PosX': [posX], 'PosY': [posY], 'PosZ': [posZ], 'B': [best_parameters1[4]], 'T': [best_parameters1[5]], 'P': [best_parameters1[6]], 'Bx': [Bx], 'By': [By], 'Bz': [Bz], 'Bi': [p[4]], 'Ti': [p[5]], 'Pi': [p[6]], 'D': [KS[0]], 'p-value': [KS[1]], 'HWHM': [H], 'theta1': [t1], 'theta2': [t2], 'theta3': [t3], 'theta4': [t4]}
            df = pd.DataFrame(data_frame, columns = ['PosX', 'PosY', 'PosZ', 'B', 'T', 'P', 'Bx', 'By', 'Bz', 'Bi', 'Ti', 'Pi', 'D', 'p-value', 'HWHM', 'theta1', 'theta2', 'theta3', 'theta4'])
            df.to_csv (fo, mode = 'a', header = False)
    elif (regionT == 1 and regionP == 4):
        if (best_parameters1[5] <= deg2 and best_parameters1[5] >= deg1 and best_parameters1[6] <= deg8 and best_parameters1[6] >= deg7):
            data_frame = {'PosX': [posX], 'PosY': [posY], 'PosZ': [posZ], 'B': [best_parameters1[4]], 'T': [best_parameters1[5]], 'P': [best_parameters1[6]], 'Bx': [Bx], 'By': [By], 'Bz': [Bz], 'Bi': [p[4]], 'Ti': [p[5]], 'Pi': [p[6]], 'D': [KS[0]], 'p-value': [KS[1]], 'HWHM': [H], 'theta1': [t1], 'theta2': [t2], 'theta3': [t3], 'theta4': [t4]}
            df = pd.DataFrame(data_frame, columns = ['PosX', 'PosY', 'PosZ', 'B', 'T', 'P', 'Bx', 'By', 'Bz', 'Bi', 'Ti', 'Pi', 'D', 'p-value', 'HWHM', 'theta1', 'theta2', 'theta3', 'theta4'])
            df.to_csv (fo, mode = 'a', header = False)
    elif (regionT == 2 and regionP == 1):
        if (best_parameters1[5] <= deg4 and best_parameters1[5] >= deg3 and best_parameters1[6] <= deg2 and best_parameters1[6] >= deg1):
            data_frame = {'PosX': [posX], 'PosY': [posY], 'PosZ': [posZ], 'B': [best_parameters1[4]], 'T': [best_parameters1[5]], 'P': [best_parameters1[6]], 'Bx': [Bx], 'By': [By], 'Bz': [Bz], 'Bi': [p[4]], 'Ti': [p[5]], 'Pi': [p[6]], 'D': [KS[0]], 'p-value': [KS[1]], 'HWHM': [H], 'theta1': [t1], 'theta2': [t2], 'theta3': [t3], 'theta4': [t4]}
            df = pd.DataFrame(data_frame, columns = ['PosX', 'PosY', 'PosZ', 'B', 'T', 'P', 'Bx', 'By', 'Bz', 'Bi', 'Ti', 'Pi', 'D', 'p-value', 'HWHM', 'theta1', 'theta2', 'theta3', 'theta4'])
            df.to_csv (fo, mode = 'a', header = False)
    elif (regionT == 2 and regionP == 2):
        if (best_parameters1[5] <= deg4 and best_parameters1[5] >= deg3 and best_parameters1[6] <= deg4 and best_parameters1[6] >= deg3):
            data_frame = {'PosX': [posX], 'PosY': [posY], 'PosZ': [posZ], 'B': [best_parameters1[4]], 'T': [best_parameters1[5]], 'P': [best_parameters1[6]], 'Bx': [Bx], 'By': [By], 'Bz': [Bz], 'Bi': [p[4]], 'Ti': [p[5]], 'Pi': [p[6]], 'D': [KS[0]], 'p-value': [KS[1]], 'HWHM': [H], 'theta1': [t1], 'theta2': [t2], 'theta3': [t3], 'theta4': [t4]}
            df = pd.DataFrame(data_frame, columns = ['PosX', 'PosY', 'PosZ', 'B', 'T', 'P', 'Bx', 'By', 'Bz', 'Bi', 'Ti', 'Pi', 'D', 'p-value', 'HWHM', 'theta1', 'theta2', 'theta3', 'theta4'])
            df.to_csv (fo, mode = 'a', header = False)
    elif (regionT == 2 and regionP == 3):
        if (best_parameters1[5] <= deg4 and best_parameters1[5] >= deg3 and best_parameters1[6] <= deg6 and best_parameters1[6] >= deg5):
            data_frame = {'PosX': [posX], 'PosY': [posY], 'PosZ': [posZ], 'B': [best_parameters1[4]], 'T': [best_parameters1[5]], 'P': [best_parameters1[6]], 'Bx': [Bx], 'By': [By], 'Bz': [Bz], 'Bi': [p[4]], 'Ti': [p[5]], 'Pi': [p[6]], 'D': [KS[0]], 'p-value': [KS[1]], 'HWHM': [H], 'theta1': [t1], 'theta2': [t2], 'theta3': [t3], 'theta4': [t4]}
            df = pd.DataFrame(data_frame, columns = ['PosX', 'PosY', 'PosZ', 'B', 'T', 'P', 'Bx', 'By', 'Bz', 'Bi', 'Ti', 'Pi', 'D', 'p-value', 'HWHM', 'theta1', 'theta2', 'theta3', 'theta4'])
            df.to_csv (fo, mode = 'a', header = False)
    elif (regionT == 2 and regionP == 4):
        if (best_parameters1[5] <= deg4 and best_parameters1[5] >= deg3 and best_parameters1[6] <= deg8 and best_parameters1[6] >= deg7):
            data_frame = {'PosX': [posX], 'PosY': [posY], 'PosZ': [posZ], 'B': [best_parameters1[4]], 'T': [best_parameters1[5]], 'P': [best_parameters1[6]], 'Bx': [Bx], 'By': [By], 'Bz': [Bz], 'Bi': [p[4]], 'Ti': [p[5]], 'Pi': [p[6]], 'D': [KS[0]], 'p-value': [KS[1]], 'HWHM': [H], 'theta1': [t1], 'theta2': [t2], 'theta3': [t3], 'theta4': [t4]}
            df = pd.DataFrame(data_frame, columns = ['PosX', 'PosY', 'PosZ', 'B', 'T', 'P', 'Bx', 'By', 'Bz', 'Bi', 'Ti', 'Pi', 'D', 'p-value', 'HWHM', 'theta1', 'theta2', 'theta3', 'theta4'])
            df.to_csv (fo, mode = 'a', header = False)
    elif (regionT == 3 and regionP == 1):
        if (best_parameters1[5] <= deg6 and best_parameters1[5] >= deg5 and best_parameters1[6] <= deg2 and best_parameters1[6] >= deg1):
            data_frame = {'PosX': [posX], 'PosY': [posY], 'PosZ': [posZ], 'B': [best_parameters1[4]], 'T': [best_parameters1[5]], 'P': [best_parameters1[6]], 'Bx': [Bx], 'By': [By], 'Bz': [Bz], 'Bi': [p[4]], 'Ti': [p[5]], 'Pi': [p[6]], 'D': [KS[0]], 'p-value': [KS[1]], 'HWHM': [H], 'theta1': [t1], 'theta2': [t2], 'theta3': [t3], 'theta4': [t4]}
            df = pd.DataFrame(data_frame, columns = ['PosX', 'PosY', 'PosZ', 'B', 'T', 'P', 'Bx', 'By', 'Bz', 'Bi', 'Ti', 'Pi', 'D', 'p-value', 'HWHM', 'theta1', 'theta2', 'theta3', 'theta4'])
            df.to_csv (fo, mode = 'a', header = False)
    elif (regionT == 3 and regionP == 2):
        if (best_parameters1[5] <= deg6 and best_parameters1[5] >= deg5 and best_parameters1[6] <= deg4 and best_parameters1[6] >= deg3):
            data_frame = {'PosX': [posX], 'PosY': [posY], 'PosZ': [posZ], 'B': [best_parameters1[4]], 'T': [best_parameters1[5]], 'P': [best_parameters1[6]], 'Bx': [Bx], 'By': [By], 'Bz': [Bz], 'Bi': [p[4]], 'Ti': [p[5]], 'Pi': [p[6]], 'D': [KS[0]], 'p-value': [KS[1]], 'HWHM': [H], 'theta1': [t1], 'theta2': [t2], 'theta3': [t3], 'theta4': [t4]}
            df = pd.DataFrame(data_frame, columns = ['PosX', 'PosY', 'PosZ', 'B', 'T', 'P', 'Bx', 'By', 'Bz', 'Bi', 'Ti', 'Pi', 'D', 'p-value', 'HWHM', 'theta1', 'theta2', 'theta3', 'theta4'])
            df.to_csv (fo, mode = 'a', header = False)
    elif (regionT == 3 and regionP == 3):
        if (best_parameters1[5] <= deg6 and best_parameters1[5] >= deg5 and best_parameters1[6] <= deg6 and best_parameters1[6] >= deg5):
            data_frame = {'PosX': [posX], 'PosY': [posY], 'PosZ': [posZ], 'B': [best_parameters1[4]], 'T': [best_parameters1[5]], 'P': [best_parameters1[6]], 'Bx': [Bx], 'By': [By], 'Bz': [Bz], 'Bi': [p[4]], 'Ti': [p[5]], 'Pi': [p[6]], 'D': [KS[0]], 'p-value': [KS[1]], 'HWHM': [H], 'theta1': [t1], 'theta2': [t2], 'theta3': [t3], 'theta4': [t4]}
            df = pd.DataFrame(data_frame, columns = ['PosX', 'PosY', 'PosZ', 'B', 'T', 'P', 'Bx', 'By', 'Bz', 'Bi', 'Ti', 'Pi', 'D', 'p-value', 'HWHM', 'theta1', 'theta2', 'theta3', 'theta4'])
            df.to_csv (fo, mode = 'a', header = False)
    elif (regionT == 3 and regionP == 4):
        if (best_parameters1[5] <= deg6 and best_parameters1[5] >= deg5 and best_parameters1[6] <= deg8 and best_parameters1[6] >= deg7):
            data_frame = {'PosX': [posX], 'PosY': [posY], 'PosZ': [posZ], 'B': [best_parameters1[4]], 'T': [best_parameters1[5]], 'P': [best_parameters1[6]], 'Bx': [Bx], 'By': [By], 'Bz': [Bz], 'Bi': [p[4]], 'Ti': [p[5]], 'Pi': [p[6]], 'D': [KS[0]], 'p-value': [KS[1]], 'HWHM': [H], 'theta1': [t1], 'theta2': [t2], 'theta3': [t3], 'theta4': [t4]}
            df = pd.DataFrame(data_frame, columns = ['PosX', 'PosY', 'PosZ', 'B', 'T', 'P', 'Bx', 'By', 'Bz', 'Bi', 'Ti', 'Pi', 'D', 'p-value', 'HWHM', 'theta1', 'theta2', 'theta3', 'theta4'])
            df.to_csv (fo, mode = 'a', header = False)
    elif (regionT == 4 and regionP == 1):
        if (best_parameters1[5] <= deg8 and best_parameters1[5] >= deg7 and best_parameters1[6] <= deg2 and best_parameters1[6] >= deg1):
            data_frame = {'PosX': [posX], 'PosY': [posY], 'PosZ': [posZ], 'B': [best_parameters1[4]], 'T': [best_parameters1[5]], 'P': [best_parameters1[6]], 'Bx': [Bx], 'By': [By], 'Bz': [Bz], 'Bi': [p[4]], 'Ti': [p[5]], 'Pi': [p[6]], 'D': [KS[0]], 'p-value': [KS[1]], 'HWHM': [H], 'theta1': [t1], 'theta2': [t2], 'theta3': [t3], 'theta4': [t4]}
            df = pd.DataFrame(data_frame, columns = ['PosX', 'PosY', 'PosZ', 'B', 'T', 'P', 'Bx', 'By', 'Bz', 'Bi', 'Ti', 'Pi', 'D', 'p-value', 'HWHM', 'theta1', 'theta2', 'theta3', 'theta4'])
            df.to_csv (fo, mode = 'a', header = False)
    elif (regionT == 4 and regionP == 2):
        if (best_parameters1[5] <= deg8 and best_parameters1[5] >= deg7 and best_parameters1[6] <= deg4 and best_parameters1[6] >= deg3):
            data_frame = {'PosX': [posX], 'PosY': [posY], 'PosZ': [posZ], 'B': [best_parameters1[4]], 'T': [best_parameters1[5]], 'P': [best_parameters1[6]], 'Bx': [Bx], 'By': [By], 'Bz': [Bz], 'Bi': [p[4]], 'Ti': [p[5]], 'Pi': [p[6]], 'D': [KS[0]], 'p-value': [KS[1]], 'HWHM': [H], 'theta1': [t1], 'theta2': [t2], 'theta3': [t3], 'theta4': [t4]}
            df = pd.DataFrame(data_frame, columns = ['PosX', 'PosY', 'PosZ', 'B', 'T', 'P', 'Bx', 'By', 'Bz', 'Bi', 'Ti', 'Pi', 'D', 'p-value', 'HWHM', 'theta1', 'theta2', 'theta3', 'theta4'])
            df.to_csv (fo, mode = 'a', header = False)
    elif (regionT == 4 and regionP == 3):
        if (best_parameters1[5] <= deg8 and best_parameters1[5] >= deg7 and best_parameters1[6] <= deg6 and best_parameters1[6] >= deg5):
            data_frame = {'PosX': [posX], 'PosY': [posY], 'PosZ': [posZ], 'B': [best_parameters1[4]], 'T': [best_parameters1[5]], 'P': [best_parameters1[6]], 'Bx': [Bx], 'By': [By], 'Bz': [Bz], 'Bi': [p[4]], 'Ti': [p[5]], 'Pi': [p[6]], 'D': [KS[0]], 'p-value': [KS[1]], 'HWHM': [H], 'theta1': [t1], 'theta2': [t2], 'theta3': [t3], 'theta4': [t4]}
            df = pd.DataFrame(data_frame, columns = ['PosX', 'PosY', 'PosZ', 'B', 'T', 'P', 'Bx', 'By', 'Bz', 'Bi', 'Ti', 'Pi', 'D', 'p-value', 'HWHM', 'theta1', 'theta2', 'theta3', 'theta4'])
            df.to_csv (fo, mode = 'a', header = False)
    elif (regionT == 4 and regionP == 4):
        if (best_parameters1[5] <= deg8 and best_parameters1[5] >= deg7 and best_parameters1[6] <= deg8 and best_parameters1[6] >= deg7):
            data_frame = {'PosX': [posX], 'PosY': [posY], 'PosZ': [posZ], 'B': [best_parameters1[4]], 'T': [best_parameters1[5]], 'P': [best_parameters1[6]], 'Bx': [Bx], 'By': [By], 'Bz': [Bz], 'Bi': [p[4]], 'Ti': [p[5]], 'Pi': [p[6]], 'D': [KS[0]], 'p-value': [KS[1]], 'HWHM': [H], 'theta1': [t1], 'theta2': [t2], 'theta3': [t3], 'theta4': [t4]}
            df = pd.DataFrame(data_frame, columns = ['PosX', 'PosY', 'PosZ', 'B', 'T', 'P', 'Bx', 'By', 'Bz', 'Bi', 'Ti', 'Pi', 'D', 'p-value', 'HWHM', 'theta1', 'theta2', 'theta3', 'theta4'])
            df.to_csv (fo, mode = 'a', header = False)
    else:
        print ("Invalid region, please select between 1 and 4")
 
    
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


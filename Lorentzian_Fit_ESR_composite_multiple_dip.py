####################### IMPORTING REQUIRED MODULES ######################
import pandas as pd
import numpy as np
import matplotlib.pyplot as pl
from scipy.optimize import leastsq # Levenberg-Marquadt Algorithm #
from scipy import optimize
############################# LOADING DATA ##############################

Data = pd.read_csv('/home/dwi/ownCloud/Anlabshared/Dwi/Data and Analysis/160627/2016_06_27_MagPar_pos5.csv')
Data.columns = ['Microwave Frequency [GHz]', 'Intensity']
x = Data['Microwave Frequency [GHz]']
y = Data['Intensity']/Data['Intensity'].max()
#HWHM = 3550000
#dip = 0.04
########################### DEFINING FUNCTIONS ##########################

#p = [8.5*10**-4, 0, 0, 7000000, 0.023]  # [B, theta, phi, HWHM, dip]
 
#def lorentzian(x,p):
    #numerator = p[3]**2
    #denumerator = ( x - (2.87*10**9 + (28*10**9 * p[0] * np.sqrt((np.sin (p[1]) * np.cos(p[2]))**2 + (np.sin(p[1]) * np.sin(p[2]))**2 + (np.cos(p[2]))**2))))**2 + p[3]**2
    #numerator2 =  p[3]**2
    #denumerator2 = ( x - (2.87*10**9 - (28*10**9 * p[0] * np.sqrt((np.sin(p[1]) * np.cos(p[2]))**2 + (np.sin(p[1]) * np.sin(p[2]))**2 + (np.cos(p[2]))**2))))**2 + p[3]**2
    #return 0.995 - p[4]*((numerator/denumerator) + (numerator2/denumerator2))

def lorentzian2(x,p):
    numerator3 = p[2]**2
    denumerator3 = (x - (2.87*10**9 - (28*10**9 * p[0] * np.cos(p[1]))))**2 + p[2]**2
    return 0.995 - p[3]*(numerator3/denumerator3)
    
#def lorentzian3(x,p):    
    #numerator4 =  p[3]**2
    #denumerator4 = ( x - (2.87*10**9 - (28*10**9 * p[0] * np.sqrt((np.sin(p[1]) * np.cos(p[2]))**2 + (np.sin(p[1]) * np.sin(p[2]))**2 + (np.cos(p[2]))**2))))**2 + p[3]**2
    #return 0.995 - p[4]*(numerator4/denumerator4)

#def residuals(p,y,x):
    #err = y - lorentzian(x,p)
    #return err

def residuals2(p,y,x):
    err = y - lorentzian2(x,p)
    return err
    
#def residuals3(p,y,x):
    #err = y - lorentzian3(x,p)    
############################## FITTING DATA ## ###########################

# initial values #
p = [5*10**-4, 10,3550000, 0.02]  # [B, theta,HWHM, dip]

#Dip search
#Dip = y.min(2.8*10**9, 2.85*10**9)

# optimization # 
#pbest = leastsq(residuals, p ,args=(y,x),full_output=1)
#best_parameters1 = pbest[0]
#y[70:143], x[70:143]
pbest2 = leastsq(residuals2, p ,args=(y[147:162], x[147:162]),full_output=1)
best_parameters2 = pbest2[0]

#pbest3 = leastsq(residuals3, p ,args=(y[0:81], x[0:81]),full_output=1)
#best_parameters3 = pbest3[0]

# fit to data #
#fit = lorentzian(x,best_parameters1)
fit2 = lorentzian2(x,best_parameters2)
#fit3 = lorentzian3(x,best_parameters3)
#comp = 1 - fit + fit2
############################## PLOTTING #################################
pl.rc('text', usetex=True)
pl.rc('font', family='serif')

pl.plot(x,y, lw = 1.5, color='g', label='Measured ESR spectrum')
#pl.plot(x, fit, lw=1.5, color = 'r', label = r"Least-square fit")
pl.plot( x, fit2, lw=1.5, color = 'k')
#pl.plot( x, fit3, lw=1.5, color = 'b')

#pl.plot(x, comp)
#pl.plot(x, lorentzian(x,p), x, lorentzian2(x,p))
pl.xlabel('Microwave frequency (GHz)', fontsize=18)
pl.ylabel('Intensity (norm.)', fontsize=18)
#pl.plot(x[143:222], y[143:222], color='r')
pl.legend()
pl.show()
#print best_parameters1
#print Dip
print best_parameters2
#print Data
print x[147:162]
#print residuals(p,y, x)
#print lorentzian (x, p)

import pandas as pd
import numpy as np
import matplotlib.pyplot as pl

Data = pd.read_csv ('/home/dwi/ownCloud/Anlabshared/Dwi/Data and Analysis/160718/2016_07_18_7.csv')
Data.columns = ['Microwave frequency', 'Intensity']
x = np.array(Data['Microwave frequency'])
y = np.array(Data['Intensity']/Data['Intensity'].max())

################# SIMULATION PARAMETERS #########################################
simp = [7*10**6, 0.03, 5.0*10**(-3), 0.0, 109.5, 109.5, 109.5] 

################# INDIVIDUAL ESR SPECTRA ########################################
num = simp[0]**2
num2 = simp[0]**2
num3 = simp[0]**2
num4 = simp[0]**2
num5 = simp[0]**2
num6 = simp[0]**2
num7 = simp[0]**2
num8 = simp[0]**2
denum = ( x - (2.87*10**9 - (28*10**9 * simp[2] * np.cos(simp[3]*np.pi/180.0))))**2 + simp[0]**2
denum2 = ( x - (2.87*10**9 - (28*10**9 * simp[2] * np.cos(simp[4]*np.pi/180.0))))**2 + simp[0]**2
denum3 = ( x - (2.87*10**9 - (28*10**9 * simp[2] * np.cos(simp[5]*np.pi/180.0))))**2 + simp[0]**2
denum4 = ( x - (2.87*10**9 - (28*10**9 * simp[2] * np.cos(simp[6]*np.pi/180.0))))**2 + simp[0]**2
denum5 = ( x - (2.87*10**9 + (28*10**9 * simp[2] * np.cos(simp[6]*np.pi/180.0))))**2 + simp[0]**2
denum6 = ( x - (2.87*10**9 + (28*10**9 * simp[2] * np.cos(simp[5]*np.pi/180.0))))**2 + simp[0]**2
denum7 = ( x - (2.87*10**9 + (28*10**9 * simp[2] * np.cos(simp[4]*np.pi/180.0))))**2 + simp[0]**2
denum8 = ( x - (2.87*10**9 + (28*10**9 * simp[2] * np.cos(simp[3]*np.pi/180.0))))**2 + simp[0]**2
Lor = (simp[1]*(num/denum)/np.pi)
Lor2 = (simp[1]*(num2/denum2)/np.pi)
Lor3 = (simp[1]*(num3/denum3)/np.pi)
Lor4 = (simp[1]*(num4/denum4)/np.pi)
Lor5 = (simp[1]*(num5/denum5)/np.pi)
Lor6 = (simp[1]*(num6/denum6)/np.pi)
Lor7 = (simp[1]*(num7/denum7)/np.pi)
Lor8 = (simp[1]*(num8/denum8)/np.pi)

Lorentzian = 1 - (Lor + Lor2 + Lor3 + Lor4 + Lor5 + Lor6 + Lor7 + Lor8)
################################################################################
pl.rc('text', usetex=True)
pl.rc('font', family='serif')
pl.plot (x, y, lw = 1.5, label = 'Experiment')
#pl.plot (x[dipind],y[dipind], 'ro')
#pl.plot(x, fit, lw=1.5, color = 'r', label = r"Least-square fit") 
#pl.plot (x, Lor, color = 'g', label = r"NV1")
#pl.plot (x, Lor8, color = 'g')
#pl.plot (x, Lor2, color = 'y', label = r"NV2")
#pl.plot (x, Lor7, color = 'y')
#pl.plot (x, Lor3, color = 'c', label = r"NV3")
#pl.plot (x, Lor6, color = 'c')
#pl.plot (x, Lor4, color = 'm', label = r"NV4")
#pl.plot (x, Lor5, color = 'm')
pl.plot (x, Lorentzian, color = "r", lw = 1.5, label = r"Simulation $B_0$ = 5 mT, $\theta_1$ = 0.0$^o$, $\theta_{2, 3, 4}$ = 109.5$^o$")
pl.xlabel('Microwave frequency (GHz)', fontsize = 18)
pl.ylabel('Intensity (norm.)', fontsize = 18)
#pl.title (r"Position 5")
pl.legend(loc = 4)
pl.show()

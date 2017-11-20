import numpy as np
import matplotlib.pyplot as pl
import scipy
from scipy import stats
x,y  = np.loadtxt('2016_05_06_Yellow_Diamond_10.csv', delimiter=',', unpack=True)
u,v = np.loadtxt('2016_05_06_Yellow_Diamond_11.csv', delimiter=',', unpack=True)
a,b = np.loadtxt('2016_05_06_Yellow_Diamond_12.csv', delimiter=',', unpack=True)
c,d = np.loadtxt('2016_05_06_Yellow_Diamond_13.csv', delimiter=',', unpack=True)
e,f = np.loadtxt('2016_05_06_Yellow_Diamond_14.csv', delimiter=',', unpack=True)
g,h = np.loadtxt('2016_05_06_Yellow_Diamond_15.csv', delimiter=',', unpack=True)
i,j = np.loadtxt('2016_05_06_Yellow_Diamond_16.csv', delimiter=',', unpack=True)
k,l = np.loadtxt('2016_05_06_Yellow_Diamond_17.csv', delimiter=',', unpack=True)
m,n = np.loadtxt('2016_05_06_Yellow_Diamond_18.csv', delimiter=',', unpack=True)
o,p = np.loadtxt('2016_05_06_Yellow_Diamond_19.csv', delimiter=',', unpack=True)
q,r = np.loadtxt('2016_05_06_Yellow_Diamond_20.csv', delimiter=',', unpack=True)
mi1 = min(y)
ma1 = max (y)
mi2 = min(v)
ma2 = max(v)
pl.plot (x, y, linewidth=2.0, color='Blue', label='singe-mode FO' )
pl.plot(u, v, linewidth=2.0, color='Red', label='multi-mode FO')
pl.plot(a, b, linewidth=2.0, color='Red', label='multi-mode FO')
pl.plot(c, d, linewidth=2.0, color='Red', label='multi-mode FO')
pl.plot(e, f, linewidth=2.0, color='Red', label='multi-mode FO')
pl.plot(g, h, linewidth=2.0, color='Red', label='multi-mode FO')
pl.plot(i, j, linewidth=2.0, color='Red', label='multi-mode FO')
pl.plot(k, l, linewidth=2.0, color='Red', label='multi-mode FO')
pl.plot(m, n, linewidth=2.0, color='Red', label='multi-mode FO')
pl.plot(o, p, linewidth=2.0, color='Red', label='multi-mode FO')
pl.plot(q, r, linewidth=2.0, color='Red', label='multi-mode FO')
pl.title ('ESR Signal of Yellow Diamond NV Center')
pl.xlabel('Microwave Frequency [Hz]')
pl.ylabel('Photon Counts per second')
pl.show()
snr1 = scipy.stats.signaltonoise(y)
snr2 = scipy.stats.signaltonoise(v)
print "Signal 1;"
print "Minimum photon counts of signal 1=" + str(mi1)
print "Maximum photon counts of signal 1=" + str(ma1)
print "Signal to noise ratio of signal 1=" + str(snr1)
print "Signal 2;"
print "Minimum photon counts of signal 2=" + str(mi2)
print "Maximum photon counts of signal 2=" + str(ma2)
print "Signal to noise ratio of signal 2=" + str(snr2)

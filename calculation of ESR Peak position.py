import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
###########################################################################
####   This program calculate ESR peak position using the angles       ####
####           between NV axis and magnetic field (theta)              ####
####         This calculation is based on next Hamiltonian             ####
####            H = hbar * D * S_z^ 2 + g * mu_B * B * S               ####
###########################################################################

####save figure############################################################
Path = 'D:/eXperiment_in_Anlab/otherexperiment/' #save directry
file_name = 'cubic.jpeg' # file name + extension(eps, png, jpeg or etc.)
save = Path + file_name
####parameter##############################################################
theta_deg = np.array([20, 40, 60, 80]) #input the values of angle(deg)
B_range = [0, 300] #input the magnetic fied range (Oe)
B_span = 1 #input the magnetic field span (Oe)
####figure profile#########################################################
#plt.figure(figsize=(6, 4)) #figure size
plt.rcParams['font.family'] = 'Times New Roman' #font
plt.xlim(B_range)
plt.ylim([2, 4]) 
ttl = 'ESR peak position'
plt.title (ttl , fontsize =20)
plt.ylabel('Frequency (GHz)', fontsize = 18)
plt.xlabel('Magnetic field (Oe)', fontsize = 18)

###constant mu_B, D
mu_B = 13.995384 / 10000#Bohr magnetoron (GHz/Oe)
D = 2.87 #Zero filed splitting (GHz)

###calculation in each theta
for i, theta in enumerate(theta_deg):
    theta_rad = theta / 180 * np.pi
    SIN = np.sin(theta_rad)
    
    listnega = []
    listpos = []
    field = range(B_range[0], B_range[1], B_span)
    ######calculation in each magnetic field
    for j in range(B_range[0], B_range[1], B_span):
        beta = 2 * j * mu_B
        ###coefficient of cubicequation
        a = 1
        b = -2 * D
        c = D ** 2 - beta ** 2
        d = beta ** 2 * D * SIN ** 2
        coef = [a, b, c, d] #coefficient of cubicequation
        solve = np.roots(coef) #Solution of cubicequation
        ##calulate the energy gap
        solvenega = solve[1]-solve[2] #m_s=0 to m_x=-1
        solvepos = solve[0]-solve[2] #m_s=0 to m_x=+1
        listnega.append(solvenega)
        listpos.append(solvepos)
        
    c = cm.gist_rainbow((float(i) )/ (len(theta_deg))) #color of plot
    plt.plot(field,listnega, color = c, label = str(theta) + 'deg')
    plt.plot(field,listpos, color = c)
    plt.legend(loc = 'upper left')

plt.savefig(save)
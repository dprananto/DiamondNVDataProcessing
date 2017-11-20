import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
###########################################################################
####   This program calculate ESR peak position using the angles       ####
####           between NV axis and magnetic field (theta)              ####
####         This calculation is based on next Hamiltonian             ####
####            H = hbar * D * S_z^ 2 + g * mu_B * B * S               ####
###########################################################################

####select diamond type [(100), (110) or (111)]############################
Diamondtype = (100)
####parameter##############################################################
B = 200. #magnetic field (Oe)
theta_range = [0, 360] #input in-plane angle of Magnetic field (deg)
theta_span = 2
Phideg = 90. #input out of plane angle of Magnetic field (deg)
####save figure############################################################
Path = 'D:/eXperiment_in_Anlab/otherexperiment/' #save directry
file_name = 'diamond100angledep.jpeg' # file name + extension(eps, png, jpeg or etc.)
save = Path + file_name
####figure profile#########################################################
#plt.figure(figsize=(6, 4)) #figure size
plt.rcParams['font.family'] = 'Times New Roman' #font
plt.xlim(theta_range)
#plt.ylim([2, 4]) 
ttl = 'ESR peak position'
plt.title (ttl , fontsize =20)
plt.ylabel('Frequency (GHz)', fontsize = 18)
plt.xlabel('theta (deg)', fontsize = 18)

###constant mu_B, D
mu_B = 13.995384 / 10000#Bohr magnetoron (GHz/Oe)
D = 2.87 #Zero filed splitting (GHz)
beta = 2 * B * mu_B
###calculation in each theta
if Diamondtype == (100):
    NV1 = np.array([[1, 1, 1]])
    NV2 = np.array([[1, -1, -1]])
    NV3 = np.array([[-1, 1, -1]])
    NV4 = np.array([[-1, -1, 1]])
if Diamondtype == (110):
    NV1 = np.array([[-1, - np.sqrt(2), 0]])
    NV2 = np.array([[-1, np.sqrt(2), 0]])
    NV3 = np.array([[1, 0, np.sqrt(2)]])
    NV4 = np.array([[1, 0, -np.sqrt(2)]])
if Diamondtype == (111):
    NV1 = np.array([[0, 0, np.sqrt(3)]])
    NV2 = np.array([[np.sqrt(2), - np.sqrt(2/3), - np.sqrt(1/3)]])
    NV3 = np.array([[0, np.sqrt(8/3), - np.sqrt(1/3)]])
    NV4 = np.array([[- np.sqrt(2), - np.sqrt(2/3), - np.sqrt(1/3)]])
NV = [NV1/ np.linalg.norm(NV1), NV2/ np.linalg.norm(NV2), NV3/ np.linalg.norm(NV3), NV4/ np.linalg.norm(NV4)]

Phi = Phideg * np.pi / 180
T_angle = range(theta_range[0], theta_range[1], theta_span)

for i, num in enumerate(NV): ##loop for each NV axis (NV1~4)
    listnega = []
    listpos = []
    ######calculation for each angle
    for theta in range(theta_range[0], theta_range[1], theta_span):
        Theta = theta / 180 * np.pi
        H = np.array([[np.cos(Theta) * np.sin(Phi), np.sin(Theta) * np.sin(Phi),  np.cos(Phi)]])
        SIN = np.sin(np.arccos(np.inner(H, num) ) ) # caluculation for the angle between NV and H
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
        
    c = cm.gist_rainbow((float(i) )/ 4) #color of plot
    plt.plot(T_angle,listnega, color = c, label = 'NV' + str(i + 1))
    plt.plot(T_angle,listpos, color = c)
    plt.legend(loc = 'upper left')

plt.savefig(save)
        
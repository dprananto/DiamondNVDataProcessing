# -*- coding: utf-8 -*-
"""
Created on Thu Oct  5 15:17:40 2017

@author: saitohlab-K

マイクロ波反射測定（カンバラプログラムのデータ解析）（仮）
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os


Path = '/home/dwi/ownCloud/Anlabshared/Dwi/researchdata/171101SSENVC_NVA/171101_SSE_NVC/'
###calculating magnetic field from applied voltage
a = 804.5
b = 38
def Magnetic(x,a,b):
    return a*x + b     

#reading file　
Base_H = 0.0000 #selecting reference file
File_name = 'number0_I='
folder_name = 'diff_by' + str(int(Magnetic(Base_H,a,b))) + 'Oe'#フォルダ名の作成に利用 差分を求めるのに用いた基準
Base_file = File_name + str("%.4f" %Base_H) + '.csv'
span = 0.01
file_num = 101

TARGET_DIR = Path + folder_name + '/'
if not os.path.isdir(TARGET_DIR):
    os.mkdir(TARGET_DIR) 
TARGET_DIR1 = TARGET_DIR + 'png/'
if not os.path.isdir(TARGET_DIR1):
    os.mkdir(TARGET_DIR1) 
TARGET_DIR2 = TARGET_DIR + 'eps/'
if not os.path.isdir(TARGET_DIR2):
    os.mkdir(TARGET_DIR2) 

Data = pd.read_csv (Path + Base_file, header=None)
num_lines = sum(1 for line in open(Path + Base_file))
Data.columns = ['freq', 'S11','phase','real','imaginal']
Base_x = np.array(Data['freq'])
Base_y = np.array(Data['S11'])

Matrix =[]
for i in range(0, file_num):
    H = i * span
    H_Oe = int(Magnetic(H,a,b))
    fig_name = str(H_Oe) +'Oe'
    file =File_name + str("%.4f" %  H) +'.csv'
    Data = pd.read_csv (Path + file, header=None)
    num_lines = sum(1 for line in open(Path + file))
    Data.columns = ['freq', 'S11','phase','real','imaginal']
    sig_x = np.array(Data['freq'])
    sig_y = np.array(Data['S11'])
    calc_y = sig_y - Base_y 
    Matrix.append(calc_y)
    ####Figure###
    plt.plot (sig_x, calc_y, color = 'r')#,label = 'Experiment')
    #label = Figure_title 
    plt.ylabel('S11(a.u.)', fontsize = 20)
    plt.xlabel('Frequency GHz', fontsize = 20)
    #plt.xlim(0,time_scale) 
    plt.title (fig_name)
    #plt.legend()
    #plt.savefig(TARGET_DIR2 + File_name + str("%.4f" %  H)  + '.eps')
    plt.savefig(TARGET_DIR1 + File_name + fig_name  + '.png')
    plt.clf()

np.savetxt(TARGET_DIR + 'H_dependence.csv', Matrix, delimiter=",")


####2次元マップの作製
Z1 = Matrix.T
x1 = Base_x
y1 = np.arange(b, a+b+0.1, span*a)
X1, Y1 = np.meshgrid(x1, y1)
plt.figure(figsize=(8, 6))
plt.axis([2*10**9,4*10**9,40,850])
plt.pcolor(X1, Y1, Z1,cmap='RdBu' ,vmax =5, vmin =-15)
plt.colorbar()
plt.xlabel('Microwave frequency (GHz)', fontsize = 12)
plt.ylabel('Magnetic field (Oe)', fontsize = 12)
plt.xticks(fontsize = 12)
plt.yticks(fontsize = 12)
#色の指定
plt.jet()

save = TARGET_DIR + 'map_short.png'
plt.savefig(save)
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 23 23:39:55 2018

@author: 朱高明
"""

import numpy as np
import pandas as pd
import math
from tqdm import tqdm

def rot(E01, E02, E03, E11, E12, E13):
#    E01 = float(E01)
#    E11 = float(E11)
#    if float(E02) < 180:
#        E03 = (float(E03) + 90) % 60
#    else:
#        E03 = 60 - (float(E03) + 90) % 60
#    
#    if float(E12) < 180:
#        E13 = (float(E13) + 90) % 60
#    else:
#        E13 = 60 - (float(E13) + 90) % 60
#        
#    E02 = float(E02) % 180
#    E12 = float(E12) % 180

    Euler01rad = E01 * math.pi / 180
    Euler02rad = E02 * math.pi / 180
    Euler03rad = E03 * math.pi / 180
    Euler11rad = E11 * math.pi / 180
    Euler12rad = E12 * math.pi / 180
    Euler13rad = E13 * math.pi / 180
    
    Euler0 = np.zeros((12, 3))
    Euler1 = np.zeros((12, 3))
    AngleM = np.zeros((12, 12))
    
    gATM = np.matrix([[0.0, 0.0, 0.0],  [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
    
    i = 0
    while i < len(Euler0):
        if i < 6:
            Euler0[i, 0] = Euler01rad
            Euler0[i, 1] = Euler02rad
            Euler0[i, 2] = Euler03rad + i * 60 * math.pi / 180
            Euler1[i, 0] = Euler11rad
            Euler1[i, 1] = Euler12rad
            Euler1[i, 2] = Euler13rad + i * 60 * math.pi / 180
        else:
            Euler0[i, 0] = Euler01rad + 180 * math.pi / 180
            Euler0[i, 1] = 180 * math.pi / 180 - Euler02rad
            Euler0[i, 2] = -Euler03rad + (i - 6) * 60 * math.pi / 180
            Euler1[i, 0] = Euler11rad + 180 * math.pi / 180
            Euler1[i, 1] = 180 * math.pi / 180 - Euler12rad
            Euler1[i, 2] = -Euler13rad + (i - 6) * 60 * math.pi / 180
        i += 1
    
    i = 0
    while i < len(Euler0):
        gATM[0, 0] = math.cos(Euler0[i, 0]) * math.cos(Euler0[i, 2]) - math.sin(Euler0[i, 0]) * math.sin(Euler0[i, 2]) * math.cos(Euler0[i, 1])
        gATM[0, 1] = -math.cos(Euler0[i, 0]) * math.sin(Euler0[i, 2]) - math.sin(Euler0[i, 0]) * math.cos(Euler0[i, 2]) * math.cos(Euler0[i, 1])
        gATM[0, 2] = math.sin(Euler0[i, 0]) * math.sin(Euler0[i, 1])
        gATM[1, 0] = math.sin(Euler0[i, 0]) * math.cos(Euler0[i, 2]) + math.cos(Euler0[i, 0]) * math.sin(Euler0[i, 2]) * math.cos(Euler0[i, 1])
        gATM[1, 1] = - math.sin(Euler0[i, 0]) * math.sin(Euler0[i, 2]) + math.cos(Euler0[i, 0]) * math.cos(Euler0[i, 2]) * math.cos(Euler0[i, 1])
        gATM[1, 2] = - math.cos(Euler0[i, 0]) * math.sin(Euler0[i, 1])
        gATM[2, 0] = math.sin(Euler0[i, 2]) * math.sin(Euler0[i, 1])
        gATM[2, 1] = math.cos(Euler0[i, 2]) * math.sin(Euler0[i, 1])
        gATM[2, 2] = math.cos(Euler0[i, 1])
        
        gBTM = np.matrix([[0.0, 0.0, 0.0],  [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]])    

        j = 0
        while j < len(Euler0):
            gBTM[0, 0] = math.cos(Euler1[j, 0]) * math.cos(Euler1[j, 2]) - math.sin(Euler1[j, 0]) * math.sin(Euler1[j, 2]) * math.cos(Euler1[j, 1])
            gBTM[0, 1] = -math.cos(Euler1[j, 0]) * math.sin(Euler1[j, 2]) - math.sin(Euler1[j, 0]) * math.cos(Euler1[j, 2]) * math.cos(Euler1[j, 1])
            gBTM[0, 2] = math.sin(Euler1[j, 0]) * math.sin(Euler1[j, 1])
            gBTM[1, 0] = math.sin(Euler1[j, 0]) * math.cos(Euler1[j, 2]) + math.cos(Euler1[j, 0]) * math.sin(Euler1[j, 2]) * math.cos(Euler1[j, 1])
            gBTM[1, 1] = - math.sin(Euler1[j, 0]) * math.sin(Euler1[j, 2]) + math.cos(Euler1[j, 0]) * math.cos(Euler1[j, 2]) * math.cos(Euler1[j, 1])
            gBTM[1, 2] = - math.cos(Euler1[j, 0]) * math.sin(Euler1[j, 1])
            gBTM[2, 0] = math.sin(Euler1[j, 2]) * math.sin(Euler1[j, 1])
            gBTM[2, 1] = math.cos(Euler1[j, 2]) * math.sin(Euler1[j, 1])
            gBTM[2, 2] = math.cos(Euler1[j, 1])
        
#        gATM = np.matrix([[gAT11, gAT12, gAT13],  [gAT21, gAT22, gAT23], [gAT31, gAT32, gAT33]])
#        gBTM = np.matrix([[gBT11, gBT12, gBT13],  [gBT21, gBT22, gBT23], [gBT31, gBT32, gBT33]])
            gB = gBTM.T
            
            gAB11 = float(gB[0,:] * gATM[:,0])
            gAB12 = float(gB[0,:] * gATM[:,1])
            gAB13 = float(gB[0,:] * gATM[:,2])
            gAB21 = float(gB[1,:] * gATM[:,0])
            gAB22 = float(gB[1,:] * gATM[:,1])
            gAB23 = float(gB[1,:] * gATM[:,2])
            gAB31 = float(gB[2,:] * gATM[:,0])
            gAB32 = float(gB[2,:] * gATM[:,1])
            gAB33 = float(gB[2,:] * gATM[:,2])
            
#        gAB = np.matrix([[gAB11, gAB12, gAB13], [gAB21, gAB22, gAB23], [gAB31, gAB32, gAB33]])
            if gAB11 <= -1 and gAB11 >= -1.001:
                gAB11 = -1
            if gAB22 <= -1 and gAB22 >= -1.001:
                gAB22 = -1
            if gAB33 <= -1 and gAB33 >= -1.001:
                gAB33 = -1
            if gAB11 >= 1 and gAB11 <= 1.001:
                gAB11 = 1
            if gAB22 >= 1 and gAB22 <= 1.001:
                gAB22 = 1
            if gAB33 >= 1 and gAB33 <= 1.001:
                gAB33 = 1
            yyy = (gAB11 + gAB22 + gAB33 -1) * 0.5
            if yyy >= 1 and yyy < 1.001:
                yyy = 1
            if yyy <= -1 and yyy > -1.001:
                yyy= -1
            AngleRad = math.acos(yyy)
            AngleM[i, j] = AngleRad * 180 / math.pi
            
            j += 1
    
        if i == 0:  #to calculate c_axis Mis
            gAT = gATM
            gBT = gBTM #np.matrix([[gBT11, gBT12, gBT13],  [gBT21, gBT22, gBT23], [gBT31, gBT32, gBT33]])
        i += 1
    
    #--------calculate c-axis Misorientation (C_angle)
    uB1 = 0  
    vB1 = 0
    tB1 = 0
    wB1 = 1
    
    m1B1 = 1.5 * uB1
    m2B1 = np.sqrt(3)/2 *(2 * vB1 + uB1)
    m3B1 = wB1 * 1.62
    m4B1 = np.sqrt(m1B1 * m1B1 + m2B1 * m2B1 + m3B1 * m3B1)
    
    m1norB1 = m1B1/m4B1
    m2norB1 = m2B1/m4B1
    m3norB1 = m3B1/m4B1
    # burgers vectors
    Burg01 = gAT[0, 0] * m1norB1 + gAT[0, 1] * m2norB1 + gAT[0, 2] * m3norB1
    Burg02 = gAT[1, 0] * m1norB1 + gAT[1, 1] * m2norB1 + gAT[1, 2] * m3norB1
    Burg03 = gAT[2, 0] * m1norB1 + gAT[2, 1] * m2norB1 + gAT[2, 2] * m3norB1
    
    Burg11 = gBT[0, 0] * m1norB1 + gBT[0, 1] * m2norB1 + gBT[0, 2] * m3norB1
    Burg12 = gBT[1, 0] * m1norB1 + gBT[1, 1] * m2norB1 + gBT[1, 2] * m3norB1
    Burg13 = gBT[2, 0] * m1norB1 + gBT[2, 1] * m2norB1 + gBT[2, 2] * m3norB1
    
    dot_BurgAB = Burg01 * Burg11 + Burg02 * Burg12 + Burg03 * Burg13
    abs_BurgAB = math.sqrt(Burg01 ** 2 + Burg02 ** 2 + Burg03 ** 2) * math.sqrt(Burg11 ** 2 + Burg12 ** 2 + Burg13 ** 2)
    xxx = dot_BurgAB/abs_BurgAB
    if xxx >= 1 and xxx < 1.0001:
        xxx = 1
    if xxx <= -1 and xxx > 1.0001:
        xxx = -1
    
    c_AngleRad = math.acos(xxx)
    c_Angle = c_AngleRad * 180 / math.pi
    if c_Angle > 90:
        c_Angle = 180 - c_Angle
    
    Min_Angle = float(min(min(AngleM[0, :]), min(AngleM[1, :]), min(AngleM[2, :]), min(AngleM[3, :]), min(AngleM[4, :]), min(AngleM[5, :]), min(AngleM[6, :]), min(AngleM[7, :]), min(AngleM[8, :]), min(AngleM[9, :]), min(AngleM[10, :]), min(AngleM[11, :])))
#    if Min_Angle > 90:
#        Min_Angle = 180 - Min_Angle
    return(Min_Angle, c_Angle)


data_file0 = pd.read_excel(r'E:\办公室diannao20180802\PythonSpyder\Wangjie\Find Grain\grainfile_0%_asone.xls') #(1)the first file
data_file1 = pd.read_excel(r'E:\办公室diannao20180802\PythonSpyder\Wangjie\Find Grain\grainfile_1%.xls') #(2)the second file

Rawdata0 = data_file0.ix[:,'Number0':'A9']
Rawdata1 = data_file1.ix[:,'Number1':'Diameter1']

Rawdata0Matrix = np.matrix(Rawdata0)
Rawdata1Matrix = np.matrix(Rawdata1)

iii = np.shape(Rawdata0Matrix)[0]-1 #(3)number of grains for the first file -1
jjj = np.shape(Rawdata1Matrix)[0]-2 #(4)number of grains for the second file -2

i = 0
c = 0

p0x = np.max(Rawdata0Matrix[:,4])
p0y = np.max(Rawdata0Matrix[:,5])
p1x = np.max(Rawdata1Matrix[:,4])
p1y = np.max(Rawdata1Matrix[:,5])

NewMatrix = np.zeros((np.shape(Rawdata0Matrix)[0], 18))


cnt = 0
for i in tqdm(range(iii)):
    dis_angle = 10
    c_angle = 10

#    j = 0
#    a = rot(Rawdata0Matrix[i, 1], Rawdata0Matrix[i, 2], Rawdata0Matrix[i, 3], Rawdata1Matrix[i, 1], Rawdata1Matrix[i, 2], Rawdata1Matrix[i, 3])
    b = 0 #create a parameter
    test = 0
    for j in range(jjj):
        if (abs(float(Rawdata0Matrix[i, 0]) - float(Rawdata1Matrix[j, 0]) * (iii / jjj))) < 100:
            if math.sqrt(((float(Rawdata1Matrix[j, 4]) * (p0x / p1x)) - float(Rawdata0Matrix[i, 4]))**2 + ((float(Rawdata1Matrix[j, 5]) * (p0y / p1y)) - float(Rawdata0Matrix[i, 5]))**2) < 50:
                if abs(float(Rawdata0Matrix[i, 7]) / float(Rawdata1Matrix[j, 7])) < 5 and abs(float(Rawdata0Matrix[i, 7]) / float(Rawdata1Matrix[j, 7])) > 0.2:
                    dis, c_mis = rot(Rawdata0Matrix[i, 1], Rawdata0Matrix[i, 2], Rawdata0Matrix[i, 3], Rawdata1Matrix[j, 1], Rawdata1Matrix[j, 2], Rawdata1Matrix[j, 3])
                    if dis < dis_angle:
        #                        a = dis
        #                        b = j
        #                        c = i
                        test = test + 1
                        for k in range(8):
                            NewMatrix[cnt, k] = Rawdata0Matrix[i, k]
                            NewMatrix[cnt, k+8] = Rawdata1Matrix[j, k]
                        NewMatrix[cnt, 16] = dis
                        NewMatrix[cnt, 17] = c_mis
                        if test == 1:
                            cnt += 1
                            print('cnt'+str(cnt))
                        dis_angle = dis
                        c_angle = c_mis                
    print(c, b) #查看进度
#    Rawdata0Matrix[i, 8] = float(Rawdata1Matrix[b,0])
#    Rawdata0Matrix[i, 9] = float(Rawdata1Matrix[b,1])
#    Rawdata0Matrix[i, 10] = float(Rawdata1Matrix[b,2])
#    Rawdata0Matrix[i, 11] = float(Rawdata1Matrix[b,3])
#    Rawdata0Matrix[i, 12] = float(Rawdata1Matrix[b,4])
#    Rawdata0Matrix[i, 13] = float(Rawdata1Matrix[b,5])
#    Rawdata0Matrix[i, 14] = float(Rawdata1Matrix[b,6])
#    Rawdata0Matrix[i, 15] = float(Rawdata1Matrix[b,7])
#    Rawdata0Matrix[i, 16] = a       
#    i += 1

OutputExcel = pd.DataFrame(NewMatrix)
OutputExcel.to_csv('OutputFindgrains(rotation_axis_0pct_1pct).csv')  #(5)


























data_file0 = pd.read_excel(r'E:\办公室diannao20180802\PythonSpyder\Wangjie\Find Grain\grainfile_1%_asone.xls') #(1)the first file
data_file1 = pd.read_excel(r'E:\办公室diannao20180802\PythonSpyder\Wangjie\Find Grain\grainfile_2%.xls') #(2)the second file

Rawdata0 = data_file0.ix[:,'Number0':'A9']
Rawdata1 = data_file1.ix[:,'Number1':'Diameter1']

Rawdata0Matrix = np.matrix(Rawdata0)
Rawdata1Matrix = np.matrix(Rawdata1)

iii = np.shape(Rawdata0Matrix)[0]-1 #(3)number of grains for the first file -1
jjj = np.shape(Rawdata1Matrix)[0]-2 #(4)number of grains for the second file -2

i = 0
c = 0

p0x = np.max(Rawdata0Matrix[:,4])
p0y = np.max(Rawdata0Matrix[:,5])
p1x = np.max(Rawdata1Matrix[:,4])
p1y = np.max(Rawdata1Matrix[:,5])

NewMatrix = np.zeros((np.shape(Rawdata0Matrix)[0], 18))


cnt = 0
for i in tqdm(range(iii)):
    dis_angle = 10
    c_angle = 10

#    j = 0
#    a = rot(Rawdata0Matrix[i, 1], Rawdata0Matrix[i, 2], Rawdata0Matrix[i, 3], Rawdata1Matrix[i, 1], Rawdata1Matrix[i, 2], Rawdata1Matrix[i, 3])
    b = 0 #create a parameter
    test = 0
    for j in range(jjj):
        if (abs(float(Rawdata0Matrix[i, 0]) - float(Rawdata1Matrix[j, 0]) * (iii / jjj))) < 100:
            if math.sqrt(((float(Rawdata1Matrix[j, 4]) * (p0x / p1x)) - float(Rawdata0Matrix[i, 4]))**2 + ((float(Rawdata1Matrix[j, 5]) * (p0y / p1y)) - float(Rawdata0Matrix[i, 5]))**2) < 50:
                if abs(float(Rawdata0Matrix[i, 7]) / float(Rawdata1Matrix[j, 7])) < 5 and abs(float(Rawdata0Matrix[i, 7]) / float(Rawdata1Matrix[j, 7])) > 0.2:
                    dis, c_mis = rot(Rawdata0Matrix[i, 1], Rawdata0Matrix[i, 2], Rawdata0Matrix[i, 3], Rawdata1Matrix[j, 1], Rawdata1Matrix[j, 2], Rawdata1Matrix[j, 3])
                    if dis < dis_angle:
        #                        a = dis
        #                        b = j
        #                        c = i
                        test = test + 1
                        for k in range(8):
                            NewMatrix[cnt, k] = Rawdata0Matrix[i, k]
                            NewMatrix[cnt, k+8] = Rawdata1Matrix[j, k]
                        NewMatrix[cnt, 16] = dis
                        NewMatrix[cnt, 17] = c_mis
                        if test == 1:
                            cnt += 1
                            print('cnt'+str(cnt))
                        dis_angle = dis
                        c_angle = c_mis                
    print(c, b) #查看进度

#    Rawdata0Matrix[i, 8] = float(Rawdata1Matrix[b,0])
#    Rawdata0Matrix[i, 9] = float(Rawdata1Matrix[b,1])
#    Rawdata0Matrix[i, 10] = float(Rawdata1Matrix[b,2])
#    Rawdata0Matrix[i, 11] = float(Rawdata1Matrix[b,3])
#    Rawdata0Matrix[i, 12] = float(Rawdata1Matrix[b,4])
#    Rawdata0Matrix[i, 13] = float(Rawdata1Matrix[b,5])
#    Rawdata0Matrix[i, 14] = float(Rawdata1Matrix[b,6])
#    Rawdata0Matrix[i, 15] = float(Rawdata1Matrix[b,7])
#    Rawdata0Matrix[i, 16] = a       
#    i += 1

OutputExcel = pd.DataFrame(NewMatrix)
OutputExcel.to_csv('OutputFindgrains(rotation_axis_1pct_2pct).csv')  #(5)

























data_file0 = pd.read_excel(r'E:\办公室diannao20180802\PythonSpyder\Wangjie\Find Grain\grainfile_2%_asone.xls') #(1)the first file
data_file1 = pd.read_excel(r'E:\办公室diannao20180802\PythonSpyder\Wangjie\Find Grain\grainfile_4%.xls') #(2)the second file

Rawdata0 = data_file0.ix[:,'Number0':'A9']
Rawdata1 = data_file1.ix[:,'Number1':'Diameter1']

Rawdata0Matrix = np.matrix(Rawdata0)
Rawdata1Matrix = np.matrix(Rawdata1)

iii = np.shape(Rawdata0Matrix)[0]-1 #(3)number of grains for the first file -1
jjj = np.shape(Rawdata1Matrix)[0]-2 #(4)number of grains for the second file -2

i = 0
c = 0

p0x = np.max(Rawdata0Matrix[:,4])
p0y = np.max(Rawdata0Matrix[:,5])
p1x = np.max(Rawdata1Matrix[:,4])
p1y = np.max(Rawdata1Matrix[:,5])

NewMatrix = np.zeros((np.shape(Rawdata0Matrix)[0], 18))


cnt = 0
for i in tqdm(range(iii)):
    dis_angle = 10
    c_angle = 10

#    j = 0
#    a = rot(Rawdata0Matrix[i, 1], Rawdata0Matrix[i, 2], Rawdata0Matrix[i, 3], Rawdata1Matrix[i, 1], Rawdata1Matrix[i, 2], Rawdata1Matrix[i, 3])
    b = 0 #create a parameter
    test = 0
    for j in range(jjj):
        if (abs(float(Rawdata0Matrix[i, 0]) - float(Rawdata1Matrix[j, 0]) * (iii / jjj))) < 100:
            if math.sqrt(((float(Rawdata1Matrix[j, 4]) * (p0x / p1x)) - float(Rawdata0Matrix[i, 4]))**2 + ((float(Rawdata1Matrix[j, 5]) * (p0y / p1y)) - float(Rawdata0Matrix[i, 5]))**2) < 50:
                if abs(float(Rawdata0Matrix[i, 7]) / float(Rawdata1Matrix[j, 7])) < 5 and abs(float(Rawdata0Matrix[i, 7]) / float(Rawdata1Matrix[j, 7])) > 0.2:
                    dis, c_mis = rot(Rawdata0Matrix[i, 1], Rawdata0Matrix[i, 2], Rawdata0Matrix[i, 3], Rawdata1Matrix[j, 1], Rawdata1Matrix[j, 2], Rawdata1Matrix[j, 3])
                    if dis < dis_angle:
        #                        a = dis
        #                        b = j
        #                        c = i
                        test = test + 1
                        for k in range(8):
                            NewMatrix[cnt, k] = Rawdata0Matrix[i, k]
                            NewMatrix[cnt, k+8] = Rawdata1Matrix[j, k]
                        NewMatrix[cnt, 16] = dis
                        NewMatrix[cnt, 17] = c_mis
                        if test == 1:
                            cnt += 1
                            print('cnt'+str(cnt))
                        dis_angle = dis
                        c_angle = c_mis                
    print(c, b) #查看进度

#    Rawdata0Matrix[i, 8] = float(Rawdata1Matrix[b,0])
#    Rawdata0Matrix[i, 9] = float(Rawdata1Matrix[b,1])
#    Rawdata0Matrix[i, 10] = float(Rawdata1Matrix[b,2])
#    Rawdata0Matrix[i, 11] = float(Rawdata1Matrix[b,3])
#    Rawdata0Matrix[i, 12] = float(Rawdata1Matrix[b,4])
#    Rawdata0Matrix[i, 13] = float(Rawdata1Matrix[b,5])
#    Rawdata0Matrix[i, 14] = float(Rawdata1Matrix[b,6])
#    Rawdata0Matrix[i, 15] = float(Rawdata1Matrix[b,7])
#    Rawdata0Matrix[i, 16] = a       
#    i += 1

OutputExcel = pd.DataFrame(NewMatrix)
OutputExcel.to_csv('OutputFindgrains(rotation_axis_2pct_4pct).csv')  #(5)










data_file0 = pd.read_excel(r'E:\办公室diannao20180802\PythonSpyder\Wangjie\Find Grain\grainfile_4%_asone.xls') #(1)the first file
data_file1 = pd.read_excel(r'E:\办公室diannao20180802\PythonSpyder\Wangjie\Find Grain\grainfile_8%.xls') #(2)the second file

Rawdata0 = data_file0.ix[:,'Number0':'A9']
Rawdata1 = data_file1.ix[:,'Number1':'Diameter1']

Rawdata0Matrix = np.matrix(Rawdata0)
Rawdata1Matrix = np.matrix(Rawdata1)

iii = np.shape(Rawdata0Matrix)[0]-1 #(3)number of grains for the first file -1
jjj = np.shape(Rawdata1Matrix)[0]-2 #(4)number of grains for the second file -2

i = 0
c = 0

p0x = np.max(Rawdata0Matrix[:,4])
p0y = np.max(Rawdata0Matrix[:,5])
p1x = np.max(Rawdata1Matrix[:,4])
p1y = np.max(Rawdata1Matrix[:,5])

NewMatrix = np.zeros((np.shape(Rawdata0Matrix)[0], 18))


cnt = 0
for i in tqdm(range(iii)):
    dis_angle = 10
    c_angle = 10

#    j = 0
#    a = rot(Rawdata0Matrix[i, 1], Rawdata0Matrix[i, 2], Rawdata0Matrix[i, 3], Rawdata1Matrix[i, 1], Rawdata1Matrix[i, 2], Rawdata1Matrix[i, 3])
    b = 0 #create a parameter
    test = 0
    for j in range(jjj):
        if (abs(float(Rawdata0Matrix[i, 0]) - float(Rawdata1Matrix[j, 0]) * (iii / jjj))) < 100000:
            if math.sqrt(((float(Rawdata1Matrix[j, 4]) * (p0x / p1x)) - float(Rawdata0Matrix[i, 4]))**2 + ((float(Rawdata1Matrix[j, 5]) * (p0y / p1y)) - float(Rawdata0Matrix[i, 5]))**2) < 50:
                if abs(float(Rawdata0Matrix[i, 7]) / float(Rawdata1Matrix[j, 7])) < 5 and abs(float(Rawdata0Matrix[i, 7]) / float(Rawdata1Matrix[j, 7])) > 0.2:
                    dis, c_mis = rot(Rawdata0Matrix[i, 1], Rawdata0Matrix[i, 2], Rawdata0Matrix[i, 3], Rawdata1Matrix[j, 1], Rawdata1Matrix[j, 2], Rawdata1Matrix[j, 3])
                    if dis < dis_angle:
        #                        a = dis
        #                        b = j
        #                        c = i
                        test = test + 1
                        for k in range(8):
                            NewMatrix[cnt, k] = Rawdata0Matrix[i, k]
                            NewMatrix[cnt, k+8] = Rawdata1Matrix[j, k]
                        NewMatrix[cnt, 16] = dis
                        NewMatrix[cnt, 17] = c_mis
                        if test == 1:
                            cnt += 1
                            print('cnt'+str(cnt))
                        dis_angle = dis
                        c_angle = c_mis                
    print(c, b) #查看进度
#    Rawdata0Matrix[i, 8] = float(Rawdata1Matrix[b,0])
#    Rawdata0Matrix[i, 9] = float(Rawdata1Matrix[b,1])
#    Rawdata0Matrix[i, 10] = float(Rawdata1Matrix[b,2])
#    Rawdata0Matrix[i, 11] = float(Rawdata1Matrix[b,3])
#    Rawdata0Matrix[i, 12] = float(Rawdata1Matrix[b,4])
#    Rawdata0Matrix[i, 13] = float(Rawdata1Matrix[b,5])
#    Rawdata0Matrix[i, 14] = float(Rawdata1Matrix[b,6])
#    Rawdata0Matrix[i, 15] = float(Rawdata1Matrix[b,7])
#    Rawdata0Matrix[i, 16] = a       
#    i += 1

OutputExcel = pd.DataFrame(NewMatrix)
OutputExcel.to_csv('OutputFindgrains(rotation_axis_4pct_8pct).csv')  #(5)























data_file0 = pd.read_excel(r'E:\办公室diannao20180802\PythonSpyder\Wangjie\Find Grain\grainfile_8%_asone.xls') #(1)the first file
data_file1 = pd.read_excel(r'E:\办公室diannao20180802\PythonSpyder\Wangjie\Find Grain\grainfile_16%.xls') #(2)the second file

Rawdata0 = data_file0.ix[:,'Number0':'A9']
Rawdata1 = data_file1.ix[:,'Number1':'Diameter1']

Rawdata0Matrix = np.matrix(Rawdata0)
Rawdata1Matrix = np.matrix(Rawdata1)

iii = np.shape(Rawdata0Matrix)[0]-1 #(3)number of grains for the first file -1
jjj = np.shape(Rawdata1Matrix)[0]-2 #(4)number of grains for the second file -2

i = 0
c = 0

p0x = np.max(Rawdata0Matrix[:,4])
p0y = np.max(Rawdata0Matrix[:,5])
p1x = np.max(Rawdata1Matrix[:,4])
p1y = np.max(Rawdata1Matrix[:,5])

NewMatrix = np.zeros((np.shape(Rawdata0Matrix)[0], 18))


cnt = 0
for i in tqdm(range(iii)):
    dis_angle = 10
    c_angle = 10

#    j = 0
#    a = rot(Rawdata0Matrix[i, 1], Rawdata0Matrix[i, 2], Rawdata0Matrix[i, 3], Rawdata1Matrix[i, 1], Rawdata1Matrix[i, 2], Rawdata1Matrix[i, 3])
    b = 0 #create a parameter
    test = 0
    for j in range(jjj):
        if (abs(float(Rawdata0Matrix[i, 0]) - float(Rawdata1Matrix[j, 0]) * (iii / jjj))) < 100000:
            if math.sqrt(((float(Rawdata1Matrix[j, 4]) * (p0x / p1x)) - float(Rawdata0Matrix[i, 4]))**2 + ((float(Rawdata1Matrix[j, 5]) * (p0y / p1y)) - float(Rawdata0Matrix[i, 5]))**2) < 50:
                if abs(float(Rawdata0Matrix[i, 7]) / float(Rawdata1Matrix[j, 7])) < 5 and abs(float(Rawdata0Matrix[i, 7]) / float(Rawdata1Matrix[j, 7])) > 0.2:
                    dis, c_mis = rot(Rawdata0Matrix[i, 1], Rawdata0Matrix[i, 2], Rawdata0Matrix[i, 3], Rawdata1Matrix[j, 1], Rawdata1Matrix[j, 2], Rawdata1Matrix[j, 3])
                    if dis < dis_angle:
        #                        a = dis
        #                        b = j
        #                        c = i
                        test = test + 1
                        for k in range(8):
                            NewMatrix[cnt, k] = Rawdata0Matrix[i, k]
                            NewMatrix[cnt, k+8] = Rawdata1Matrix[j, k]
                        NewMatrix[cnt, 16] = dis
                        NewMatrix[cnt, 17] = c_mis
                        if test == 1:
                            cnt += 1
                            print('cnt'+str(cnt))
                        dis_angle = dis
                        c_angle = c_mis                
    print(c, b) #查看进度

#    Rawdata0Matrix[i, 8] = float(Rawdata1Matrix[b,0])
#    Rawdata0Matrix[i, 9] = float(Rawdata1Matrix[b,1])
#    Rawdata0Matrix[i, 10] = float(Rawdata1Matrix[b,2])
#    Rawdata0Matrix[i, 11] = float(Rawdata1Matrix[b,3])
#    Rawdata0Matrix[i, 12] = float(Rawdata1Matrix[b,4])
#    Rawdata0Matrix[i, 13] = float(Rawdata1Matrix[b,5])
#    Rawdata0Matrix[i, 14] = float(Rawdata1Matrix[b,6])
#    Rawdata0Matrix[i, 15] = float(Rawdata1Matrix[b,7])
#    Rawdata0Matrix[i, 16] = a       
#    i += 1

OutputExcel = pd.DataFrame(NewMatrix)
OutputExcel.to_csv('OutputFindgrains(rotation_axis_8pct_16pct).csv')  #(5)







# Copyright (C) 2021 Edgar Sutawika - All Rights Reserve
# For educational purposes

# frequently used calculation functions
import numpy as np
import pandas as pd

def prettyMatVect2(matVect):
    prettyMatVect = pd.DataFrame(matVect, columns =
                                ['Rx', 'Ry', 'Rz', 't0', 't1',
                                't2', 't3', 'c1'])
    return prettyMatVect

def ATrans (theta0, theta1, theta2, theta3): #A
    theta0 = float(theta0)
    theta1 = float(theta1)
    theta2 = float(theta2)
    theta3 = float(theta3)

    a11 = 1-2*theta2**2-2*theta3**2
    a12 = 2*(theta1*theta2-theta0*theta3)
    a13 = 2*(theta1*theta3+theta0*theta2)

    a21 = 2*(theta1*theta2+theta0*theta3)
    a22 = 1-2*theta1**2-2*theta3**2
    a23 = 2*(theta3*theta2-theta0*theta1)

    a31 = 2*(theta1*theta3-theta0*theta2)
    a32 = 2*(theta3*theta2+theta0*theta1)
    a33 = 1-2*theta1**2-2*theta2**2

    Atrans = np.array([[a11, a12, a13],
                       [a21, a22, a23],
                       [a31, a32, a33]], dtype = float)
    return Atrans

def uBarSkew(ubiP):
    ubiPx = float(ubiP[0])
    ubiPy = float(ubiP[1])
    ubiPz = float(ubiP[2])

    ubarSkew = np.array([[      0, -ubiPz,  ubiPy],
                         [  ubiPz,      0, -ubiPx],
                         [ -ubiPy,  ubiPx,      0]], dtype = float)

    return ubarSkew

def GBarMat(link, qi):

    t0 = float(qi[link2index(link, "t0")])
    t1 = float(qi[link2index(link, "t1")])
    t2 = float(qi[link2index(link, "t2")])
    t3 = float(qi[link2index(link, "t3")])

    EBarMat = np.array([ [-t1,  t0,  t3, -t2],
                         [-t2, -t3,  t0,  t1],
                         [-t3,  t2, -t1,  t0] ], dtype = float)
          
    GbarMat = 2*EBarMat
    return GbarMat

def GBarMatDot(link, qiDot):
    t0Dot = float(qiDot[link2index(link, "t0")])
    t1Dot = float(qiDot[link2index(link, "t1")])
    t2Dot = float(qiDot[link2index(link, "t2")])
    t3Dot = float(qiDot[link2index(link, "t3")])

    EBarMatDot = np.array([ [-t1Dot,  t0Dot,  t3Dot, -t2Dot],
                            [-t2Dot, -t3Dot,  t0Dot,  t1Dot],
                            [-t3Dot,  t2Dot, -t1Dot,  t0Dot] ], dtype = float)
          
    GbarMatDot = 2*EBarMatDot
    return GbarMatDot


def local2global3D(qi, u_bar_iP, link_i):
    # To calculate Point of Interest positions in terms of global coordinates
    id_RR    = link2index(link_i, "x")
    id_Theta = link2index(link_i, "t0")

    Ri = np.array([qi[  id_RR], 
                   qi[1+id_RR], 
                   qi[2+id_RR]], dtype = float) 

    t0, t1, t2, t3 = qi[id_Theta], qi[id_Theta+1], qi[id_Theta+2], qi[id_Theta+3]

    A_matrix = ATrans(t0, t1, t2, t3)
    riP      = Ri + A_matrix @ u_bar_iP
    
    return riP

def local2globalDot3D(qi, qiDot, u_bar_iP, link_i):
    # To calculate Point of Interest positions in terms of global coordinates
    id_RR = link2index(link_i, "x")
    id_Theta = link2index(link_i, "t0")
    
    Ri_Dot = np.array([qiDot[  id_RR], 
                       qiDot[1+id_RR], 
                       qiDot[2+id_RR]], dtype = float)

    t0, t1, t2, t3 = qi[id_Theta], qi[id_Theta+1], qi[id_Theta+2], qi[id_Theta+3]
    
    A_matrix = ATrans(t0, t1, t2, t3)
    ubarSkewMat = uBarSkew(u_bar_iP)
    GBar = GBarMat(1, qi)
    thetaDot = np.array([ qiDot[id_RR], 
                          qiDot[id_RR+1], 
                          qiDot[id_RR+2], 
                          qiDot[id_RR+3] ], dtype = float)

    riP_Dot = Ri_Dot - A_matrix @ ubarSkewMat @ GBar @ thetaDot

    return riP_Dot

def link2index(link, string):
    if string == "x":
        index = 7*(link-1)
    elif string == "y":
        index = 7*(link-1)+1
    elif string == "z":
        index = 7*(link-1)+2
    elif string == "t0":
        index = 7*(link-1)+3
    elif string == "t1":
        index = 7*(link-1)+4
    elif string == "t2":
        index = 7*(link-1)+5
    elif string == "t3":
        index = 7*(link-1)+6

    index = int(index)
    return index


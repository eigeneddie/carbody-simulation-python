# Copyright (C) 2021 Edgar Sutawika - All Rights Reserve
# For educational purposes

import numpy as np
from calcModule3D import link2index as l2i, uBarSkew as uBSkew

def QeThetaEP(GBarMat, ATrans, uBarSkew, F_R_i):
    uSkew = np.dot(ATrans, uBarSkew)
    matrixA = np.dot(uSkew, GBarMat)
    F_theta_i = -np.transpose(matrixA)@F_R_i

    return F_theta_i

def linSprg3D(GBarMat, uBariP, ATrans, stiffness, riP, rjP, lo):
    Ls       = riP-rjP
    LsMag    = np.sqrt(float(np.dot(np.transpose(Ls),Ls))) # Scalar
    
    FsMag    = stiffness*(LsMag-lo) # Scalar
    
    
    if LsMag == 0: #avoid division by zero
        LsUnit = np.zeros((3,1))
    else:
        LsUnit   = Ls/LsMag # unit vector 3x1

    uBarSkew = uBSkew(uBariP)

    Fs_i     = -FsMag*LsUnit # vector 3x1
    Fs_j     =  FsMag*LsUnit # vector 3x1
    QTheta_i = -np.transpose(ATrans@uBarSkew@GBarMat)@Fs_i # vector 4x1
    QTheta_j = -np.transpose(ATrans@uBarSkew@GBarMat)@Fs_j # vector 4x1

    return Fs_i, Fs_j, QTheta_i, QTheta_j

def linDamp3D(GBarMat, uBariP, ATrans, damping, riP_Dot, rjP_Dot):
    LsDot       = riP_Dot-rjP_Dot
    LsDotMag    = np.sqrt(float(np.dot(np.transpose(LsDot),LsDot)))
    FdMag       = damping*LsDotMag

    if LsDotMag ==0: #avoid division by zero
        LsDotUnit = np.zeros((3,1))
    else:
        LsDotUnit   = LsDot/LsDotMag# unit vector 3x1

    uBarSkew    = uBSkew(uBariP)
    '''print("rip")
    print(riP_Dot)
    print("rjp")
    print(rjP_Dot)'''
    Fd_i     = -FdMag*LsDotUnit
    Fd_j     =  FdMag*LsDotUnit
    QTheta_i =  -np.transpose(ATrans@uBarSkew@GBarMat)@Fd_i
    QTheta_j =  -np.transpose(ATrans@uBarSkew@GBarMat)@Fd_j

    return Fd_i, Fd_j, QTheta_i, QTheta_j


#=========================================
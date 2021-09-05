# Copyright (C) 2021 Edgar Sutawika - All Rights Reserve
# For educational purposes

import numpy as np
from calcModule3D import link2index as l2i, uBarSkew as uBSkew

def QeThetaEP(GBarMat, ATrans, uBarSkew, F_R_i): #random external force
    uSkew = np.dot(ATrans, uBarSkew)
    matrixA = np.dot(uSkew, GBarMat)
    F_theta_i = -np.transpose(matrixA)@F_R_i

    return F_theta_i

def linFS3D(GBarMat, uBariP, ATrans, stiffness, damping, riP, rjP, riPDot, rjPDot, lo):
    Ls          = riP-rjP
    LsMag       = np.linalg.norm(Ls)
    LsDot       = riPDot-rjPDot
    LsDotMag    = np.linalg.norm(LsDot)
    FsMag       = stiffness*(LsMag-lo) + damping*LsDotMag # Scalar
    
    if LsMag == 0: #avoid division by zero
        LsUnit = np.zeros((3,1))
    else:
        LsUnit   = Ls/LsMag # unit vector 3x1

    uBarSkew = uBSkew(uBariP)

    Fs_i     = -FsMag*LsUnit # vector 3x1
    Fs_j     =  FsMag*LsUnit # vector 3x1
    QTheta_i = -np.transpose(ATrans@uBarSkew@GBarMat)@Fs_i # vector 4x1
    QTheta_j = -np.transpose(ATrans@uBarSkew@GBarMat)@Fs_j # vector 4x1

    return Fs_i, Fs_j, QTheta_i, QTheta_j, LsDotMag

#=========================================
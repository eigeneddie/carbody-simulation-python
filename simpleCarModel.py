# Copyright (C) 2021 Edgar Sutawika - All Rights Reserve
# For educational purposes

# CAR BODY SIMULATION 

# python module
import numpy as np  
import matplotlib.pyplot  as plt

# custom module
from calcModule3D import GBarMat, GBarMatDot, link2index as l2i
from calcModule3D import local2global3D as l2g, prettyMatVect2
from calcModule3D import local2globalDot3D as l2gDot

import calcModule3D       as calMod
import forceModule3D      as fMod
import constraintModule3D as conMod


# 1. === USER INPUT PARAMETERS (GLOBAL VARIABLES) ======

# 1.1. simulation parameter
#---------------------------
timeStart, timeEnd, stepSize = 0, 60,0.05 # [s] 
time = np.arange(timeStart, timeEnd, stepSize)
epsilon = 0.000000000001
gravity = 9.81 # [m/s^2]
n, nc = 7, 1 # gen. coord., const. eq.
             # simulation only consist of 1 body

# 1.2. car model parameters
#---------------------------
mass_body = 1500 # [kg]

axleLength   = 1.8 # [m]
axleDistance = 2.8 # [m] distance betwen front and rear
bodyHeight   = 1.4 # [m]

staticSpringLength = 0.2 #[m] (lo)
springK            = 21000 #[N/m]
damperC            = 200 # [Ns/m]

# 1.3. Uneven road  parameters
#---------------------------
# assuming road unevenness is defined by a sine wave

velocityCar = 2.78      # [m/s], 10 km/hr Constant velocity
lamdaWave   = 15         # [m] wave length of uneven road
phaseLR     = np.pi/3   # [phase difference between left and right wheel]
roadAmp     = 0.05       # [m] uneven road amplitude

# 2. ===== DERIVED PARAMETERS =========

# 2.1. Total spring length
springLength = staticSpringLength + mass_body*gravity/4/springK

# 2.2 Mass matrix properties (MASS MODULE)

Ixx = 1/12*mass_body*(  axleLength**2+bodyHeight**2)
Iyy = 1/12*mass_body*(axleDistance**2+axleLength**2)
Izz = 1/12*mass_body*(axleDistance**2+bodyHeight**2)

massRR = np.array([ [mass_body,         0,         0],
                    [        0, mass_body,         0],
                    [        0,         0, mass_body] ])

Itheta2 = np.array([ [Ixx,   0,   0],
                     [  0, Iyy,   0],
                     [  0,   0, Izz] ])

# 2.3. Initial conditions of body
Ry1 =  staticSpringLength+bodyHeight/2  # [m] everything else coincides 
                                      #     with global coordinate
theta0Init = 1
#===============================

# 3. ====POINTS OF INTEREST, LOCAL JOINTS=====
# 3.1. car body point of interest (local) CONSTANT THROUGHOUT
uBar_FLW  = np.array([[ axleDistance/2], [-bodyHeight/2], [-axleLength/2] ]) 
uBar_FRW  = np.array([[ axleDistance/2], [-bodyHeight/2], [ axleLength/2] ])
uBar_RLW  = np.array([[-axleDistance/2], [-bodyHeight/2], [-axleLength/2] ])
uBar_RRW  = np.array([[-axleDistance/2], [-bodyHeight/2], [ axleLength/2] ])

# 3.2. ground point of interest (global) CONSTANT THROUGHOUT

# Position
rG_FLW = np.array([[ axleDistance/2], [0], [-axleLength/2] ]) 
rG_FRW = np.array([[ axleDistance/2], [0], [ axleLength/2] ])
rG_RLW = np.array([[-axleDistance/2], [0], [-axleLength/2] ])
rG_RRW = np.array([[-axleDistance/2], [0], [ axleLength/2] ])

# Velocity: initially ZERO
rG_FLW_Dot = np.zeros((3,1))
rG_FRW_Dot = np.zeros((3,1))
rG_RLW_Dot = np.zeros((3,1))
rG_RRW_Dot = np.zeros((3,1))

# 3.3. ====Memory reservation===============
#   a. Generalized coordinates and derivatives
qi             = np.zeros((n,1))    # gen. position
qiDot          = np.zeros((n,1))    # gen. velocity
qiDotDot_lamda = np.zeros((n+nc,1)) # gen. acceleration

qi[l2i(1, "y")]  = Ry1
qi[l2i(1, "t0")] = theta0Init

#    b. Save generalized coordinates for all time
q_allTime       = np.zeros((np.size(time), n))
v_allTime       = np.zeros((np.size(time), n))
a_allTime       = np.zeros((np.size(time), n))
checkvar1        = np.zeros((np.size(time), 3))
checkvar2        = np.zeros((np.size(time), 3))
# 4. ====Main Program====
def mainProg():
    print(" ")
    print(" ")
    print("========START SIMULATION=======")
    print(" ")

    global qi, qiDot, qiDotDot_lamda
    global rG_FRW, rG_RRW, rG_FLW, rG_RLW
    global rG_FLW_Dot, rG_FRW_Dot, rG_RLW_Dot, rG_RRW_Dot

    timeNow  = timeStart
    loopCount = 0

    #rGAllTime = np.zeros((np.size(time),1))
    #r1AllTime = np.zeros((np.size(time),1))

    for timeID in range(np.size(time)): # loop @ t

        # Activate road unevenness at t = timeEnd/3
        if timeNow > timeEnd/3:

            omega_timeFRW = 2*np.pi*velocityCar/lamdaWave*timeNow
            omega_timeRRW = 2*np.pi*velocityCar/lamdaWave*timeNow - axleDistance/lamdaWave*2*np.pi
            omega_timeFLW = 2*np.pi*velocityCar/lamdaWave*timeNow - phaseLR
            omega_timeRLW = 2*np.pi*velocityCar/lamdaWave*timeNow - axleDistance/lamdaWave*2*np.pi - phaseLR

            # ===Global position of rG_P

            # Right wheels
            rG_FRW[1] = roadAmp*np.sin(omega_timeFRW)
            rG_RRW[1] = roadAmp*np.sin(omega_timeRRW)

            # Left wheels
            rG_FLW[1] = roadAmp*np.sin(omega_timeFLW)
            rG_RLW[1] = roadAmp*np.sin(omega_timeRLW)
            

            # ===Global Velocity of rG_P

            # Right wheels
            rG_FRW_Dot[1] = roadAmp*2*np.pi*velocityCar/lamdaWave*np.cos(omega_timeFRW)
            rG_RRW_Dot[1] = roadAmp*2*np.pi*velocityCar/lamdaWave*np.cos(omega_timeRRW)

            # Left wheels
            rG_FLW_Dot[1] = roadAmp*2*np.pi*velocityCar/lamdaWave*np.cos(omega_timeFLW)
            rG_RLW_Dot[1] = roadAmp*2*np.pi*velocityCar/lamdaWave*np.cos(omega_timeRLW)
        
        Cq, _ = config(qi) # Cq @ t, ignore constraintVect

        # 5. ====FIND ACCELERATION @ t =====
        qiDotDot_lamda, checkValue, checkValue2 = systemEquation(timeNow, Cq, qi, qiDot)

        # 6. ====STORE EVERYTHING @ t ========
        q_allTime[timeID,:]      = qi.T
        v_allTime[timeID,:]      = qiDot.T
        a_allTime[timeID,:]      = qiDotDot_lamda[0:n].T
        checkvar1 [timeID]        = checkValue.T
        checkvar2 [timeID]        = checkValue2.T
        print("Loop Count:")
        print(loopCount)
        print(" ")

        print("Generalized Coordinate")
        print(qi.T)
        print(" ")

        print("Velocity")
        print(qiDot.T)
        print(" ")  
        '''
        print("Acceleration:")
        print(qiDotDot_lamda[0:n].T)
        print(" ")'''

        print("Check value:")
        #print(prettyMatVect2(massAug))
        print(checkValue)
        print(" ")

        # 7.====CALCULATE q, qdot  @ t+1 ======
        qi, qiDot = rungeKutta4_AtTimeNow( qi, qiDot, systemEquation, 
                                            stepSize, timeNow)

        loopCount = loopCount +1
        timeNow = timeNow + stepSize

    #print(np.size(time))
    #print(q_allTime[:, l2i(1, "y")])
    
    # 8. =====PLOT, BABY, PLOT!! ======
    plt.figure(1)
    plt.plot(time, q_allTime[:, l2i(1, "x")])
    plt.title('y')
    plt.ylabel('position y')
    plt.xlabel('time [s]')
    plt.xlim((0,timeEnd))
    plt.grid(True)

    plt.figure(2)
    plt.plot(time, checkvar1[:,1])
    plt.plot(time, checkvar2[:,1])
    plt.title('y')
    plt.ylabel('position y')
    plt.xlabel('time [s]')
    plt.xlim((0,timeEnd))
    plt.grid(True)
    plt.legend(["r1FLW", "rGFLW"])

    plt.show()

# =======IMPORTANT CALCULATION FUNCTIONS
def config(qi): #OKAY!

    constraintVect = conMod.constraintEquation(qi, 1)
    Cq = conMod.jacobianMatrix(qi, 1)

    return Cq, constraintVect

def systemEquation(t, Cq, qi, qiDot):
    #===== I. CONSTRUCT MCq matrix (MASS MODULE)====

    massAugmented = np.zeros((n+nc, n+nc))

    th0 = qi[l2i(1, "t0")]
    th1 = qi[l2i(1, "t1")]
    th2 = qi[l2i(1, "t2")]
    th3 = qi[l2i(1, "t3")]
    
    th0Dot = qiDot[l2i(1, "t0")]
    th1Dot = qiDot[l2i(1, "t1")]
    th2Dot = qiDot[l2i(1, "t2")]
    th3Dot = qiDot[l2i(1, "t3")]

    GBar = GBarMat(1, qi)
    A_Matrix = calMod.ATrans(th0, th1, th2, th3)
    
    massAugmented[0:3,0:3] = massRR
    massAugmented[3:n,3:n] = np.transpose(GBar)@Itheta2@GBar 
    massAugmented[n:n+nc,0:n] = Cq
    massAugmented[0:n,n:n+nc] = np.transpose(Cq)

    #=====II. CONSTRUCT Qe & Qd vector ===
    Qe = np.zeros((n,1))
    
    # 1. External Force from Weight
    QeY1 = -mass_body*gravity

    # 2. External Force from spring

    # Position POI in Global Coordinate
    r1_FLW = l2g(qi, uBar_FLW, 1)
    r1_FRW = l2g(qi, uBar_FRW, 1)
    r1_RLW = l2g(qi, uBar_RLW, 1)
    r1_RRW = l2g(qi, uBar_RRW, 1)

    # Velocity POI in GLobal Coordinate
    r1_FLW_Dot = l2gDot(qi, qiDot, uBar_FLW, 1)
    r1_FRW_Dot = l2gDot(qi, qiDot, uBar_FRW, 1)
    r1_RLW_Dot = l2gDot(qi, qiDot, uBar_RLW, 1)
    r1_RRW_Dot = l2gDot(qi, qiDot, uBar_RRW, 1)

    # -FLW
    Fs1_FLW, FsG_FLW, QsTheta1_FLW, _ = fMod.linSprg3D(GBar, uBar_FLW, A_Matrix, 
                                            springK, r1_FLW, rG_FLW, springLength)
    
    # -FRW
    Fs1_FRW, _, QsTheta1_FRW, _ = fMod.linSprg3D(GBar, uBar_FRW, A_Matrix, 
                                            springK, r1_FRW, rG_FRW, springLength)
    # -RLW
    Fs1_RLW, _, QsTheta1_RLW, _ = fMod.linSprg3D(GBar, uBar_RLW, A_Matrix, 
                                            springK, r1_RLW, rG_RLW, springLength)
    # -RRW
    Fs1_RRW, _, QsTheta1_RRW, _ = fMod.linSprg3D(GBar, uBar_RRW, A_Matrix, 
                                            springK, r1_RRW, rG_RRW, springLength)
    
    springForceTotal  = Fs1_FLW + Fs1_FRW + Fs1_RLW + Fs1_RRW
    springMomentTotal = QsTheta1_FLW + QsTheta1_FRW + QsTheta1_RLW + QsTheta1_RRW

    # 3. External Force from damper
    # -FLW
    Fd1_FLW, _, QdTheta1_FLW, _ = fMod.linDamp3D(GBar, uBar_FLW, A_Matrix, 
                                            damperC, r1_FLW_Dot, rG_FLW_Dot)
    # -FRW
    Fd1_FRW, _, QdTheta1_FRW, _= fMod.linDamp3D(GBar, uBar_FRW, A_Matrix, 
                                            damperC, r1_FRW_Dot, rG_FRW_Dot)
    # -RLW
    Fd1_RLW, _, QdTheta1_RLW, _ = fMod.linDamp3D(GBar, uBar_RLW, A_Matrix, 
                                            damperC, r1_RLW_Dot, rG_RLW_Dot)
    # -RRW
    Fd1_RRW, _, QdTheta1_RRW, _ = fMod.linDamp3D(GBar, uBar_RRW, A_Matrix, 
                                            springK, r1_RRW_Dot, rG_RRW_Dot)
    
    dampingForceTotal  = 0#Fd1_FLW + Fd1_FRW + Fd1_RLW + Fd1_RRW
    dampingMomentTotal = 0#QdTheta1_FLW + QdTheta1_FRW + QdTheta1_RLW + QdTheta1_RRW

    Qe[l2i(1,"x") :l2i(1,"z")+1] = springForceTotal + dampingForceTotal
    Qe[l2i(1, "y")] = Qe[l2i(1, "y")] + QeY1
    
    # 4. Centrifugal forces
    thetaDot =  np.array([th0Dot, th1Dot, th2Dot, th3Dot])
    omegaBar =  GBar@thetaDot
    GBarDot  =  GBarMatDot(1, qiDot)
    IthetaOmega = np.transpose(Itheta2@omegaBar)
    dummyMat =  np.transpose(np.cross(np.transpose(omegaBar), IthetaOmega)) + Itheta2@GBarDot@thetaDot
    Qv1Theta = -np.matmul(np.transpose(GBar),dummyMat)

    Qe[l2i(1,"t0"):l2i(1,"t3")+1]= springMomentTotal  + Qv1Theta + dampingMomentTotal

    # 5. Qd vector
    Qd = np.array([[conMod.QdEulPar1(qiDot, 1)]], dtype = object)

    # 6. Qe & Qd vector
    QeAug = np.concatenate((Qe, Qd), axis = 0) # vector 8x1
    #print(QeAug)

    #=====III. SOLVE QIDOTDOT_LAMDA ===
    mass_MatInverse = np.linalg.inv(massAugmented)
    qiDotDot_lamda  = np.dot(mass_MatInverse, QeAug)

    return qiDotDot_lamda, r1_FLW_Dot, rG_FLW_Dot

def rungeKutta4_AtTimeNow(qi, qiDot, systemFunction, stepSize, timeNow):
    # This function works with ANY number of DOF
    x   = np.array([qi[l2i(1,  "x")], 
                    qi[l2i(1,  "y")],
                    qi[l2i(1,  "z")],
                    qi[l2i(1, "t0")],
                    qi[l2i(1, "t1")],
                    qi[l2i(1, "t2")],
                    qi[l2i(1, "t3")]])

    xDot = np.array([qiDot[l2i(1,  "x")], 
                     qiDot[l2i(1,  "y")],
                     qiDot[l2i(1,  "z")],
                     qiDot[l2i(1, "t0")],
                     qiDot[l2i(1, "t1")],
                     qiDot[l2i(1, "t2")],
                     qiDot[l2i(1, "t3")]])

    y = np.concatenate((x, xDot), axis = 0)
    numberOfDOF = int(np.size(y)/2)
    
    # RungeKutta4
    t1  = timeNow
    Cq,_ = config(qi)
    f_1,_,_ = systemFunction(t1, Cq, qi, qiDot)
    k1  = np.zeros((np.size(y), 1))        
    for x    in range(numberOfDOF):
        k1[x]             = y[x+numberOfDOF]
        k1[x+numberOfDOF] = f_1[x]
    
    t2  = t1+ 0.5*stepSize
    y2  = y + 0.5*k1*stepSize
    for i in range(numberOfDOF):
        qi   [i] = y2[i]
        qiDot[i] = y2[i+numberOfDOF]
    Cq,_ = config(qi)
    f_2,_,_ = systemFunction(t2, Cq, qi, qiDot)
    k2  = np.zeros((np.size(y), 1))           
    for x    in range(numberOfDOF):
        k2[x]             = y2[x+numberOfDOF]
        k2[x+numberOfDOF] = f_2[x]
    
    t3  = t1+0.5*stepSize
    y3  = y + 0.5*k2*stepSize
    for i in range(numberOfDOF):
        qi   [i] = y3[i]
        qiDot[i] = y3[i+numberOfDOF]
    Cq,_  = config(qi)
    f_3,_,_ = systemFunction(t3, Cq, qi, qiDot)
    k3  = np.zeros((np.size(y), 1))           
    for x    in range(numberOfDOF):
        k3[x]             = y3[x+numberOfDOF]
        k3[x+numberOfDOF] = f_3[x]
    
    t4  = t1+stepSize
    y4  = y + k3*stepSize
    for i in range(numberOfDOF):
        qi   [i] = y4[i]
        qiDot[i] = y4[i+numberOfDOF]
    Cq,_  = config(qi)
    f_4,_,_ = systemFunction(t4, Cq, qi, qiDot)
    k4  = np.zeros((np.size(y), 1))
    for x    in range(numberOfDOF):
        k4[x]             = y4[x+numberOfDOF]
        k4[x+numberOfDOF] = f_4[x]

    RKFunct = (k1 + 2*k2 + 2*k3 + k4)/6

    yNew = y + stepSize*RKFunct

    for i in range(numberOfDOF):
        qi   [i] = yNew[i]
        qiDot[i] = yNew[i+numberOfDOF]

    return qi, qiDot

# Run main program
if __name__=="__main__":
    mainProg()

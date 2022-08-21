#!/usr/bin/env python2
# -*- coding: utf-8 -*-
from __future__ import division
from threading import Thread
import numpy as np
import math
import matplotlib.pyplot as plt


def SCurvePara(Tf, v, a):
    T = np.zeros(7)
    for i in range(0,100):
        J = (a**2 * v) / (Tf*v*a - v**2 - a)
        T[0] = a / J
        T[1] = v / a - a / J #t2 = v / a - t1
        T[2] = T[0]
        T[3] = Tf - 2 * a / J - 2 * v / a    # t4 = Tf - 4*t1 - 2*t2
        T[4] = T[2]
        T[5] = T[1]
        T[6] = T[0]
        if T[1] < -1e-6:
            a = math.sqrt(v*J)
        elif T[3] < -1e-6:
            v = Tf*a/2 - a*a/J
        elif J < -1e-6:
            Tf = (v**2 + a) / (v*a) + 1e-1
        else:
            break
    A = a
    V = v
    TF = Tf
    return TF,V,A,J,T

def SCurveScaling(t,V,A,J,T,Tf):
    if A <= V*V:
        s = 0
        print("input error")
        return s
    else:
        if A <= 2*V*V:
            if (Tf > 2/V) or (Tf <= V/A + 1/V):    # Tf <= V/A + 1/V  ʱ J<=0
                s = 0
                print("input error")
                return s
        else:
            if (Tf > 2*V / A + 1/V) or (Tf <= V/A + 1/V):
                s = 0
                print("input error")
                return s

    J = (A**2 * V) / (Tf*V*A - V**2 - A)
    T[0] = A / J
    T[1] = V / A - A / J  # T(2) = V / A - T(1)
    T[2] = T[0]
    T[3] = Tf - 2 * A / J - 2 * V / A    # T(4) = Tf - 4*T(1) - 2*T(2)
    T[4] = T[2]
    T[5] = T[1]
    T[6] = T[0]

    if t >= 0 and t <= T[0]:
        s = 1/6 * J * t**3
    elif t > T[0] and t <= T[0]+T[1]:
        dt = t - T[0]
        s = 1/2 * A * dt**2 + A**2/(2*J) * dt+ A**3/(6*J**2)
    elif t > T[0]+T[1] and t <= T[0]+T[1]+T[2]:
        dt = t - T[0] - T[1]
        s = -1/6*J*dt**3 + 1/2*A*dt**2 + (A*T[1] + A**2/(2*J))*dt + 1/2*A*T[1]**2 + A**2/(2*J)*T[1] + A**3/(6*J**2)
    elif t > T[0]+T[1]+T[2] and t <= T[0]+T[1]+T[2]+T[3]:
        dt = t - T[0] - T[1] - T[2]
        s = V*dt + (-1/6*J*T[2]**3) + 1/2*A*T[2]**2 + (A*T[1] + A**2/(2*J))*T[2] + 1/2*A*T[1]**2 + A**2/(2*J)*T[1] + A**3/(6*J**2)
    elif t > T[0]+T[1]+T[2]+T[3] and t <= T[0]+T[1]+T[2]+T[3]+T[4]:
        t_temp = Tf - t 
        dt = t_temp - T[0] - T[1]
        s = -1/6*J*dt**3 + 1/2*A*dt**2 + (A*T[1] + A**2/(2*J))*dt + 1/2*A*T[1]**2 + A**2/(2*J)*T[1] + A**3/(6*J**2)
        s = 1 - s
    elif t > T[0]+T[1]+T[2]+T[3]+T[4] and t <= T[0]+T[1]+T[2]+T[3]+T[4]+T[5]:
        t_temp = Tf - t 
        dt = t_temp - T[0]
        s = 1/2 * A * dt**2 + A**2/(2*J) * dt + A**3/(6*J**2)
        s = 1 - s  
    elif t > T[0]+T[1]+T[2]+T[3]+T[4]+T[5] and t <= T[0]+T[1]+T[2]+T[3]+T[4]+T[5]+T[6] + 1e5:
        t_temp = Tf - t 
        s = 1/6 * J * t_temp**3
        s = 1 - s     
    return s
    
def computeTrj(Start, End, Vel, Acc, Time):

    if(Vel == 0):
       Vel = 2 * (Start - End) /Time

    if(Acc == 0):
        Acc = 2 * (Start - End) / (Time ** 2)
        
    N = 2000
    ThetaStart = Start
    ThetaEnd = End
    VTheta = Vel
    ATheta = Acc
    Tf = Time

    v = VTheta/(ThetaEnd - ThetaStart)
    a = ATheta/(ThetaEnd - ThetaStart)
    v = abs(v)
    a = abs(a)

    Theta = np.zeros(N)
    s = np.zeros(N)
    sd = np.zeros(N)
    sdd = np.zeros(N)

    TF,V,A,J,T = SCurvePara(Tf, v, a)

    dTf = TF - Tf
    dV = V - v
    dA = A - a

    t = np.linspace(0, TF, N)
    dt = t[1] - t[0]

    for i in range(0, N):
        if 1 == N:
            a = a
            break
        s[i] = SCurveScaling(t[i],V,A,J,T,TF)
        Theta[i] = ThetaStart + s[i] * (ThetaEnd - ThetaStart)
        if i>0:
            sd[i-1] = (s[i] - s[i-1]) / dt
        if i>1:
            sdd[i-2] = (sd[i-1] - sd[i-2]) / dt

    # plt.subplot(2,2,1)
    # plt.plot(t, Theta, linewidth=1)

    # plt.subplot(2,2,2)
    # plt.plot(t, s, linewidth=1)

    # plt.subplot(2,2,3)
    # plt.plot(t, sd, linewidth=1)

    # plt.subplot(2,2,4)
    # plt.plot(t, sdd, linewidth=1)

    # plt.show()


    return Theta, s, sd, sdd

# if __name__ == "__main__":
#     computeTrj( 0,0.14 ,1, 1, 1)
    # computeTrj(10,90,2,4,10)

# -*- coding: utf-8 -*-
"""
Created on Mon Jan 31 16:49:29 2022

@author: Marshall
"""

from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt

'''
dL/dt = -kon1*L*P1+koff1*LP1-kon2*L*P2+koff*LP2
dLP1/dt = kon1*L*P1-koff1*LP1
dLP2/dt = kon2*L*P2-koff1*LP2
dP1/dt = -kon1*L*P1+koff1*LP1
dP2/dt = -kon2*L*P2+koff2*LP2
'''

#rate constants
kon1 = 49e3 #1/M/s
koff1 = 0.015e-3 #1/s
kon2 = 150e3 #1/M/s
koff2 = 5e-3 #1/s

#initial values
L_0 = 1e10/(6.022e23*25e-6) #M
P1_0 = 100e-9 #M
P2_0 = 1000e-9 #M
LP1_0 = 0 #M
LP2_0 = 0 #M

#start time
t_0 = 0 #s
t_end = 3600*5 #s

#number of time steps
n=1000

def rhs(t,c):
    return [-kon1*c[0]*c[1]+koff1*c[3]-kon2*c[0]*c[2]+koff2*c[4],
            -kon1*c[0]*c[1]+koff1*c[3],
            -kon2*c[0]*c[2]+koff2*c[4],
            kon1*c[0]*c[1]-koff1*c[3],
            kon2*c[0]*c[2]-koff2*c[4]]

res = solve_ivp(rhs,[t_0,t_end],[L_0,P1_0,P2_0,LP1_0,LP2_0],t_eval=np.linspace(t_0,t_end,n),method='BDF')

plt.plot(res.t,res.y.T[:,3],label='LP1')
plt.plot(res.t,res.y.T[:,4],label='LP2')
plt.plot(res.t,res.y.T[:,0],label='L')
plt.legend()
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 10 15:56:32 2022

@author: marsh
"""
from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from library_ode import * 
# =============================================================================
# plot sigmoid with fit
# ===============================================5==============================
#rate constants

# #solution phase:
# kon1 = 5e5 #1/M/s
# koff1 = 1.5e-3 #1/s

# #BSD phase:
# kon2 = 5e5 #1/M/s
# koff2 = 1.5e-5 #1/s

# #concentrations of solution peptide
# lower = -10
# upper = -6
# n_points = 8
# peptide_concs = np.logspace(lower,upper,n_points)

# #concentration of mcl-1
# mcl1_conc = 1e-8

# BSD_conc = 1e10/(6.022e23*100e-6)


#
# =============================================================================
# analyze data
# =============================================================================
#concentration of mcl-1
mcl1_conc = 1e-8
BSD_conc = 1e10/(6.022e23*100e-6)
lower = -10
upper = -6
n_points = 8
peptide_concs = np.logspace(lower,upper,n_points)

##MCL1
#solution phase:
kon1 = 5e5 #1/M/s
koff1 = 1.5e-2 #1/s

#BSD phase:
kon2 = 5e5 #1/M/s
koff2 = 1.5e-5 #1/s

fig,axes = plt.subplots(nrows=2,ncols=3,figsize=(20,10))
axes = np.ravel(axes)

L_1,P1_1,_,LP1_1,_,_=plot(kon1,koff1,kon2,koff2,mcl1_conc,peptide_concs[-1],0,0,0,plots=['L','P1','LP1'],labels=['L','P1','LP1'],figure=True,ax=axes[0])
plot(kon1,koff1,kon2,koff2,L_1,P1_1,BSD_conc,LP1_1,0,plots=['L','P1','LP1','P2','LP2'],labels=['L','P1','LP1','P2','LP2'],figure=True,ax=axes[1])
plt.plot([900,900],[1e-13,1e-6],linestyle='--')
plotKd(kon1,koff1,kon2,koff2,peptide_concs,BSD_conc,mcl1_conc,ax=axes[2])


##BFL1
#solution phase:
kon1 = 5e5 #1/M/s
koff1 = 0.75e-1 #1/s

#BSD phase:
kon2 = 5e5 #1/M/s
koff2 = 0.75e-2 #1/s

L_1,P1_1,_,LP1_1,_,_=plot(kon1,koff1,kon2,koff2,mcl1_conc,peptide_concs[-1],0,0,0,plots=['L','P1','LP1'],labels=['L','P1','LP1'],figure=True,ax=axes[3])
plot(kon1,koff1,kon2,koff2,L_1,P1_1,BSD_conc,LP1_1,0,plots=['L','P1','LP1','P2','LP2'],labels=['L','P1','LP1','P2','LP2'],figure=True,ax=axes[4])
plt.plot([900,900],[1e-13,1e-6],linestyle='--')
plotKd(kon1,koff1,kon2,koff2,peptide_concs,BSD_conc,mcl1_conc,ax=axes[5])
plt.tight_layout()

# #display level sensitivity analysis
# plt.figure()
# plt.scatter(np.logspace(8,13,6),[plotKd(kon1,koff1,kon2,koff2,peptide_concs,d_l/6.022e23/25e-6,mcl1_conc,toPlot=False) for d_l in np.logspace(8,13,6)])
# plt.xscale('log')
# plt.yscale('log')
# plt.xlabel('# of peptides displayed on bacteria')
# plt.ylabel('fit Kd')

# #mcl1 concentration sensitivity analysis
# plt.figure()
# plt.scatter(np.logspace(-12,-6,10),[plotKd(kon1,koff1,kon2,koff2,peptide_concs,BSD_conc,mc,toPlot=False) for mc in np.logspace(-12,-6,10)])
# plt.xscale('log')
# plt.yscale('log')
# plt.xlabel('Mcl-1 conc [M]')
# plt.ylabel('fit Kd')

# #bsd affinity sensitivity analysis
# plt.figure()
# plt.scatter(np.logspace(-6,-1,12),[plotKd(kon1,koff1,kon2,k2,peptide_concs,BSD_conc,mcl1_conc,toPlot=False) for k2 in np.logspace(-6,-1,12)])
# plt.xscale('log')
# plt.yscale('log')
# plt.xlabel('BSD-peptide k_off [1/s]')
# plt.ylabel('fit Kd')

# =============================================================================
# plot raw concentrations
# =============================================================================
# for peptide_conc in peptide_concs:
    
#     L_1,P1_1,P2_1,LP1_1,LP2_1 = plot(kon1,koff1,kon2,koff2,mcl1_conc,peptide_conc,0,0,0,toPlot=False)
#     print('P1_0: ' + str(peptide_conc) + ' P1_1: ' + str(P1_1))
    
#     plot(kon1,koff1,kon2,koff2,L_1,P1_1,1e10/(6.022e23*25e-6),LP1_1,LP2_1,label=peptide_conc)

# plt.legend()

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
import pandas as pd


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
# reproduce experimental observation where mcl-1 doesn't block due to very high BSD affinity
# =============================================================================
# #concentration of mcl-1
# mcl1_conc = 1e-9
# BSD_conc = 1e10/(6.022e23*100e-6)
# lower = -10
# upper = -6
# n_points = 8
# peptide_concs = np.logspace(lower,upper,n_points)

# ##MCL1
# #solution phase:
# kon1 = 5e5 #1/M/s
# koff1 = 1.5e-2 #1/s

# #BSD phase:
# kon2 = 5e5 #1/M/s
# koff2 = 1.5e-5 #1/s

# fig,axes = plt.subplots(nrows=2,ncols=3,figsize=(20,10))
# axes = np.ravel(axes)

# L_1,P1_1,_,LP1_1,_,_=plot(kon1,koff1,kon2,koff2,mcl1_conc,peptide_concs[-1],0,0,0,plots=['L','P1','LP1'],labels=['L','P1','LP1'],figure=True,ax=axes[0])
# plot(kon1,koff1,kon2,koff2,L_1,P1_1,BSD_conc,LP1_1,0,plots=['L','P1','LP1','P2','LP2'],labels=['L','P1','LP1','P2','LP2'],figure=True,ax=axes[1])
# plt.plot([900,900],[1e-13,1e-6],linestyle='--')
# plotKd(kon1,koff1,kon2,koff2,peptide_concs,BSD_conc,mcl1_conc,ax=axes[2])


# ##BFL1
# #solution phase:
# kon1 = 5e5 #1/M/s
# koff1 = 0.75e-1 #1/s

# #BSD phase:
# kon2 = 5e5 #1/M/s
# koff2 = 0.75e-2 #1/s

# L_1,P1_1,_,LP1_1,_,_=plot(kon1,koff1,kon2,koff2,mcl1_conc,peptide_concs[-1],0,0,0,plots=['L','P1','LP1'],labels=['L','P1','LP1'],figure=True,ax=axes[3])
# plot(kon1,koff1,kon2,koff2,L_1,P1_1,BSD_conc,LP1_1,0,plots=['L','P1','LP1','P2','LP2'],labels=['L','P1','LP1','P2','LP2'],figure=True,ax=axes[4])
# plt.plot([900,900],[1e-13,1e-6],linestyle='--')
# plotKd(kon1,koff1,kon2,koff2,peptide_concs,BSD_conc,mcl1_conc,ax=axes[5])
# plt.tight_layout()


# =============================================================================
# fit experimental data for SORT-1-PE and SORT-5-PE
# =============================================================================
folder = 'C:\\Users\\marsh\\Dropbox (University of Michigan)\\Michigan\\assorted antigens\\bcl2\\flow experiments\\9aug22\\'
filename = 'data_10aug22.xlsx'
data = pd.read_excel(folder+filename)

# #concentration of mcl-1
mcl1_conc = 1e-8
BSD_conc = 1e10/(6.022e23*100e-6)
lower = -10
upper = -6
n_points = 8
peptide_concs = data['SORT-1-PE']*1e-9
peptide_concs = peptide_concs.values

##MCL1
#solution phase:
kon1 = 5e5 #1/M/s
koff1 = 1.5e-2 #1/s

#BSD phase:
kon2 = 5e5 #1/M/s
koff2 = 1.5e-5 #1/s

#BSD affinities
BSD_affinities = [30e-12,5e-9,5e-9,15e-9,15e-9]

fig,axes = plt.subplots(ncols=5,nrows=2,figsize=(20,10),sharex=True,sharey=True)
axes = np.ravel(axes)
for index,ax,col in zip(range(10),axes,np.concatenate((data.columns[1:6],data.columns[7:12]))):
    fitKd,ax = fitToExperiment(peptide_concs,data[col],ax=ax)
    ax.set_title(col + ', fit Ki: ' + str(np.round(getKiFromIC50(0.5,BSD_affinities[index % 5],BSD_conc,fitKd),9)) + 'nM')
plt.tight_layout()



# =============================================================================
# see how mcl-1 could be measured using a lower concentration of mcl-1
# =============================================================================
fig,axes = plt.subplots(nrows=2,ncols=3,figsize=(20,10))
axes = np.ravel(axes)

L_1,P1_1,_,LP1_1,_,_=plot(kon1,koff1,kon2,koff2,mcl1_conc,peptide_concs[0],0,0,0,plots=['L','P1','LP1'],labels=['L','P1','LP1'],figure=True,ax=axes[0])
plot(kon1,koff1,kon2,koff2,L_1,P1_1,BSD_conc,LP1_1,0,plots=['L','P1','LP1','P2','LP2'],labels=['L','P1','LP1','P2','LP2'],figure=True,ax=axes[1])
plotKd(kon1,koff1,kon2,koff2,peptide_concs,BSD_conc,mcl1_conc,ax=axes[2])

L_1,P1_1,_,LP1_1,_,_=plot(kon1,koff1,kon2,koff2,mcl1_conc/10,peptide_concs[0],0,0,0,plots=['L','P1','LP1'],labels=['L','P1','LP1'],figure=True,ax=axes[3])
plot(kon1,koff1,kon2,koff2,L_1,P1_1,BSD_conc,LP1_1,0,plots=['L','P1','LP1','P2','LP2'],labels=['L','P1','LP1','P2','LP2'],figure=True,ax=axes[4])
plotKd(kon1,koff1,kon2,koff2,peptide_concs,BSD_conc,mcl1_conc/10,ax=axes[5])

# =============================================================================
# produce plot to observe effects of lower protein concentration
# =============================================================================
fig,axes = plt.subplots(nrows=3,ncols=3,figsize=(25,15))
axes = np.ravel(axes)

mcl1_concs = np.array([10,1,0.1])*1e-9

for index,mcl1_conc in enumerate(mcl1_concs):
    L_1,P1_1,_,LP1_1,_,_=plot(kon1,koff1,kon2,koff2,mcl1_conc,peptide_concs[0],0,0,0,plots=['L','P1','LP1'],labels=['L','P1','LP1'],figure=True,ax=axes[(index)*3])
    plot(kon1,koff1,kon2,koff2,L_1,P1_1,BSD_conc,LP1_1,0,plots=['L','P1','LP1','P2','LP2'],labels=['L','P1','LP1','P2','LP2'],figure=True,ax=axes[(index)*3+1])
    plotKd(kon1,koff1,kon2,koff2,peptide_concs,BSD_conc,mcl1_conc,ax=axes[index*3+2])

plt.tight_layout()

# =============================================================================
# sensitivity analysis
# =============================================================================
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

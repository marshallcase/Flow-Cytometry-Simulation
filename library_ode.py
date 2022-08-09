# -*- coding: utf-8 -*-
"""
Created on Mon Jan 31 16:49:29 2022

@author: Marshall
"""

from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def rhs(t,c):
    '''
    dL/dt = -kon1*L*P1+koff1*LP1-kon2*L*P2+koff*LP2
    dLP1/dt = kon1*L*P1-koff1*LP1
    dLP2/dt = kon2*L*P2-koff1*LP2
    dP1/dt = -kon1*L*P1+koff1*LP1
    dP2/dt = -kon2*L*P2+koff2*LP2
    '''
    return [-kon1*c[0]*c[1]+koff1*c[3]-kon2*c[0]*c[2]+koff2*c[4], #free ligand
            -kon1*c[0]*c[1]+koff1*c[3], #protein 1
            -kon2*c[0]*c[2]+koff2*c[4], #protein 2
            kon1*c[0]*c[1]-koff1*c[3], #ligand-protein1
            kon2*c[0]*c[2]-koff2*c[4]] #ligand-protein2

def plot(kon1,koff1,kon2,koff2,L_0,P1_0,P2_0,LP1_0,LP2_0,label=None,t_0=0,t_end=3600*5,n=1000,toPlot=True):
    '''
    Inputs:
        kon1 1/M/s
        koff1 1/s
        kon2 1/M/s
        koff 1/s
        L_0: L | t=0
        P1_0: P1 | t=0
        P2_0,LP1_0, LP2_0 etc.
        label: label for plot (in the loop, use concentration)
        t_0: default 0 (seconds)
        t_end: t | t=end, default 3600*5
        n: number of points, default 1000
        toPlot: plots if True
    Outputs:
        L
        P1
        P2
        LP1
        LP2
    '''
    res = solve_ivp(rhs,[t_0,t_end],[L_0,P1_0,P2_0,LP1_0,LP2_0],t_eval=np.linspace(t_0,t_end,n),method='BDF')
    if toPlot:
        # plt.figure()
        # plt.plot(res.t,res.y.T[:,1],label='P1')
        #plt.plot(res.t,res.y.T[:,2],label=b2)
        # plt.plot(res.t,res.y.T[:,3],label='LP1')
        plt.plot(res.t,res.y.T[:,4],label=label)
        # plt.plot(res.t,res.y.T[:,0],label='L')
        #plt.yscale('log')
        # plt.legend()

    return res.y.T[-1,0],res.y.T[-1,1],res.y.T[-1,2],res.y.T[-1,3],res.y.T[-1,4]

def sigmoid(t,a,b,c):
    return  a/(b+np.exp(-c*t))

def sigmoid2(x,top,bottom,kd):
    return bottom + x*(top-bottom)/(kd+x)

def fitToExperiment(exp_x,exp_y):
    '''
    Inputs:
        exp_x: concentrations of experimental data
        exp_y: % block of experimental data
    Outputs:
        fit_kd: fit Kd
    '''
    plt.figure()
    top = 0
    bottom = 1
    popt, pcov = curve_fit(lambda x,kd: sigmoid2(x,top,bottom,kd), exp_x, exp_y,method='trf')
    plt.scatter(exp_x,exp_y,label='exp')
    plt.xscale('log')
    x = np.logspace(np.log10(min(exp_x)),np.log10(max(exp_x)*10),n_points*100)
    y=sigmoid2(x,top,bottom,popt[0])
    plt.plot(x,y,label='exp fit')
    plt.legend()
    fit_kd = x[np.argmin(np.abs(y-0.5))]
    plt.plot([fit_kd,fit_kd],[0,0.5],linestyle='--',color='tab:blue')
    plt.plot([min(x),x[np.argmin(np.abs(y-0.5))]],[0.5,0.5],color='tab:blue',linestyle='--')
    plt.text(x[np.argmin(np.abs(y-0.96))],y[np.argmin(np.abs(y-0.4))],s='fit IC50: ' +str(np.around(fit_kd,9)))
    plt.xlim([x[0]*0.5,x[-1]*1.5])
    plt.ylim([-0.1,1.1])
    return fit_kd

def plotKd(kon1,koff1,kon2,koff2,peptide_concs,BSD_conc,mcl1_conc,toPlot=True):
    '''
    Inputs
        kon1
        koff1
        kon2
        koff2
        peptide_concs
        BSD_conc
        mcl1_conc
        toPlot: plot if True
    
    Outputs:
        fit_kd: fit Kd
    '''
    #get phase 1 concentrations (blocking phase) based on peptide concentrations
    y1 = np.hstack([plot(kon1,koff1,kon2,koff2,mcl1_conc,peptide_conc,0,0,0,toPlot=False) for peptide_conc in peptide_concs])
    y1 = np.reshape(y1,(8,5))
    L_1 = y1[:,0]
    P1_1 = y1[:,1]
    P2_1 = y1[:,2]
    LP1_1 = y1[:,3]
    LP2_1 = y1[:,4]
    
    #get BSD-peptide-mcl1 concentrations (phase 2)
    PL2 = [plot(kon1,koff1,kon2,koff2,L,P1,BSD_conc,LP1,LP2,toPlot=False) for L,P1,P2,LP1,LP2 in zip(L_1,P1_1,P2_1,LP1_1,LP2_1)]
    PL2 = np.hstack(PL2)
    PL2 = np.reshape(PL2,(8,5))
    PL2 = PL2[:,4]
    PL2 = PL2/max(PL2)
    
    #calculate and plot Kd's of BSD and solution peptide
    Kd_BSD = koff2/kon2
    Kd_solution = koff1/kon1
    
    #curve_fit and  assume sigmoid goes from y=1 as x-> 0 and y=0 as x-> infinity
    top = 0
    bottom = 1
    lower = np.log10(min(peptide_concs))
    upper = np.log10(max(peptide_concs))
    curve_fit_x = np.logspace(lower,upper,n_points*100)
    popt, pcov = curve_fit(lambda x,kd: sigmoid2(x,top,bottom,kd), peptide_concs, PL2,method='trf')
    y = sigmoid2(curve_fit_x,top,bottom,popt[0])

    fit_kd = curve_fit_x[np.argmin(np.abs(y-0.5))]
    #Kd_obs = Kd_solution*(1+BSD_conc/Kd_BSD)
    
    if toPlot:
        #set up scatter plot with simulated experimental data points
        plt.figure()
        plt.scatter(peptide_concs,PL2)
        plt.xscale('log')
        plt.xlabel('conc of solution-peptide [M]')
        plt.ylabel('fraction BSD-peptide bound')
        #get Kd's of BSD and solution peptide
        plt.title('Kd,BSD: ' + str(Kd_BSD) + ', Kd,solution: ' + str(Kd_solution))
        #plot the curve_fit
        plt.plot(curve_fit_x,y)
        #set up axes and label Kd on graph
        plt.plot([fit_kd,fit_kd],[0,0.5],linestyle='--',color='tab:blue')
        plt.plot([min(curve_fit_x),curve_fit_x[np.argmin(np.abs(y-0.5))]],[0.5,0.5],color='tab:blue',linestyle='--')
        plt.text(curve_fit_x[np.argmin(np.abs(y-0.96))],y[np.argmin(np.abs(y-0.4))],s='fit IC50: ' +str(np.around(fit_kd,9)))
        plt.xlim([curve_fit_x[0]*0.5,curve_fit_x[-1]*1.5])
        plt.ylim([-0.1,1.1])
        
    return fit_kd

def getKiFromIC50(f0,Kd_BSD,BSD_conc,IC50):
    '''

    Parameters
    ----------
    f0 : float
        fraction bound. At IC50, f0=0.5.
    Kd_BSD : float
        Kd of the BSD peptide-protein interaction. Ususally around 1nM.
    BSD_conc : float
        Concentration of BSD peptide. Usually 1e10/(6.023e23*25e-6)
    IC50 : float
        fit IC50

    Returns
    -------
    Ki: float
        concentration of inhibtor for 50% inhibition 

    '''
    return (IC50/((f0*Kd_BSD)/(1-f0)/(2-f0)+f0*BSD_conc/(2-f0))-1)*Kd_BSD*f0/(2-f0)
#rate constants



# #initial values
# L_0 =  3.52e-10#M
# # L_0 = 1e-6 #M
# P1_0 = 9.003e-7 #M
# P2_0 = 1e10/(6.022e23*25e-6) #M
# LP1_0 = 9.96e-8 #M
# LP2_0 = 0 #M

# #start time
# t_0 = 0 #s
# t_end = 3600*5 #s

# #number of time steps
# n=1000





# =============================================================================
# plot raw concentrations
# =============================================================================
# for b2 in b2_conc:
    
#     L_1,P1_1,P2_1,LP1_1,LP2_1 = plot(kon1,koff1,kon2,koff2,1e-8,b2,0,0,0,toPlot=False)
#     print('P1_0: ' + str(b2) + ' P1_1: ' + str(P1_1))
    
#     plot(kon1,koff1,kon2,koff2,L_1,P1_1,1e10/(6.022e23*25e-6),LP1_1,LP2_1,label=b2)

# plt.legend()


# =============================================================================
# plot sigmoid with fit
# ===============================================5==============================
#rate constants

#solution phase:
kon1 = 5e5 #1/M/s
koff1 = 1.5e-4 #1/s

#BSD phase:
kon2 = 5e5 #1/M/s
koff2 = 1.5e-4 #1/s

#concentrations of solution peptide
lower = -10
upper = -5
n_points = 8
peptide_concs = np.logspace(lower,upper,n_points)

#concentration of mcl-1
mcl1_conc = 1e-8

BSD_conc = 1e10/(6.022e23*25e-6)

#display level sensitivity analysis
plt.figure()
plt.scatter(np.logspace(8,13,6),[plotKd(kon1,koff1,kon2,koff2,peptide_concs,d_l/6.022e23/25e-6,mcl1_conc,toPlot=False) for d_l in np.logspace(8,13,6)])
plt.xscale('log')
plt.yscale('log')
plt.xlabel('# of peptides displayed on bacteria')
plt.ylabel('fit Kd')

#mcl1 concentration sensitivity analysis
plt.figure()
plt.scatter(np.logspace(-12,-6,10),[plotKd(kon1,koff1,kon2,koff2,peptide_concs,BSD_conc,mc,toPlot=False) for mc in np.logspace(-12,-6,10)])
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Mcl-1 conc [M]')
plt.ylabel('fit Kd')

#bsd affinity sensitivity analysis
plt.figure()
plt.scatter(np.logspace(-6,-1,12),[plotKd(kon1,koff1,kon2,k2,peptide_concs,BSD_conc,mcl1_conc,toPlot=False) for k2 in np.logspace(-6,-1,12)])
plt.xscale('log')
plt.yscale('log')
plt.xlabel('BSD-peptide k_off [1/s]')
plt.ylabel('fit Kd')

# -*- coding: utf-8 -*-
"""
Created on Mon Jan 31 16:49:29 2022

@author: Marshall
"""

from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def rhs(t,c,kon1,koff1,kon2,koff2):
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

def plot(kon1,koff1,kon2,koff2,L_0,P1_0,P2_0,LP1_0,LP2_0,plots=[],labels=[],figure=False,t_0=0,t_end=3600*5,n=1000,ax=None):
    '''
    Inputs:
        kon1 1/M/s
        koff1 1/s
        kon2 1/M/s
        koff 1/s
        L_0: L | t=0
        P1_0: P1 | t=0
        P2_0,LP1_0, LP2_0 etc.
        plots: which plots to plot (array subset of 'L','P1','P2','LP1','LP2')
        labels: labels for corresponding plot
        figure: make new figure (True)
        t_0: default 0 (seconds)
        t_end: t | t=end, default 3600*5
        n: number of points, default 1000
        ax: substitute ax object for plotting in subplots, defaults to None (create new plot)
    Outputs:
        L at t=end
        P1 at t=end
        P2 at t=end
        LP1 at t=end
        LP2 at t=end
    '''
    res = solve_ivp(lambda t,c: rhs(t,c,kon1,koff1,kon2,koff2),[t_0,t_end],[L_0,P1_0,P2_0,LP1_0,LP2_0],t_eval=np.linspace(t_0,t_end,n),method='BDF')
    if ax is None:
        if len(plots) > 0:
            if figure:
                plt.figure()
                plt.yscale('log')
    
    
            for plot,label in zip(plots,labels):
                if plot == 'L':
                    plt.plot(res.t,res.y.T[:,0],label=label)
                elif plot =='P1':
                    plt.plot(res.t,res.y.T[:,1],label=label)
                elif plot =='P2':
                    plt.plot(res.t,res.y.T[:,2],label=label)
                elif plot == 'LP1':
                    plt.plot(res.t,res.y.T[:,3],label=label)
                elif plot == 'LP2':
                    plt.plot(res.t,res.y.T[:,4],label=label)
                
            if figure:
                plt.legend(loc='lower right')
                #calculate and plot Kd's of BSD and solution peptide
                Kd_BSD = koff2/kon2
                Kd_solution = koff1/kon1
                plt.title('Kd,BSD: ' + str(Kd_BSD) + ', Kd,solution: ' + str(Kd_solution))
                plt.ylim([1e-13,5e-6])
                plt.ylabel('Concentration [M]')
                plt.xlabel('Time (sec)')
            
        return res.y.T[-1,0],res.y.T[-1,1],res.y.T[-1,2],res.y.T[-1,3],res.y.T[-1,4]
    else:
        if len(plots) > 0:
            if figure:
                ax.set_yscale('log')
    
    
            for plot,label in zip(plots,labels):
                if plot == 'L':
                    ax.plot(res.t,res.y.T[:,0],label=label)
                elif plot =='P1':
                    ax.plot(res.t,res.y.T[:,1],label=label)
                elif plot =='P2':
                    ax.plot(res.t,res.y.T[:,2],label=label)
                elif plot == 'LP1':
                    ax.plot(res.t,res.y.T[:,3],label=label)
                elif plot == 'LP2':
                    ax.plot(res.t,res.y.T[:,4],label=label)
                
            if figure:
                ax.legend(loc='lower right')
                #calculate and plot Kd's of BSD and solution peptide
                Kd_BSD = koff2/kon2
                Kd_solution = koff1/kon1
                ax.set_title('Kd,BSD: ' + str(Kd_BSD) + ', Kd,solution: ' + str(Kd_solution))
                ax.set_ylim([1e-13,5e-6])
                ax.set_ylabel('Concentration [M]')
                ax.set_xlabel('Time (sec)')
        return res.y.T[-1,0],res.y.T[-1,1],res.y.T[-1,2],res.y.T[-1,3],res.y.T[-1,4],ax

    

def sigmoid(t,a,b,c):
    return  a/(b+np.exp(-c*t))

def sigmoid2(x,top,bottom,kd):
    return bottom + x*(top-bottom)/(kd+x)

def fitToExperiment(exp_x,exp_y,ax=None):
    '''
    Inputs:
        exp_x: concentrations of experimental data
        exp_y: % block of experimental data
        ax: substitute ax object for plotting in subplots, defaults to None (create new plot)
    Outputs:
        fit_kd: fit Kd
    '''
    top = 0
    bottom = 1
    popt, pcov = curve_fit(lambda x,kd: sigmoid2(x,top,bottom,kd), exp_x, exp_y,method='trf')
    n_points = 8
    x = np.logspace(np.log10(min(exp_x)),np.log10(max(exp_x)*10),n_points*100)
    y=sigmoid2(x,top,bottom,popt[0])
    fit_kd = x[np.argmin(np.abs(y-0.5))]

    if ax is None:
        plt.figure()
        plt.scatter(exp_x,exp_y,label='exp')
        plt.xscale('log')
        plt.plot(x,y,label='exp fit')
        plt.legend()
        plt.plot([fit_kd,fit_kd],[0,0.5],linestyle='--',color='tab:blue')
        plt.plot([min(x),x[np.argmin(np.abs(y-0.5))]],[0.5,0.5],color='tab:blue',linestyle='--')
        plt.text(x[np.argmin(np.abs(y-0.96))],y[np.argmin(np.abs(y-0.4))],s='fit IC50: ' +str(np.around(fit_kd,9)))
        plt.xlim([x[0]*0.5,x[-1]*1.5])
        plt.ylim([-0.2,1.2])
        plt.xlabel('Fraction Blocked')
        plt.ylabel('Concentration [M]')
        return fit_kd
        
    else:
        ax.scatter(exp_x,exp_y,label='exp')
        ax.set_xscale('log')
        ax.plot(x,y,label='exp fit')
        ax.legend()
        ax.plot([fit_kd,fit_kd],[0,0.5],linestyle='--',color='tab:blue')
        ax.plot([min(x),x[np.argmin(np.abs(y-0.5))]],[0.5,0.5],color='tab:blue',linestyle='--')
        ax.text(x[np.argmin(np.abs(y-0.96))],y[np.argmin(np.abs(y-0.4))],s='fit IC50: ' +str(np.around(fit_kd,9)))
        ax.set_xlim([x[0]*0.5,x[-1]*1.5])
        ax.set_ylim([-0.2,1.2])
        ax.set_xlabel('Fraction Blocked')
        ax.set_ylabel('Concentration [M]')
        return fit_kd, ax

def plotKd(kon1,koff1,kon2,koff2,peptide_concs,BSD_conc,mcl1_conc,ax=None):
    '''
    Inputs
        kon1
        koff1
        kon2
        koff2
        peptide_concs
        BSD_conc
        mcl1_conc
        ax: substitute ax object for plotting in subplots, defaults to None (create new plot)
    Outputs:
        fit_kd: fit Kd
    '''
    #get phase 1 concentrations (blocking phase) based on peptide concentrations
    y1 = np.hstack([plot(kon1,koff1,kon2,koff2,mcl1_conc,peptide_conc,0,0,0) for peptide_conc in peptide_concs])
    y1 = np.reshape(y1,(len(peptide_concs),5))
    L_1 = y1[:,0]
    P1_1 = y1[:,1]
    P2_1 = y1[:,2]
    LP1_1 = y1[:,3]
    LP2_1 = y1[:,4]
    
    #get BSD-peptide-mcl1 concentrations (phase 2)
    PL2 = [plot(kon1,koff1,kon2,koff2,L,P1,BSD_conc,LP1,LP2) for L,P1,P2,LP1,LP2 in zip(L_1,P1_1,P2_1,LP1_1,LP2_1)]
    PL2 = np.hstack(PL2)
    PL2 = np.reshape(PL2,(len(peptide_concs),5))
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
    n_points = 8
    curve_fit_x = np.logspace(lower,upper,n_points*100)
    popt, pcov = curve_fit(lambda x,kd: sigmoid2(x,top,bottom,kd), peptide_concs, PL2,method='trf')
    y = sigmoid2(curve_fit_x,top,bottom,popt[0])

    fit_kd = curve_fit_x[np.argmin(np.abs(y-0.5))]
    #Kd_obs = Kd_solution*(1+BSD_conc/Kd_BSD)
    

    #set up scatter plot with simulated experimental data points
    if ax is None:
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
    else:
        ax.scatter(peptide_concs,PL2)
        ax.set_xscale('log')
        ax.set_xlabel('conc of solution-peptide [M]')
        ax.set_ylabel('fraction BSD-peptide bound')
        #get Kd's of BSD and solution peptide
        ax.set_title('Kd,BSD: ' + str(Kd_BSD) + ', Kd,solution: ' + str(Kd_solution))
        #plot the curve_fit
        ax.plot(curve_fit_x,y)
        #set up axes and label Kd on graph
        ax.plot([fit_kd,fit_kd],[0,0.5],linestyle='--',color='tab:blue')
        ax.plot([min(curve_fit_x),curve_fit_x[np.argmin(np.abs(y-0.5))]],[0.5,0.5],color='tab:blue',linestyle='--')
        ax.text(curve_fit_x[np.argmin(np.abs(y-0.96))],y[np.argmin(np.abs(y-0.4))],s='fit IC50: ' +str(np.around(fit_kd,9)))
        ax.set_xlim([curve_fit_x[0]*0.5,curve_fit_x[-1]*1.5])
        ax.set_ylim([-0.1,1.1])

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
    return (IC50/((f0*Kd_BSD)/(1-f0)/(2-f0)+f0*BSD_conc/2)-1)*Kd_BSD*f0/(2-f0)



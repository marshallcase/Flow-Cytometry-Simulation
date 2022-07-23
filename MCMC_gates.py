# -*- coding: utf-8 -*-
"""
Created on Tue Oct 19 11:38:15 2021

@author: marsh
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.path as mplPath
import math
from scipy.ndimage import gaussian_filter1d


def sigmoid(x):
    return 1 / (1 + math.exp(-x))


def perturbGates(yl,yr,xs,intensity=0.05,fraction=0.5):
    '''
    yl,yr,xs : coordinates of gates
    intensity: amount of perturbation (as a percentage)
    fraction: fraction of gates to perturb
    '''
    #concatenated coordinate vector
    coords = np.concatenate((yl,yr,xs))
    #binary perturbation vector
    overlap=True
    while overlap == True:
        bin_per = np.random.choice([0,1],size=(len(coords),),p=[fraction,1-fraction])
        #amount of perturbation vector
        am_per = np.random.uniform(1-intensity,1+intensity,len(coords))
        perturb = bin_per*am_per
        perturb[perturb==0]=1
        test = coords*perturb
        if (test[0] < test [1]) and (test[1] < test [2]) and (test[2] < test[3])\
            and (test[3] < test[4]):
            if (test[5] < test [6]) and (test[6] < test [7]) and (test[7] < test[8])\
                and (test[8] < test[9]):
                if test[10] < test[11]:
                    overlap=False
                    coords=test
    yl = coords[:5]
    yr = coords[5:10]
    xs = coords[10:12]
    return yl,yr,xs

def convertCoordsToGates(yl,yr,xs):
    '''
    convert yl, yr, and xs to gates array
    '''
    gates = [[[xs[0],yl[i]],[xs[1],yr[i]],[xs[1],yr[i+1]],[xs[0],yl[i+1]],[xs[0],yl[i]]] for i in range(len(yl)-1)]
    return gates


def getBinScore(gates,exp_signal,bind_signal):
    '''
    this function calculates the average bin number of a flow population
    gates = shape of the current gates
    exp_signal = vector of expression signals
    bind_signal = vector of binding signals
    '''
    scores = np.array([])
    interpolation_num = 9
    #caclulate the number of points that fall inside all the gates for the
    #denominator
    all_gate = [gates[0][0],gates[0][1],gates[-1][2],gates[-1][3],gates[0][0]]
    all_gate = getInterpolatedGates(all_gate,interpolation_num)
    bbPath_all = mplPath.Path(np.array(all_gate))
    num_points = bbPath_all.contains_points(np.array([exp_signal,bind_signal]).T).sum()
    
    #get the negative population, use as many positive, points
    exp_signal_neg = np.random.normal(100,100,num_points)
    bind_signal_neg = np.random.normal(100,100,num_points)
    num_points_neg = bbPath_all.contains_points(np.array([exp_signal_neg,bind_signal_neg]).T).sum()
    neg_score = num_points_neg / num_points
    for index,gate in enumerate(gates):
        bbPath = mplPath.Path(np.array(getInterpolatedGates(gate[:4],interpolation_num)))
        score = bbPath.contains_points(np.array([exp_signal,bind_signal]).T).sum()/num_points
        scores = np.append(scores,score)
    
    bin_score = sum([scores[i]*(i+1) for i in range(len(gates))])
    return bin_score,neg_score

def getEquidistantPoints(p1, p2, parts):
    '''
    linspace between (x1,y1) and (x2,y2)
    '''
    return np.linspace(p1[0], p2[0], parts+1),np.linspace(p1[1], p2[1], parts+1)

def getLogEquidistantPoints(p1, p2, parts):
    '''
    logspace between (x1,y1) and (x2,y2)
    '''
    return np.logspace(np.log10(p1[0]), np.log10(p2[0]), parts+1),np.logspace(np.log10(p1[1]), np.log10(p2[1]), parts+1)

def plotGates(gates):
    '''
    plot gates and label axes
    '''
    for gate in gates:
        x, y = zip(*gate)
        plt.plot(x,y,color='black')
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel('expression')
        plt.ylabel('binding')
        plt.xlim([10,1e7])
        plt.ylim([10,1e7])
                
def plotFlowData(ligand_concs,gates,Kd,num_points=100):
    '''
    plot a series of populations corresponding to the Kd and ligand concentrations
    vector
    '''
    labels = ['c='+str(ligand_conc) for ligand_conc in ligand_concs]
    bin_scores_df = pd.DataFrame(index=ligand_concs,columns=['score'])
    for ligand_conc,label in zip(ligand_concs,labels):
        frac_bound = (ligand_conc)/(ligand_conc+Kd)
        exp_level = np.random.lognormal(np.log(7000),np.log(5),num_points)
        
        exp_signal = exp_level + np.random.normal(100,100,num_points)
        bind_signal = np.random.lognormal(0,0.5,num_points)*exp_level*frac_bound+np.random.normal(100,100,num_points)
        bin_score = np.around(getBinScore(gates,exp_signal,bind_signal),2)
        bin_scores_df.loc[ligand_conc,'score']=bin_score
        plt.scatter(exp_signal,bind_signal,label=label+' score='+str(bin_score),alpha=0.4)
        plt.legend()
            
        
def getFlowData(ligand_conc,Kd,num_points):
    '''
    given a vector of ligand concentrations and the Kd, generate vectors of 
    expression and binding signals that are realistic distributions of flow data
    '''
    frac_bound = (ligand_conc)/(ligand_conc+Kd)
    #7000 is the scaling factor for bacterial surface display
    #5 is the variance in the display level
    exp_level = np.random.lognormal(np.log(7000),np.log(5),num_points)
    #100 is flow cytometer noise
    exp_signal = exp_level + np.random.normal(100,100,num_points)
    #lognormal adds noise to the dependence between binding signal and expression signal
    bind_signal = np.random.lognormal(0,0.5,num_points)*exp_level*frac_bound+np.random.normal(100,100,num_points)
    return exp_signal,bind_signal

def plotNegative(num_points=100):
    '''
    plot a negative population with realistic noise for a flow cytometer
    '''
    plt.scatter(np.random.normal(100,100,num_points),np.random.normal(100,100,num_points),label='negative',alpha=0.2)

def getMCMCScore(bin_scores_df,Kd):
    '''
    this function calculates the error between the average bin scores of the ligand
    concentration series and the "idealized" distribution, which is where the bin
    number directly correlates to the fraction bound. This is mathematically,
    average bin number = (max bin - min bin)(sigmoid(conc)) + (min bin)
    '''
    #TODO: penalize how many negatives end up in the gates
    #TODO: penalize how many positives end up not in the gates
    bin_scores_df['ideal'] = 3*(bin_scores_df.index/ (bin_scores_df.index+Kd))+1
    error = ((bin_scores_df['ideal']-bin_scores_df['score'])**2).sum()
    neg_error = ((bin_scores_df['neg_score']-0).sum())**2
    total_error = error + neg_error
    return total_error, neg_error

def getMinCoordinates(metrics_df):
    '''
    given a database of MCMC iteration data, get the gate and error information
    from the best iteration
    '''
    yl,yr,xs=metrics_df.loc[metrics_df[['gateError']].astype('float').idxmin()].iloc[0,2:]
    return yl,yr,xs

def getInterpolatedGates(gate,i):
    '''
    given a series of gates, use the logarithmic spacing function to get more
    granular gates. Without this, the gate will be calculated in linear space and
    not as it appears in logspace, which will create errors in bin score assignment
    high "i" parameter is more accurate but more computationally expensive
    '''
    interpolated_gate = np.zeros((0,2))
    for j in range(len(gate)-1):
        x,y = getLogEquidistantPoints((gate[j][0],gate[j][1]),(gate[j+1][0],gate[j+1][1]),i)
        interpolated_gate = np.concatenate((interpolated_gate,np.hstack((x.reshape(i+1,1),y.reshape(i+1,1)))[:-1]))
    return interpolated_gate

def plotMetrics(metrics_df,sigma=3):
    '''
    given a database of MCMC iterations, plot the error vs iterations. The sigma parameter
    also controls a smoothing parameter that plots a smooth line on top of the
    rough error values.
    '''
    plt.figure()
    if sigma==0:
        plt.plot(metrics_df['gateError'])
    else:
        plt.plot(metrics_df['gateError'],alpha=0.1)
        plt.plot(gaussian_filter1d(metrics_df['gateError'].astype('float'), sigma=sigma))
    plt.xlabel('MCMC Iterations')
    plt.ylabel('FACS error')
    
def plotMinGates(metrics_df,ligand_concs,Kd):
    '''
    given a database of MCMC iterations, plot the gates that gave the minimum error
    '''
    yl,yr,xs = getMinCoordinates(metrics_df)
    gates = convertCoordsToGates(yl, yr, xs)
    plotGates(gates)
    plotFlowData(ligand_concs, gates,Kd)
    plotNegative()
    plt.title('gate error: ' + str(np.around(metrics_df[['gateError']].astype('float').min(0)[0],3)))


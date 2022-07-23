# -*- coding: utf-8 -*-
"""
Created on Wed Feb  9 15:28:49 2022

@author: Marshall
"""

import numpy as np
import matplotlib.pyplot as plt
import matlab.engine
import os

# #define display level
# display_level_median = 5e3
# display_level_variance = 0.3e1
# num_points=1000
# display_level = np.random.lognormal(np.log(display_level_median),np.log(display_level_variance),num_points)

# #define negative population
# neg_display_level_median = 5e1
# neg_display_level_variance = 0.2e1
# neg_num_points = 1000
# neg_display_level = np.random.lognormal(np.log(neg_display_level_median),np.log(neg_display_level_variance),neg_num_points)

# #plot +/- display level distribution
# plt.figure()
# plt.hist(display_level,bins = 10**np.linspace(np.log10(min(display_level)),np.log10(max(display_level)),50),alpha=0.4,label='positive')
# plt.hist(neg_display_level,bins = 10**np.linspace(np.log10(min(neg_display_level)),np.log10(max(neg_display_level)),50),alpha=0.4,label='negative')
# plt.gca().set_xscale('log')
# plt.legend()
# plt.xlabel('Expression Signal')
# plt.show()

# #define binding signal
# noise_median = 7.5e1
# noise_variance = 2.5e1
# fraction_bound = 1
# noise = np.random.lognormal(np.log10(noise_median),np.log10(noise_variance),num_points)
# binding_signal = display_level*np.random.normal(1,0.5,1000)*fraction_bound+noise
# neg_binding_signal = np.random.lognormal(np.log(neg_display_level_median),np.log(neg_display_level_variance),neg_num_points)

# #plot binding vs expression
# plt.figure()
# plt.scatter(display_level,binding_signal,label='frac='+str(fraction_bound),alpha=0.4)
# plt.scatter(neg_display_level,neg_binding_signal,label='neg',alpha=0.4)
# plt.gca().set_xscale('log')
# plt.gca().set_yscale('log')
# plt.xlabel('Expression Signal')
# plt.ylabel('Binding Signal')
# plt.legend()
# plt.show()

# #start and define matlab engine
# eng = matlab.engine.start_matlab()

def getDisplayLevel(n=1000,median=5e3,variance=0.3e1,negative=False):
    if negative:
        median = 2.5e1
        variance = 0.2e1
    display_level = np.random.lognormal(np.log(median),np.log(variance),n)
    return display_level

def getFractionBound(eng,kon1,koff1,kon2,koff2,p1,p2):
    if os.getcwd() != eng.pwd():
        eng.cd(os.getcwd())
    kon1 = float(kon1)
    kon2 = float(kon2)
    return eng.getFracBound(kon1,koff1,kon2,koff2,p1,p2)

def getBindingLevel(display_level,frac=1,n=1000,negative=False,median=1,variance=0.5):
    if negative:
        median = 2.5e1
        variance = 0.2e1
        binding_level = np.random.lognormal(np.log(median),np.log(variance),n)
    else:
        binding_level = display_level*np.random.normal(median,variance,n)*frac
    return binding_level

def getSignal(level,median=1e2,variance=4e1,n=1000):
    signal = level+np.random.normal(median,variance,n)
    return signal

def plotHistogram(signal,label='',ax=None,**kwargs):
    if ax is None:
        ax = plt.gca()
    ax.hist(signal,bins = 10**np.linspace(np.log10(min(signal)),np.log10(max(signal)),50),label=label,**kwargs)
    return ax
    
def plotBindVsExpression(exp,bind,label='',ax=None,alpha=0.4,**kwargs):
    if ax is None:
        ax = plt.gca()
    ax.scatter(exp,bind,label=label,alpha=0.4)
    return ax
    
def plotQuadrants(negative_display,negative_bind):
    x_max = np.percentile(negative_display,99)
    y_max = np.percentile(negative_bind,99)
    plt.plot([x_max,x_max],[1,4e5],color='black')
    plt.plot([1,4e5],[y_max,y_max],color='black')
    
def setLog(label=''):
    plt.gca().set_xscale('log')
    plt.xlabel(str(label)+' Signal')
    plt.legend()
    plt.xlim([10,3e5])
    plt.ylim([10,3e5])
    
def setLogLog(xlabel='Expression',ylabel='Binding'):
    plt.gca().set_xscale('log')
    plt.gca().set_yscale('log')
    plt.xlabel(str(xlabel)+' Signal')
    plt.ylabel(str(ylabel)+' Signal')
    plt.legend()
    plt.xlim([10,3e5])
    plt.ylim([10,3e5])
    
def plot(Kd1=[],Kd2=[],p1=[],p2=[],label=None,kon1=3.17e5,kon2=3.17e5,title=''):
    disp_neg = getDisplayLevel(negative=True)
    bind_neg = getBindingLevel(disp_neg,negative=True)
    disp_neg = getSignal(disp_neg)
    bind_neg = getSignal(bind_neg)
    
    if len(Kd1)+len(Kd2)+len(p1)+len(p2)==4:
        #single plot
        koff1 = Kd1[0]*kon1
        koff2 = Kd2[0]*kon2
        eng = matlab.engine.start_matlab()
        frac = getFractionBound(eng,kon1,koff1,kon2,koff2,p1[0],p2[0])
        
        disp = getDisplayLevel()
        bind = getBindingLevel(disp,frac)
        disp = getSignal(disp)
        bind = getSignal(bind)
        
        if label is None:
            label = [np.around(frac,2)]
        plt.figure()
        plotBindVsExpression(disp_neg, bind_neg,label='negative')
        plotBindVsExpression(disp, bind,label=label[0])
        setLogLog()
        plt.show()
        
    else:
        #plot with kd1,kd2,p1, or p2 vector
        eng = matlab.engine.start_matlab()
        plt.figure()
        plotBindVsExpression(disp_neg, bind_neg,label='negative')
        if len(Kd1) != 1:
            Kd2 = [Kd2[0] for i in range(len(Kd1))]
            p1 = [p1[0] for i in range(len(Kd1))]
            p2 = [p2[0] for i in range(len(Kd1))]
        elif len(Kd2) != 1:
            Kd1 = [Kd1[0] for i in range(len(Kd2))]
            p1 = [p1[0] for i in range(len(Kd2))]
            p2 = [p2[0] for i in range(len(Kd2))]
        elif len(p1) != 1:
            Kd1 = [Kd1[0] for i in range(len(p1))]
            Kd2 = [Kd2[0] for i in range(len(p1))]
            p2 = [p2[0] for i in range(len(p1))]
        elif len(p2) != 1:
            Kd1 = [Kd1[0] for i in range(len(p2))]
            Kd2 = [Kd2[0] for i in range(len(p2))]
            p1 = [p1[0] for i in range(len(p2))]
        index=0
        for Kd_1,Kd_2,p_1,p_2 in zip(Kd1,Kd2,p1,p2):
            koff1 = Kd_1*kon1
            koff2 = Kd_2*kon2
            disp = getDisplayLevel()
            frac = getFractionBound(eng,kon1,koff1,kon2,koff2,p_1,p_2)
            bind = getBindingLevel(disp,frac)
            disp = getSignal(disp)
            bind = getSignal(bind)
            if label is None:
                label = np.around(frac,2)
                plotBindVsExpression(disp, bind,label=label)
                label=None
            else:
                plotBindVsExpression(disp, bind,label=label[index])
            index += 1
        setLogLog()
        plt.title(title)
        plt.show()
            
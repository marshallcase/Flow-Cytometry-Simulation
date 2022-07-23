# -*- coding: utf-8 -*-
"""
Created on Thu Oct 21 12:01:44 2021

@author: marsh
"""

import numpy as np
import pandas as pd
from MCMC_gates import *

# #diagonal gates
# yl = [2e1,5e1,2e2,1.4e3,3e3] #coordinates of left-y gate coordinates
# yr = [2e1,2e3,2e5,1.63e6,4e6] #coordinates of right-y gate coordinates
# xs = [4e2,2e6] #coordinates of x coordinates (left and right)

#sortcery diagonal gates
yl = np.logspace(np.log10(5),np.log10(2.5e3),5) #coordinates of left-y gate coordinates
xs = [4e2,4e6] #coordinates of x coordinates (left and right)
yr = np.logspace(np.log10(5e4),np.log10(2.5e7),5) #coordinates of right-y gate coordinates

# #rectangular gates
# yl = [1e1,1e3,1e4,1e5,1e6]  #coordinates of left-y gate coordinates
# yr = [1e1,1e3,1e4,1e5,1e6] #coordinates of right-y gate coordinates
# xs = [4e2,2e6] #coordinates of x coordinates (left and right)


coords = np.concatenate((yl,yr,xs))
intensity=0.05 #amount of gate coordinate perturbation
fraction = 0.5 #fraction of gates to perturb with each iteration
Kd = 1 #kd of the clone that the gates are tested on (nM)
num_points=1000 #number of clones to simulate at every step
ligand_concs = [0,0.01,0.1,1,10,1000] #range of concentrations to titrate (nM)
gates = convertCoordsToGates(yl,yr,xs)
iterations = 10000 #number of MCMC iterations
metrics_df = pd.DataFrame(index=range(iterations),columns=['gateError','neg_score','yl','yr','xs'])

# for i in range(iterations):
    
#     bin_scores_df = pd.DataFrame(index=ligand_concs,columns=['score'])
#     gates = convertCoordsToGates(yl, yr, xs)
#     for ligand_conc in ligand_concs:
#         exp_signal,bind_signal = getFlowData(ligand_conc,Kd,num_points)
#         bin_score,neg_score = getBinScore(gates,exp_signal,bind_signal)
#         bin_scores_df.loc[ligand_conc,'score']=bin_score
#         bin_scores_df.loc[ligand_conc,'neg_score']=neg_score
        
#     metrics_df.loc[i,['gateError','neg_score']]=getMCMCScore(bin_scores_df,Kd)
#     metrics_df.loc[i,'yl']=yl
#     metrics_df.loc[i,'yr']=yr
#     metrics_df.loc[i,'xs']=xs
#     if i == 0:
#         yl,yr,xs = perturbGates(yl,yr,xs)
#     else:
#         if metrics_df.loc[i,'gateError']<metrics_df.loc[i-1,'gateError']:
#             yl,yr,xs = perturbGates(yl,yr,xs)
#         elif metrics_df.loc[i,'gateError']>=metrics_df.loc[i-1,'gateError']:
#             frac_error = sigmoid(metrics_df.loc[i,'gateError']-metrics_df.loc[i-1,'gateError'])
#             prob_to_perturb = np.random.choice([0,1],p=[1-frac_error,frac_error])
#             if prob_to_perturb == 1:
#                 yl,yr,xs = perturbGates(yl,yr,xs,intensity=0.1)
#             else:
#                 #go back to previous values
#                 yl,yr,xs = metrics_df.loc[i,['yl','yr','xs']]


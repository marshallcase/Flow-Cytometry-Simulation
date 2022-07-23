# -*- coding: utf-8 -*-
"""
Created on Mon Oct 18 12:37:24 2021

@author: marsh
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.path as mplPath

# =============================================================================
# compute bin score of population
# =============================================================================
def getBinScore(gates,exp_signal,bind_signal):
    scores = np.array([])
    for index,gate in enumerate(gates):
        bbPath = mplPath.Path(np.array(gate[:4]))
        score = bbPath.contains_points(np.array([exp_signal,bind_signal]).T).sum()/num_points
        scores = np.append(scores,score)
    bin_score = sum([scores[i]*(i+1) for i in range(len(gates))])
    return bin_score

def plotGates(gates):
    for gate in gates:
        x, y = zip(*gate)
        plt.plot(x,y,color='black')

def getGates(gate_type='titeSeq'):
    if gate_type=='titeSeq':
        ycoords = [2e1,2e2,2e3,2e4,2e5]
        xcoords = [2e1,2e6]
        gates = [[[xcoords[0],ycoords[index]],[xcoords[1],ycoords[index]],[xcoords[1],ycoords[index+1]],
                  [xcoords[0],ycoords[index+1]],[xcoords[0],ycoords[index]]] for index,ycoord in enumerate(ycoords[:4])]
    elif gate_type=='expression':
        ycoords = [2e1,2e2,2e3,2e4,2e5]
        xcoords = [4e2,2e6]
        gates = [[[xcoords[0],ycoords[index]],[xcoords[1],ycoords[index]],[xcoords[1],ycoords[index+1]],
                  [xcoords[0],ycoords[index+1]],[xcoords[0],ycoords[index]]] for index,ycoord in enumerate(ycoords[:4])]
    elif gate_type=='diagonal':
        yl = [2e1,5e1,1e2,6e2,3e3]
        yr = [2e1,2e3,1e5,7e5,4e6]
        xs = [4e2,2e6]
        gates = [[[xs[0],yl[i]],[xs[1],yr[i]],[xs[1],yr[i+1]],[xs[0],yl[i+1]],[xs[0],yl[i]]] for i in range(len(yl)-1)]
    return gates


#parameters
ligand_concs = [0,0.1,1000]
#ligand_concs = np.logspace(-3,3,10)
labels = ['c='+str(ligand_conc) for ligand_conc in ligand_concs]
Kd = 1 #nM
num_points=1000
bin_scores_df = pd.DataFrame(index=ligand_concs,columns=['score'])
plt.figure()

for ligand_conc,label in zip(ligand_concs,labels):
    frac_bound = (ligand_conc)/(ligand_conc+Kd)
    exp_level = np.random.lognormal(np.log(7000),np.log(5),num_points)
    
    exp_signal = exp_level + np.random.normal(100,100,num_points)
    bind_signal = np.random.lognormal(0,0.5,num_points)*exp_level*frac_bound+np.random.normal(100,100,num_points)
    gates = getGates('diagonal')
    bin_score = np.around(getBinScore(gates,exp_signal,bind_signal),2)
    bin_scores_df.loc[ligand_conc,'score']=bin_score
    plt.scatter(exp_signal,bind_signal,label=label+' score='+str(bin_score),alpha=0.4)
    plotGates(gates)
    
#plot negative population
plt.scatter(np.random.normal(100,100,num_points),np.random.normal(100,100,num_points),label='negative',alpha=0.2)
    
plt.yscale('log')
plt.xscale('log')
plt.xlabel('expression')
plt.ylabel('binding')
plt.xlim([10,1e7])
plt.ylim([10,1e7])
plt.legend()
plt.title('titration of Kd=1nM clone')


plt.show()

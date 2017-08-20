# -*- coding: utf-8 -*-
"""
Created on Sat Nov  5 21:43:31 2016

@author: Achin Jain, UPenn

visualization script

"""

import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np
    
def plot_true_predicted(ypred, ytrue, argStr):
    
    sns.set()
    plt.figure()
    l = ytrue.shape[0]
    plt.plot(range(l), ytrue, linewidth=2, label='true')
    plt.plot(range(l), ypred, linewidth=2, label='predicted', linestyle='--')
    plt.xlabel('time index')
    plt.ylabel(argStr)
    plt.show()
    plt.legend(loc=1)
    
def plot_gp(ytrue, ymu, ystd):
    
    sns.set_style('white')
    numSamples = ytrue.shape[0]
    plt.figure(figsize=(11, 8))
    gs = gridspec.GridSpec(3,1)
    plt.subplot(gs[:-1,:])
    xplot = np.linspace(1, numSamples, numSamples)
    plt.plot(xplot, ymu, '#990000', ls='-', lw=1.5, zorder=9, 
             label='GP prediction')
    plt.fill_between(xplot, (ymu+2*ystd), (ymu-2*ystd),
                     alpha=0.2, color='m', label='+-2sigma')
    plt.plot(xplot, ytrue, '#e68a00', ls='-', lw=1, zorder=9, 
             label='ground truth')
    plt.legend(loc='upper right')
    plt.title('true vs predicted')    
    plt.xlim(0,numSamples)          
    plt.subplot(gs[2,:])
    xplot = np.linspace(1, numSamples, numSamples)
    plt.plot(xplot, np.abs(np.array(ytrue).flatten()-ymu), '#990000', 
             ls='-', lw=0.5, zorder=9)
    plt.fill_between(xplot, 0*xplot, 2*ystd,
                     alpha=0.2, color='m')
    plt.title("model error and predicted variance")
    plt.xlim(0,numSamples)
    
    plt.tight_layout()    
    plt.show()
    
def plot_experiment_design(X, y_mean, y_std, X_true, y_true):

    sns.set()
    plt.figure()
    plt.plot(X, y_mean, '#990000', ls='-', lw=1.5, zorder=9, 
             label='GP prediction')
    plt.fill_between(X, (y_mean+2*y_std), (y_mean-2*y_std),
                     alpha=0.2, color='m', label='+-2sigma')
    plt.scatter(X_true, y_true, s = 100, marker = '*', color='g', label='observation')
    plt.legend(loc=0)
    plt.tight_layout()    
    plt.show()
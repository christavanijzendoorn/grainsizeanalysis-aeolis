# -*- coding: utf-8 -*-
"""
Created on Wed Sep 21 11:27:45 2022

@author: cijzendoornvan
"""

import matplotlib.pyplot as plt
import os
import numpy as np

def prep_visualize():
    # Prep visualisation
    S = 14
    M = 18
    L = 20
    
    plt.rc('font', size=S)          # controls default text sizes
    plt.rc('axes', titlesize=S)     # fontsize of the axes title
    plt.rc('axes', labelsize=M)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=S)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=S)    # fontsize of the tick labels
    plt.rc('legend', fontsize=S)    # legend fontsize
    plt.rc('figure', titlesize=L)  # fontsize of the figure title
    
def get_cases(directory):
    cases = []
    for txtfile in os.listdir(directory):
        name, ext = os.path.splitext(txtfile)
        if ext == '.nc':
            cases.append(txtfile)  
    return cases

def calc_ut(d):
    A = 0.085
    rhoa = 1.225
    rhos = 2650
    g = 9.81
    z0    = d / 30.
    kappa = 0.41
    z = 10
    
    ut = A * np.sqrt(d*g*(rhos-rhoa)/rhoa)
    utw = ut / ((kappa) / np.log(z/z0))
    return utw

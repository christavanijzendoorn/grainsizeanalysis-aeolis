# -*- coding: utf-8 -*-
"""
Created on Fri Feb 24 16:38:13 2023

@author: cijzendoornvan
"""

# Load packages
import netCDF4
import numpy as np
import matplotlib.pyplot as plt

import os
import sys
sys.path.append(os.getcwd())
from grainsize_calculation import calc_median_grain_size, calc_mean_grain_size
from visualization_tools import *

prep_visualize() # sets fonts to larger default sizes

model_directory1 = r"C:\Users\cijzendoornvan\Documents\DuneForce\AEOLIS\grainsizeanalysis-aeolis\scenarios\B_singlefraction_5050mix_1day"
# model_directory1 = r"C:\Users\cijzendoornvan\Documents\DuneForce\AEOLIS\grainsizeanalysis-aeolis\scenarios\B_singlefraction_5050mix_1day_Bagnoldgs"

#%% Define visualization function

def fluxes_overview_withwindandflux(model_directory, cases, labels):
    
    colors = ['silver', 'grey', 'k', 'red']
    fig, ax = plt.subplots(3, 1, figsize=(8,7))
    fig.subplots_adjust(left = 0.15, right=0.95, bottom = 0.1, hspace=0.5, wspace = 0.5)
    for i, case in enumerate(cases):
        # Load parameters resulting from model run
        with netCDF4.Dataset(os.path.join(model_directory+'/', case), 'r') as ds:
            t = ds.variables['time'][:]# / 3600.
            qs = ds.variables['qs'][:]
            uw = ds.variables['uw'][:]
            fractions = ds.variables['fractions'][:]
            delta_t = t[1] - t[0]
            if 'day' in case:
                t = t/(3600) # convert time to hours
            elif 'year' in case:
                t = t/(3600*24) # convert time to days            
            if len(fractions) == 1:          
                # Calculate cumulative flux over time at the end of the domain
                qs_cumsum = np.cumsum(qs[:, 0, -2, 0])*delta_t
                ax[1].plot(t, qs[:, 0, -2, 0], color = colors[i], label = labels[i])
            else:
                # Calculate cumulative flux of all fractions
                qs_sum = np.sum(qs[:, 0, :, :], axis =2)
                ax[1].plot(t, qs_sum[:, -2], label = labels[i], color = colors[i], linestyle = 'dotted')
            
                # Calculate cumulative flux over time at the end of the domain
                qs_cumsum = np.cumsum(qs_sum[:,-2])*delta_t
            
            if i == 0:
                qs_cumsum_all = qs_cumsum
            else:
                qs_cumsum_all = np.vstack((qs_cumsum_all, qs_cumsum))
                
    rhog = 2650
    porosity = 0.4
    qs_cumsum_all = qs_cumsum_all/(rhog*(1-porosity))
    for i in range(len(qs_cumsum_all)):
        ax[2].set_ylabel('Cumul. sediment\nflux ($m^3$/m)')
        if i == 3:
            ax[2].plot(t, qs_cumsum_all[i], color = colors[i], label = labels[i], linestyle = 'dotted') 
        else:
            ax[2].plot(t, qs_cumsum_all[i], color = colors[i], label = labels[i], linewidth=2)
    ax[2].set_xlim([0, 24])
    ax[2].set_xlabel('Time (hours)')
    
    ax[1].set_ylabel('Sediment\nflux (kg/m/s)')
    ax[1].legend(ncol = 2)
    ax[1].set_xlim([0, 24])
    
    ax[0].plot(t, uw[:,0,0], linewidth=2)
    ax[0].hlines(calc_ut(0.000250), 0, np.max(t), color = 'silver', label = '250 $\mu$m')
    ax[0].hlines(calc_ut(0.000375), 0, np.max(t), color = 'grey', label = '375 $\mu$m')
    ax[0].hlines(calc_ut(0.000500), 0, np.max(t), color = 'black', label = '500 $\mu$m')
    ax[0].set_ylabel('Wind\nspeed (m/s)')
    ax[0].set_xlim([0, 24])
    ax[0].set_ylim([0, 30])
    
    plt.gcf().text(0.16, 0.9, 'a)', fontsize=16, weight = 'bold')
    plt.gcf().text(0.16, 0.61, 'b)', fontsize=16, weight = 'bold')
    plt.gcf().text(0.16, 0.32, 'c)', fontsize=16, weight = 'bold')
        
    plt.show()
    plt.savefig(model_directory + '/../../analysis/figures/' + 'singlefraction_5050mix_1day.png') #Bagnold_gs
    # plt.close()
    
    return uw

#%% Visualize overview of single fraction and 50-50% mix behavior during 1 day variable wind

cases = get_cases(model_directory1)        

labels = ['250 $\mu$m', '375 $\mu$m', '500 $\mu$m', '250/500\n$\mu$m mix']        

uw = fluxes_overview_withwindandflux(model_directory1, cases, labels)

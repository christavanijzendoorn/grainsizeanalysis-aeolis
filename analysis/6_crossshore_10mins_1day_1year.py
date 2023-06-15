# -*- coding: utf-8 -*-
"""
Created on Fri Feb 24 15:50:30 2023

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
from matplotlib.colors import ListedColormap
from visualization_tools import *

prep_visualize() # sets fonts to larger default sizes

model_directory1 = r"C:\Users\cijzendoornvan\Documents\DuneForce\AEOLIS\grainsizeanalysis-aeolis\scenarios\F_crossshore_10mins"
# model_directory1 = r"C:\Users\cijzendoornvan\Documents\DuneForce\AEOLIS\grainsizeanalysis-aeolis\scenarios\F_crossshore_10mins_Bagnoldgs"
model_directory2 = r"C:\Users\cijzendoornvan\Documents\DuneForce\AEOLIS\grainsizeanalysis-aeolis\scenarios\F_crossshore_1day"
# model_directory2 = r"C:\Users\cijzendoornvan\Documents\DuneForce\AEOLIS\grainsizeanalysis-aeolis\scenarios\F_crossshore_1day_Bagnoldgs"
model_directory3 = r"C:\Users\cijzendoornvan\Documents\DuneForce\AEOLIS\grainsizeanalysis-aeolis\scenarios\F_crossshore_1year"

labels = ['coarse-fine', 'fine-coarse', 'fine-coarse-fine']
colors = ['royalblue', 'forestgreen', 'goldenrod']

#%% Define figure creation

def fluxes_overview_crossshore(model_directory, cases, labels, colors, vmin, vmax):
    
    rhog = 2650
    porosity = 0.4
    with netCDF4.Dataset(os.path.join(model_directory+'/', cases[3]), 'r') as ds:
        qs = ds.variables['qs'][:]
        t = ds.variables['time'][:]# / 3600.
        delta_t = t[1] - t[0]
        qs_cumsum = np.cumsum(qs[:, 0, -2, 0])*delta_t
        qs_normalize = qs_cumsum/(rhog*(1-porosity))
    qs_cumsum = []
    
    fig, ax = plt.subplots(5, 1, figsize=(10,9))
    fig.subplots_adjust(top = 0.95, bottom=0.1, hspace=1.0)
    
    for i, case in enumerate(cases[0:3]):
        # Load parameters resulting from model run
        with netCDF4.Dataset(os.path.join(model_directory+'/', case), 'r') as ds:
            t = ds.variables['time'][:]# / 3600.
            x = ds.variables['x'][...]
            qs = ds.variables['qs'][:]
            uw = ds.variables['uw'][:]
            fractions = ds.variables['fractions'][:]
            mass = ds.variables['mass'][:,:,:,:,:]              
            pickup = ds.variables['pickup'][...]
            delta_t = t[1] - t[0]
            
            
            if len(fractions) == 1:          
                # Calculate cumulative flux over time at the end of the domain
                qs_cumsum = np.cumsum(qs[:, 0, -2, 0])*delta_t
            else:
                # Calculate cumulative flux of all fractions
                qs_sum = np.sum(qs[:, 0, :, :], axis =2)
                # Calculate cumulative flux over time at the end of the domain
                qs_cumsum = np.cumsum(qs_sum[:,-2])*delta_t
            
            if i == 0:
                qs_cumsum_all = qs_cumsum
            else:
                qs_cumsum_all = np.vstack((qs_cumsum_all, qs_cumsum))
                
            if i < 3:
                perc_fine = mass[:,0,:,0,0]/np.sum(mass[:,0,:,0,:], axis = 2)*100
                gs = ax[2+i].pcolormesh(perc_fine[:, :].T, vmax = vmax, vmin = vmin, cmap = 'viridis')
                
                pickup_fine =  pickup[:, 0, :, 0]
                pickup_fine[pickup_fine <= 0] = np.nan
                pickup_fine[pickup_fine > 0] = 1
                pickup_coarse =  pickup[:, 0, :, 1]
                pickup_coarse[pickup_coarse <= 0] = np.nan
                pickup_coarse[pickup_coarse > 0] = 1          
                
                puf = ax[2+i].pcolor(pickup_fine.T, cmap = ListedColormap(['none']), hatch = '\\\\\\', edgecolor='white', lw = 0., label = 'pickup fine') # (t, y, x, frac)
                puc = ax[2+i].pcolor(pickup_coarse.T, cmap = ListedColormap(['none']), hatch = '///', edgecolor='white', lw = 0., label = 'pickup coarse')

                # print(np.nansum(pickup_fine))
                if '1day' in cases[i]:
                    ax[2+i].set_xticks([0, 30, 60, 90, 120])
                    ax[2+i].set_xticklabels([0, 5, 10, 15, 20])
                if '10ms' in cases[i]:
                    ax[2+i].set_xticks([0, 11, 21, 31, 41, 51, 61])
                    ax[2+i].set_xticklabels([0, 100, 200, 300, 400, 500, 600])
                
    ax[2].legend(bbox_to_anchor=(1.38, 1.12))
    ax[2].set_title('Case: coarse - fine')
    ax[3].set_title('Case: fine - coarse')
    ax[4].set_title('Case: fine - coarse - fine')
    ax[3].set_ylabel('Cross-shore location (m)')        
    cbar = fig.colorbar(gs, ax=ax[:], shrink=0.33, location='right', anchor=(0, 0)) 
    cbar.set_label('Percentage fine fraction')

    if '10ms' in cases[0]:
        ax[4].set_xlabel('Time (s)', labelpad=2)    
        ax[0].set_xlim([0,600])
        ax[1].set_xlim([0,600])
        ax[1].set_ylim([-250,250])
        qs_cumsum_all = qs_cumsum_all/(rhog*(1-porosity))
        ax[1].set_ylabel('Diff. in cumul.\nsed. flux (%)')
    elif '1day' in cases[0]:
        t = t/3600
        qs_cumsum_all = qs_cumsum_all/(rhog*(1-porosity))
        ax[1].set_ylabel('Diff. in cumul.\nsed. flux (%)')
        ax[4].set_xlabel('Time (hours)')
        ax[0].set_xlim([0,24])
        ax[1].set_xlim([0,24])
        ax[1].set_ylim([-22,22])
    elif '1year' in cases[0]: 
        t = t/(3600*24)
        qs_cumsum_all = qs_cumsum_all/(rhog*(1-porosity))
        ax[1].set_ylabel('Diff. in cumul.\nsed. flux (%)')
        ax[4].set_xlabel('Time (days)')
        ax[0].set_xlim([0,365])
        ax[1].set_xlim([0,365])
        ax[1].set_ylim([-22,22])

    if len(qs_cumsum_all.shape) == 1:
        ax[1].plot(t, (qs_cumsum_all-qs_normalize)/np.max(qs_normalize)*100, label = labels[i], linewidth=2, color = colors[i])
    else:
        for i in range(len(qs_cumsum_all)):
            ax[1].plot(t, (qs_cumsum_all[i]-qs_normalize)/np.max(qs_normalize)*100, label = labels[i], linewidth=2, color = colors[i])
             
    ax[1].plot(t, qs_normalize-qs_normalize, label = 'reference D$_{50}$', linewidth=1.5, color = 'k')
    ax[1].legend(bbox_to_anchor=(1.02, 1.12))
    ax[0].plot(t, uw[:,0,0])
    ax[0].set_ylabel('Wind (m/s)')
    
    plt.gcf().text(0.14, 0.965, 'a)', fontsize=16, weight = 'bold')
    plt.gcf().text(0.14, 0.78, 'b)', fontsize=16, weight = 'bold')
    plt.gcf().text(0.14, 0.59, 'c)', fontsize=16, weight = 'bold')
    plt.gcf().text(0.14, 0.4, 'd)', fontsize=16, weight = 'bold')
    plt.gcf().text(0.14, 0.215, 'e)', fontsize=16, weight = 'bold')
    
    plt.show()
    plt.savefig(model_directory + '/../../analysis/figures/' + model_directory[86:] + '.png') #_Bagnoldgs

#%% Visualization on 10-minute scale

cases = get_cases(model_directory1)
        
fluxes_overview_crossshore(model_directory1, cases, labels, colors, vmin = 0, vmax = 100)

#%% Visualization on 1 day scale
        
cases = get_cases(model_directory2)        
        
fluxes_overview_crossshore(model_directory2, cases, labels, colors, vmin = 0, vmax = 100)

#%% Visualization on 1 year scale

cases = get_cases(model_directory3)       
        
fluxes_overview_crossshore(model_directory3, cases, labels, colors, vmin = 0, vmax = 100)
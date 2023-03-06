# -*- coding: utf-8 -*-
"""
Created on Fri Feb 24 16:54:50 2023

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

model_directory1 = r"C:\Users\cijzendoornvan\Documents\DuneForce\AEOLIS\GrainSize\scenarios\C_8020mixes_1year"

#%% Define visualization function

def fluxes_overview_1yr_D50mix(model_directory, cases, labels, colors, vmin, vmax):
    
    rhog = 2650
    porosity = 0.4
    
    with netCDF4.Dataset(os.path.join(model_directory+'/', cases[0]), 'r') as ds:
        pickup_mix = ds.variables['Cu'][...]       
        x = ds.variables['x'][...]       
    with netCDF4.Dataset(os.path.join(model_directory+'/', cases[3]), 'r') as ds:
        pickup_single = ds.variables['Cu'][...]       

    with netCDF4.Dataset(os.path.join(model_directory+'/', cases[-2]), 'r') as ds:
        qs = ds.variables['qs'][:]
        t = ds.variables['time'][:]# / 3600.
        delta_t = t[1] - t[0]
        qs_cumsum = np.cumsum(qs[:, 0, -2, 0])*delta_t
        qs_normalize = qs_cumsum/(rhog*(1-porosity))
    qs_cumsum = []
    
    fig, ax = plt.subplot_mosaic([['top', 'top', 'top', 'blank'], ['middle', 'middle', 'middle', 'blank'], ['bottom', 'bottom', 'bottom', 'blank']], empty_sentinel='blank',
                              constrained_layout=True, figsize=(10,8))
    fig.subplots_adjust(left = 0.12, right=0.95, bottom = 0.1, top = 0.95, hspace = 0.5)
        
    for i, case in enumerate(cases):
        # Load parameters resulting from model run
        with netCDF4.Dataset(os.path.join(model_directory+'/', case), 'r') as ds:
            t = ds.variables['time'][:]# / 3600.
            qs = ds.variables['qs'][:]
            uw = ds.variables['uw'][:]
            fractions = ds.variables['fractions'][:]
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
                
            # Include plot of grain size in bed surface    
            if '125_2000' in case:
                mass = ds.variables['mass'][:,:,:,:,:]
                perc_fine = mass[:,0,:,0,0]/np.sum(mass[:,0,:,0,:], axis = 2)*100
                gs = ax['bottom'].pcolormesh(perc_fine[:, :].T, vmax = vmax, vmin = vmin, cmap = 'viridis')
                
                ax['bottom'].set_title('Case: 125/2000 $\mu$m')
                ax['bottom'].set_ylabel('Cross-shore\nlocation (m)')
                ax['bottom'].set_xlim([0,365])
                
                from mpl_toolkits.axes_grid1.inset_locator import inset_axes
                axins_cbar = ax['bottom'].inset_axes([1.05, 0, 0.05, 1])
                cbar = fig.colorbar(gs, cax=axins_cbar)
                cbar.set_label('Percentage\nfine fraction')
            
    if 'day' in case:
        t = t/(3600) # convert time to hours
        ax[0].set_ylabel('Sediment flux (kg/m)')
        ax[2].set_xlabel('Time (hours)')
    elif 'year' in case:
        t = t/(3600*24) # convert time to days       
        ax['middle'].set_ylabel('Difference in cumul.\nsediment flux (%)')
        ax['bottom'].set_xlabel('Time (days)')
        
    for i in range(len(qs_cumsum_all)):
         # convert to m3/m
        qs_cumsum_single = qs_cumsum_all[i,:]/(rhog*(1-porosity))
        if i<3:
            ax['middle'].plot(t, (qs_cumsum_single-qs_normalize)/np.max(qs_normalize)*100, color=colors[i], linewidth=2.5, label = labels[i])    
        else:    
            ax['middle'].plot(t, (qs_cumsum_single-qs_normalize)/np.max(qs_normalize)*100, color=colors[i], linewidth=1, label = labels[i])
        
        ax['middle'].legend(loc = (1.05, -0.05)) # loc = 'upper left' bbox_to_anchor=(1.0, 1.03)
        ax['middle'].set_xlim([0,365])
        ax['middle'].set_ylim([-20,40])

    grainsizes = [0.000125, 0.000500, 0.002000]
    colors = ['silver', 'grey', 'k']
    for i in range(0, len(grainsizes)):
        ax['top'].hlines(calc_ut(grainsizes[i]), 0, 365, linewidth = 1.5, color = colors[i], label = str(int(grainsizes[i]*1e6)) + ' $\mu$m')
    ax['top'].plot(t, uw[:,0,0], linewidth=1.5)
    ax['top'].legend(loc = (1.05, 0.25), title = 'Threshold velocity') #bbox_to_anchor=(1, 1.03), , bbox_to_anchor=(0.90, 0.4, 1.0, 0.4)
    
    ax['top'].set_ylabel('Wind\nspeed (m/s)')
    ax['top'].set_xlim([0,365])
    ax['top'].set_ylim([0,31])
    
    plt.gcf().text(0.13, 0.97, 'a)', fontsize=16, weight = 'bold')
    plt.gcf().text(0.13, 0.65, 'b)', fontsize=16, weight = 'bold')
    plt.gcf().text(0.13, 0.33, 'c)', fontsize=16, weight = 'bold')

    plt.show()
    plt.savefig(model_directory + '/../../analysis/figures/' + '8020mixes_1year.png')    

#%% Go through all files in directory and visualize for each run
cases = get_cases(model_directory1)       

labels = ['125/2000 $\mu$m', '300/1300 $\mu$m', '375/1000 $\mu$m', '125 $\mu$m', '500 $\mu$m', '2000 $\mu$m']
colors = ['royalblue', 'forestgreen', 'goldenrod', 'silver', 'grey', 'k']
      
fluxes_overview_1yr_D50mix(model_directory1, cases, labels, colors, vmin=0, vmax=100)


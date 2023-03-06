# -*- coding: utf-8 -*-
"""
Created on Mon Jul 18 14:46:15 2022

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

model_directory1 = r"C:\Users\cijzendoornvan\Documents\DuneForce\AEOLIS\GrainSize\scenarios\D_variedpercentagemixes_1year"

#%% Define visualizetion function

def fluxes_overview_1yr_percmix(model_directory, cases, labels, colors, vmin, vmax):
    
    rhog = 2650
    porosity = 0.4
    
    with netCDF4.Dataset(os.path.join(model_directory+'/', cases[-1]), 'r') as ds:
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
            if '90-10' in case:
                mass = ds.variables['mass'][:,:,:,:,:]
                pickup = ds.variables['Cu'][...]  
                
                pickup_fine =  pickup[:, 0, :, 0]
                pickup_fine[pickup_fine <= 0] = np.nan
                pickup_fine[pickup_fine > 0] = 1
                pickup_coarse =  pickup[:, 0, :, 1]
                pickup_coarse[pickup_coarse <= 0] = np.nan
                pickup_coarse[pickup_coarse > 0] = 1
                
                perc_fine = mass[:,0,:,0,0]/np.sum(mass[:,0,:,0,:], axis = 2)*100
                gs = ax['bottom'].pcolormesh(perc_fine[:, :].T, vmax = vmax, vmin = vmin, cmap = 'viridis')

                ax['bottom'].set_title('Case: 90% - 10%')
                ax['bottom'].set_ylabel('Cross-shore\nlocation (m)')
                ax['bottom'].set_xlim([0,365])
                
                axins_cbar = ax['bottom'].inset_axes([1.05, 0, 0.05, 1])
                cbar = fig.colorbar(gs, cax=axins_cbar)
                cbar.set_label('Percentage\nfine fraction')
            
    if 'day' in case:
        t = t/(3600) # convert time to hours
        ax[0].set_ylabel('Sediment flux (kg/m)')
        ax[1].set_xlabel('Time (hours)')
    elif 'year' in case:
        t = t/(3600*24) # convert time to days       
        ax['middle'].set_ylabel('Difference in cumul.\nsediment flux (%)')
        ax['bottom'].set_xlabel('Time (days)')
        
    for i in range(len(qs_cumsum_all)):
         # convert to m3/m
        qs_cumsum_single = qs_cumsum_all[i,:]/(rhog*(1-porosity))
        if i<4:
            ax['middle'].plot(t, (qs_cumsum_single-qs_normalize)/np.max(qs_normalize)*100, color=colors[i], linewidth=2.5, label = labels[i])    
        else:    
            ax['middle'].plot(t, (qs_cumsum_single-qs_normalize)/np.max(qs_normalize)*100, color=colors[i], linewidth=1, label = labels[i])
        
        ax['middle'].legend(loc = (1.05, 0.1)) # loc = 'upper left' bbox_to_anchor=(1.0, 1.03)
        ax['middle'].set_xlim([0,365])
        ax['middle'].set_ylim([-50,6])

    grainsizes = [0.000125, 0.002000]
    colors = ['silver', 'k']
    for i in range(0, len(grainsizes)):
        ax['top'].hlines(calc_ut(grainsizes[i]), 0, 365, linewidth = 1.5, color = colors[i], label = str(int(grainsizes[i]*1e6)) + ' $\mu$m')
    ax['top'].plot(t, uw[:,0,0], linewidth=1.5)
    ax['top'].legend(loc = (1.05, 0.4), title = 'Threshold velocity') #bbox_to_anchor=(1, 1.03), , bbox_to_anchor=(0.90, 0.4, 1.0, 0.4)
    
    ax['top'].set_ylabel('Wind\nspeed (m/s)')
    ax['top'].set_xlim([0,365])
    ax['top'].set_ylim([0,31])
    
    plt.gcf().text(0.13, 0.97, 'a)', fontsize=16, weight = 'bold')
    plt.gcf().text(0.13, 0.65, 'b)', fontsize=16, weight = 'bold')
    plt.gcf().text(0.13, 0.33, 'c)', fontsize=16, weight = 'bold')
    
    plt.show()
    plt.savefig(model_directory + '/../../analysis/figures/' + 'variedpercentagemixes_1year.png')      

#%% Go through all files in directory and visualize for each run
cases = get_cases(model_directory1)      

labels = ['80% - 20%', '90% - 10%', '95% - 5%', '99% - 1%', '125 $\mu$m']
colors = ['goldenrod', 'forestgreen', 'royalblue', 'grey', 'silver']
      
fluxes_overview_1yr_percmix(model_directory1, cases, labels, colors, vmin=0, vmax=100)


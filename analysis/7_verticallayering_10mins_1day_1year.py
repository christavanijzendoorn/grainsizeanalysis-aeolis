# -*- coding: utf-8 -*-
"""
Created on Fri Feb 24 16:02:35 2023

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

model_directory1 = r"C:\Users\cijzendoornvan\Documents\DuneForce\AEOLIS\grainsizeanalysis-aeolis\scenarios\M_verticallayering_10mins"
# model_directory2 = r"C:\Users\cijzendoornvan\Documents\DuneForce\AEOLIS\grainsizeanalysis-aeolis\scenarios\M_verticallayering_1day"
model_directory2 = r"C:\Users\cijzendoornvan\Documents\DuneForce\AEOLIS\grainsizeanalysis-aeolis\scenarios\M_verticallayering_1day_Bagnoldgs"
model_directory3 = r"C:\Users\cijzendoornvan\Documents\DuneForce\AEOLIS\grainsizeanalysis-aeolis\scenarios\M_verticallayering_1year"

labels = ['1$^{st}$ coarse', '4$^{th}$ coarse']
colors = ['royalblue', 'forestgreen']

#%% Define visualization functions

def verticallayering_10mins(model_directory, cases):

    rhog = 2650
    porosity = 0.4
    with netCDF4.Dataset(os.path.join(model_directory+'/', cases[3]), 'r') as ds:
        qs = ds.variables['qs'][:]
        t = ds.variables['time'][:]# / 3600.
        delta_t = t[1] - t[0]
        qs_cumsum = np.cumsum(qs[:, 0, -2, 0])*delta_t
        qs_normalize_300_10ms = qs_cumsum/(rhog*(1-porosity))
    with netCDF4.Dataset(os.path.join(model_directory+'/', cases[4]), 'r') as ds:
        qs = ds.variables['qs'][:]
        t = ds.variables['time'][:]# / 3600.
        delta_t = t[1] - t[0]
        qs_cumsum = np.cumsum(qs[:, 0, -2, 0])*delta_t
        qs_normalize_300_20ms = qs_cumsum/(rhog*(1-porosity))
    with netCDF4.Dataset(os.path.join(model_directory+'/', cases[5]), 'r') as ds:
        qs = ds.variables['qs'][:]
        t = ds.variables['time'][:]# / 3600.
        delta_t = t[1] - t[0]
        qs_cumsum = np.cumsum(qs[:, 0, -2, 0])*delta_t
        qs_normalize_300_30ms = qs_cumsum/(rhog*(1-porosity))        

    fig, axs = plt.subplots(1,3, figsize=(12,4), sharex=True)
    fig.subplots_adjust(top=0.8)
    fig.subplots_adjust(bottom=0.2)
    fig.subplots_adjust(left=0.15)
    fig.subplots_adjust(right=0.84)
    fig.subplots_adjust(wspace=0.28)
    
    labels = ['250 $\mu$m', '250 $\mu$m', '250 $\mu$m', '300 $\mu$m', '300 $\mu$m', '300 $\mu$m', '500 $\mu$m', '500 $\mu$m', '500 $\mu$m', '1$^{st}$ coarse', '1$^{st}$ coarse', '1$^{st}$ coarse', '4$^{th}$ coarse', '4$^{th}$ coarse', '4$^{th}$ coarse']
    colors = ['lightgrey', 'lightgrey', 'lightgrey', 'grey', 'grey', 'grey', 'k', 'k', 'k', 'royalblue', 'royalblue', 'royalblue', 'forestgreen', 'forestgreen', 'forestgreen']
    plotindex = [0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2]
    
    for i, case in enumerate(cases):    
    
        # Load parameters resulting from model run
        with netCDF4.Dataset(os.path.join(model_directory+'/', case), 'r') as ds:
            t = ds.variables['time'][:]# / 3600.
            qs = ds.variables['qs'][:]
            delta_t = t[1] - t[0]

            qs = qs/(rhog*(1-porosity))
            qs_cumsum = np.cumsum(qs[:, 0, -2, 0])*delta_t

            axs[plotindex[i]].tick_params(axis='both', which='major', labelsize=16)
            axs[plotindex[i]].yaxis.offsetText.set_fontsize(16)
            
        if '250' in case:
            qs_cumsum = np.cumsum(qs[:, 0, -2, 0])*delta_t
            if '10ms' in case:
                axs[plotindex[i]].plot(t[:], (qs_cumsum-qs_normalize_300_10ms)/np.max(qs_normalize_300_10ms)*100, label=labels[i], linewidth=2.0, color=colors[i]) 
            elif '20ms' in case:
                axs[plotindex[i]].plot(t[:], (qs_cumsum-qs_normalize_300_20ms)/np.max(qs_normalize_300_20ms)*100, label=labels[i], linewidth=2.0, color=colors[i]) 
            elif '30ms' in case:
                axs[plotindex[i]].plot(t[:], (qs_cumsum-qs_normalize_300_30ms)/np.max(qs_normalize_300_30ms)*100, label=labels[i], linewidth=2.0, color=colors[i]) 
        elif '300' in case:
            qs_cumsum = np.cumsum(qs[:, 0, -2, 0])*delta_t
            if '10ms' in case:
                axs[plotindex[i]].plot(t[:], (qs_cumsum-qs_normalize_300_10ms)/np.max(qs_normalize_300_10ms)*100, label=labels[i], linewidth=2.0, color=colors[i]) 
            elif '20ms' in case:
                axs[plotindex[i]].plot(t[:], (qs_cumsum-qs_normalize_300_20ms)/np.max(qs_normalize_300_20ms)*100, label=labels[i], linewidth=2.0, color=colors[i]) 
            elif '30ms' in case:
                axs[plotindex[i]].plot(t[:], (qs_cumsum-qs_normalize_300_30ms)/np.max(qs_normalize_300_30ms)*100, label=labels[i], linewidth=2.0, color=colors[i]) 
        elif '500' in case:
            qs_cumsum = np.cumsum(qs[:, 0, -2, 0])*delta_t
            if '10ms' in case:
                axs[plotindex[i]].plot(t[:], (qs_cumsum-qs_normalize_300_10ms)/np.max(qs_normalize_300_10ms)*100, label=labels[i], linewidth=2.0, color=colors[i]) 
            elif '20ms' in case:
                axs[plotindex[i]].plot(t[:], (qs_cumsum-qs_normalize_300_20ms)/np.max(qs_normalize_300_20ms)*100, label=labels[i], linewidth=2.0, color=colors[i]) 
            elif '30ms' in case:
                axs[plotindex[i]].plot(t[:], (qs_cumsum-qs_normalize_300_30ms)/np.max(qs_normalize_300_30ms)*100, label=labels[i], linewidth=2.0, color=colors[i]) 
        elif '10ms' in case:
            axs[plotindex[i]].set_title('a) 10 m/s',fontweight='bold')
            qs_cumsum = np.cumsum(qs[:, 0, -2, 0] + qs[:, 0, -2, 1])*delta_t
            axs[plotindex[i]].plot(t[:], (qs_cumsum-qs_normalize_300_10ms)/np.max(qs_normalize_300_10ms)*100, label=labels[i], linewidth=2.0, color=colors[i]) 
            axs[plotindex[i]].set_xlim([0,600])
            axs[plotindex[i]].set_ylim([-100,75])
        elif '20ms' in case:
            axs[plotindex[i]].set_title('b) 20 m/s',fontweight='bold')
            qs_cumsum = np.cumsum(qs[:, 0, -2, 0] + qs[:, 0, -2, 1])*delta_t
            axs[plotindex[i]].plot(t[:], (qs_cumsum-qs_normalize_300_20ms)/np.max(qs_normalize_300_20ms)*100, label=labels[i], linewidth=2.0, color=colors[i]) 
            axs[plotindex[i]].set_xlim([0,600])
            axs[plotindex[i]].set_ylim([-30,22.5])
        elif '30ms' in case:
            axs[plotindex[i]].set_title('c) 30 m/s',fontweight='bold')
            qs_cumsum = np.cumsum(qs[:, 0, -2, 0] + qs[:, 0, -2, 1])*delta_t
            axs[plotindex[i]].plot(t[:], (qs_cumsum-qs_normalize_300_30ms)/np.max(qs_normalize_300_30ms)*100, label=labels[i], linewidth=2.0, color=colors[i]) 
            axs[plotindex[i]].set_xlim([0,600])
            axs[plotindex[i]].set_ylim([-15,11.25])
        
    axs[2].legend(bbox_to_anchor=(1.03, 1.05))
    axs[1].set_xlabel('Time (s)', fontsize=18)
    axs[0].set_ylabel('Difference in cumul.\nsediment flux (%)', fontsize=18, labelpad=25)

    plt.savefig(model_directory + '/../../analysis/figures/' + 'verticallayering_10mins.png') 
    
def fluxes_overview_verticallayering(model_directory, cases, labels, colors, vmin, vmax):
    
    rhog = 2650
    porosity = 0.4
    with netCDF4.Dataset(os.path.join(model_directory+'/', cases[2]), 'r') as ds:
        qs = ds.variables['qs'][:]
        t = ds.variables['time'][:]
        delta_t = t[1] - t[0]
        qs_cumsum = np.cumsum(qs[:, 0, -2, 0])*delta_t
        qs_normalize = qs_cumsum/(rhog*(1-porosity))
    qs_cumsum = []    
    
    fig, ax = plt.subplots(4, 1, figsize=(10,9))
    fig.subplots_adjust(bottom = 0.1, top = 0.95, hspace = 0.7)
    
    for i, case in enumerate(cases[0:2]):
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
                
                
            #d50, sorting and skewness
            if i < 2:
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
                if '1day' in cases[i]:
                    ax[2+i].set_xticks([0, 30, 60, 90, 120])
                    ax[2+i].set_xticklabels([0, 5, 10, 15, 20])
            
    ax[2].legend(bbox_to_anchor=(1.38, 1.12))
    ax[2].set_title('1$^{st}$ layer coarse') 
    ax[3].set_title('4$^{th}$ layer coarse')
    fig.text(0.05, 0.28, 'Cross-shore location (m)', va='center', ha='center', rotation='vertical', fontsize=18)       
     
    cbar = fig.colorbar(gs, ax=ax[:], shrink=0.33, location='right', anchor=(0, -0.045)) 
    cbar.set_label('Percentage\nfine fraction')
    
    if '1day' in cases[0]:  
        t = t/3600
        qs_cumsum_all = qs_cumsum_all/(rhog*(1-porosity))
        ax[1].set_ylabel('Diff. in cumul.\nsed. flux (%)')
        ax[3].set_xlabel('Time (hours)')
        ax[0].set_xlim([0,24])
        ax[1].set_xlim([0,24])
        ax[1].set_ylim([-8,8])
    elif '1year' in cases[0]:
        t = t/(3600*24)
        qs_cumsum_all = qs_cumsum_all/(rhog*(1-porosity))
        ax[1].set_ylabel('Diff. in cumul.\nsed. flux (%)')
        ax[3].set_xlabel('Time (days)')
        ax[0].set_xlim([0,365])
        ax[1].set_xlim([0,365])
        ax[1].set_ylim([-1.5,1.5])

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
    plt.gcf().text(0.14, 0.73, 'b)', fontsize=16, weight = 'bold')
    plt.gcf().text(0.14, 0.495, 'c)', fontsize=16, weight = 'bold')
    plt.gcf().text(0.14, 0.255, 'd)', fontsize=16, weight = 'bold')
    
    plt.show()
    plt.savefig(model_directory + '/../../analysis/figures/' + model_directory[88:] + '.png')     

#%% Visualization on 10-minute scale
cases = get_cases(model_directory1)

verticallayering_10mins(model_directory1, cases)    

#%% Visualization on 1 day scale
cases = get_cases(model_directory2)

fluxes_overview_verticallayering(model_directory2, cases, labels, colors, vmin = 0, vmax = 100) 

#%% Visualization on 1 year scale
cases = get_cases(model_directory3)
      
fluxes_overview_verticallayering(model_directory3, cases, labels, colors, vmin = 0, vmax = 100)

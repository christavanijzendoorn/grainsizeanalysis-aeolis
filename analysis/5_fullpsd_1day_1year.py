# -*- coding: utf-8 -*-
"""
Created on Fri Feb 24 16:42:47 2023

@author: cijzendoornvan
"""

# Load packages
import netCDF4
import numpy as np
import matplotlib.pyplot as plt

import os
import sys
sys.path.append(os.getcwd())
from visualization_tools import *

prep_visualize() # sets fonts to larger default sizes

model_directory1 = r"C:\Users\cijzendoornvan\Documents\DuneForce\AEOLIS\GrainSize\scenarios\E_fullPSDs_1day"
model_directory2 = r"C:\Users\cijzendoornvan\Documents\DuneForce\AEOLIS\GrainSize\scenarios\E_fullPSDs_1year"

#%% Define visualization functions

def fluxes_overview_1day_shapesnorm_detailed(model_directory, cases, labels, colors):
    rhog = 2650
    porosity = 0.4
    with netCDF4.Dataset(os.path.join(model_directory+'/', cases[0]), 'r') as ds:
        qs = ds.variables['qs'][:]
        t = ds.variables['time'][:]# / 3600.
        delta_t = t[1] - t[0]
        qs_cumsum = np.cumsum(qs[:, 0, -2, 0])*delta_t
        qs_normalize230 = qs_cumsum/(rhog*(1-porosity))
    with netCDF4.Dataset(os.path.join(model_directory+'/', cases[1]), 'r') as ds:
        qs = ds.variables['qs'][:]
        t = ds.variables['time'][:]# / 3600.
        delta_t = t[1] - t[0]
        qs_cumsum = np.cumsum(qs[:, 0, -2, 0])*delta_t
        qs_normalize240 = qs_cumsum/(rhog*(1-porosity))
    with netCDF4.Dataset(os.path.join(model_directory+'/', cases[2]), 'r') as ds:
        qs = ds.variables['qs'][:]
        t = ds.variables['time'][:]# / 3600.
        delta_t = t[1] - t[0]
        qs_cumsum = np.cumsum(qs[:, 0, -2, 0])*delta_t
        qs_normalize250 = qs_cumsum/(rhog*(1-porosity))        
        
    with netCDF4.Dataset(os.path.join(model_directory+'/', cases[3]), 'r') as ds:
        qs = ds.variables['qs'][:]
        t = ds.variables['time'][:]# / 3600.
        delta_t = t[1] - t[0]
        qs_cumsum = np.cumsum(qs[:, 0, -2, 0])*delta_t
        qs_normalize470 = qs_cumsum/(rhog*(1-porosity))
    with netCDF4.Dataset(os.path.join(model_directory+'/', cases[4]), 'r') as ds:
        qs = ds.variables['qs'][:]
        t = ds.variables['time'][:]# / 3600.
        delta_t = t[1] - t[0]
        qs_cumsum = np.cumsum(qs[:, 0, -2, 0])*delta_t
        qs_normalize480 = qs_cumsum/(rhog*(1-porosity))
    with netCDF4.Dataset(os.path.join(model_directory+'/', cases[5]), 'r') as ds:
        qs = ds.variables['qs'][:]
        t = ds.variables['time'][:]# / 3600.
        delta_t = t[1] - t[0]
        qs_cumsum = np.cumsum(qs[:, 0, -2, 0])*delta_t
        qs_normalize490 = qs_cumsum/(rhog*(1-porosity))
    qs_cumsum = []
    
    fig, ax = plt.subplots(3, 1, figsize=(10,6))
    fig.subplots_adjust(left = 0.14, right=0.75, bottom = 0.1, hspace=0.5, wspace = 0.5)
    for i, case in enumerate(cases):
        # Load parameters resulting from model run
        with netCDF4.Dataset(os.path.join(model_directory+'/', case), 'r') as ds:
            t = ds.variables['time'][:]# / 3600.
            qs = ds.variables['qs'][:]
            uw = ds.variables['uw'][:]
            fractions = ds.variables['fractions'][:]
            delta_t = t[1] - t[0]
            mass = ds.variables['mass'][:,:,:,:,:]
            # print(np.sum(mass[0, 0, :, 0]))
            
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
                # print(qs_cumsum_all.shape, qs_cumsum.shape)
                qs_cumsum_all = np.vstack((qs_cumsum_all, qs_cumsum))
    if 'day' in case:
        t = t/(3600) # convert time to hours
        # ax[1].set_ylabel('Difference in cumul.\nsediment flux (%)')
        ax[2].set_xlabel('Time (hours)')
        ax[0].set_xlim([0,24])
        ax[1].set_xlim([0,24])
        ax[2].set_xlim([0,24])
        ax[1].set_ylim([-15,3])
        ax[2].set_ylim([-15,3])
    for i in range(len(qs_cumsum_all)):
        qs_cumsum_single = qs_cumsum_all[i,:]/(rhog*(1-porosity))
        if i == 9:
            ax[1].plot(t, (qs_cumsum_single-qs_normalize230)/np.max(qs_normalize230)*100, color=colors[i], linewidth=2, label = labels[i])
        elif i == 10:
            ax[1].plot(t, (qs_cumsum_single-qs_normalize240)/np.max(qs_normalize240)*100, color=colors[i], linewidth=2, label = labels[i])
        elif i == 11:
            ax[1].plot(t, (qs_cumsum_single-qs_normalize250)/np.max(qs_normalize250)*100, color=colors[i], linewidth=2, label = labels[i])
        elif i == 6:
            ax[2].plot(t, (qs_cumsum_single-qs_normalize470)/np.max(qs_normalize470)*100, color=colors[i], linewidth=2, label = labels[i])
        elif i == 7:
            ax[2].plot(t, (qs_cumsum_single-qs_normalize480)/np.max(qs_normalize480)*100, color=colors[i], linewidth=2, label = labels[i])
        elif i == 8:
            ax[2].plot(t, (qs_cumsum_single-qs_normalize490)/np.max(qs_normalize490)*100, color=colors[i], linewidth=2, label = labels[i])            
        elif i == 0: 
            ax[1].plot(t, (qs_cumsum_single-qs_normalize230)/np.max(qs_normalize230)*100, color=colors[i], linewidth=2, label = labels[i])
        elif i == 3: 
            ax[2].plot(t, (qs_cumsum_single-qs_normalize470)/np.max(qs_normalize470)*100, color=colors[i], linewidth=2, label = labels[i])
        ax[1].legend(bbox_to_anchor=(1.03, 1.09))
        ax[2].legend(bbox_to_anchor=(1.03, 1.09))
        
    fig.text(0.05, 0.33, 'Difference in cumul.\nsediment flux (%)', va='center', ha='center', rotation='vertical', fontsize=18)
        
    ax[0].plot(t, uw[:,0,0], linewidth=1.5)
    ax[0].hlines(calc_ut(0.000250), 0, 365, color = 'grey', label = '250 $\mu$m' )
    ax[0].hlines(calc_ut(0.000500), 0, 365, color = 'k', label = '500 $\mu$m')
    ax[0].set_ylabel('Wind\nspeed (m/s)')
    ax[0].legend(loc = (1.045, 0.18), title = 'Threshold velocity') 
    
    plt.gcf().text(0.15, 0.89, 'a)', fontsize=16, weight = 'bold')
    plt.gcf().text(0.15, 0.60, 'b)', fontsize=16, weight = 'bold')
    plt.gcf().text(0.15, 0.31, 'c)', fontsize=16, weight = 'bold')
    
    plt.show()
    plt.savefig(model_directory + '/../../analysis/figures/' + model_directory[73:] + '.png')         
    
def fluxes_overview_1year_shapesnorm_detailed(model_directory, cases, labels, colors):
    rhog = 2650
    porosity = 0.4
    with netCDF4.Dataset(os.path.join(model_directory+'/', cases[0]), 'r') as ds:
        qs = ds.variables['qs'][:]
        t = ds.variables['time'][:]# / 3600.
        delta_t = t[1] - t[0]
        qs_cumsum = np.cumsum(qs[:, 0, -2, 0])*delta_t
        qs_normalize230 = qs_cumsum/(rhog*(1-porosity))
    with netCDF4.Dataset(os.path.join(model_directory+'/', cases[1]), 'r') as ds:
        qs = ds.variables['qs'][:]
        t = ds.variables['time'][:]# / 3600.
        delta_t = t[1] - t[0]
        qs_cumsum = np.cumsum(qs[:, 0, -2, 0])*delta_t
        qs_normalize240 = qs_cumsum/(rhog*(1-porosity))
    with netCDF4.Dataset(os.path.join(model_directory+'/', cases[2]), 'r') as ds:
        qs = ds.variables['qs'][:]
        t = ds.variables['time'][:]# / 3600.
        delta_t = t[1] - t[0]
        qs_cumsum = np.cumsum(qs[:, 0, -2, 0])*delta_t
        qs_normalize250 = qs_cumsum/(rhog*(1-porosity))        
        
    with netCDF4.Dataset(os.path.join(model_directory+'/', cases[3]), 'r') as ds:
        qs = ds.variables['qs'][:]
        t = ds.variables['time'][:]# / 3600.
        delta_t = t[1] - t[0]
        qs_cumsum = np.cumsum(qs[:, 0, -2, 0])*delta_t
        qs_normalize470 = qs_cumsum/(rhog*(1-porosity))
    with netCDF4.Dataset(os.path.join(model_directory+'/', cases[4]), 'r') as ds:
        qs = ds.variables['qs'][:]
        t = ds.variables['time'][:]# / 3600.
        delta_t = t[1] - t[0]
        qs_cumsum = np.cumsum(qs[:, 0, -2, 0])*delta_t
        qs_normalize480 = qs_cumsum/(rhog*(1-porosity))
    with netCDF4.Dataset(os.path.join(model_directory+'/', cases[5]), 'r') as ds:
        qs = ds.variables['qs'][:]
        t = ds.variables['time'][:]# / 3600.
        delta_t = t[1] - t[0]
        qs_cumsum = np.cumsum(qs[:, 0, -2, 0])*delta_t
        qs_normalize490 = qs_cumsum/(rhog*(1-porosity))
    qs_cumsum = []
    
    fig, ax = plt.subplots(3, 1, figsize=(10,6))
    fig.subplots_adjust(left = 0.14, right=0.75, bottom = 0.1, hspace=0.5, wspace = 0.5)
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
    if 'year' in case:
        t = t/(3600*24) # convert time to days  
        # ax[2].set_ylabel('Difference in cumulative\nsediment flux (%)')
        ax[2].set_xlabel('Time (days)')
        ax[0].set_xlim([0,365])
        ax[1].set_xlim([0,365])
        ax[2].set_xlim([0,365])
        ax[1].set_ylim([-8,2])
        ax[2].set_ylim([-8,2])
    for i in range(len(qs_cumsum_all)):
        qs_cumsum_single = qs_cumsum_all[i,:]/(rhog*(1-porosity))
        if i == 9:
            ax[1].plot(t, (qs_cumsum_single-qs_normalize230)/np.max(qs_normalize230)*100, color=colors[i], linewidth=2, label = labels[i])
        elif i == 10:
            ax[1].plot(t, (qs_cumsum_single-qs_normalize240)/np.max(qs_normalize240)*100, color=colors[i], linewidth=2, label = labels[i])
        elif i == 11:
            ax[1].plot(t, (qs_cumsum_single-qs_normalize250)/np.max(qs_normalize250)*100, color=colors[i], linewidth=2, label = labels[i])
        elif i == 6:
            ax[2].plot(t, (qs_cumsum_single-qs_normalize470)/np.max(qs_normalize470)*100, color=colors[i], linewidth=2, label = labels[i])
        elif i == 7:
            ax[2].plot(t, (qs_cumsum_single-qs_normalize480)/np.max(qs_normalize480)*100, color=colors[i], linewidth=2, label = labels[i])
        elif i == 8:
            ax[2].plot(t, (qs_cumsum_single-qs_normalize490)/np.max(qs_normalize490)*100, color=colors[i], linewidth=2, label = labels[i])            
        elif i == 0: 
            ax[1].plot(t, (qs_cumsum_single-qs_normalize230)/np.max(qs_normalize230)*100, color=colors[i], linewidth=2, label = labels[i])
        elif i == 3: 
            ax[2].plot(t, (qs_cumsum_single-qs_normalize470)/np.max(qs_normalize470)*100, color=colors[i], linewidth=2, label = labels[i])
        ax[1].legend(bbox_to_anchor=(1.03, 1.09))
        ax[2].legend(bbox_to_anchor=(1.03, 1.09))
    
    fig.text(0.05, 0.33, 'Difference in cumul.\nsediment flux (%)', va='center', ha='center', rotation='vertical', fontsize=18)
        
    ax[0].plot(t, uw[:,0,0], linewidth=1.5)
    ax[0].hlines(calc_ut(0.000250), 0, 365, color = 'grey', label = '250 $\mu$m' )
    ax[0].hlines(calc_ut(0.000500), 0, 365, color = 'k', label = '500 $\mu$m')
    ax[0].set_ylabel('Wind\nspeed (m/s)')
    ax[0].legend(loc = (1.045, 0.18), title = 'Threshold velocity') 
    
    plt.gcf().text(0.15, 0.89, 'a)', fontsize=16, weight = 'bold')
    plt.gcf().text(0.15, 0.60, 'b)', fontsize=16, weight = 'bold')
    plt.gcf().text(0.15, 0.31, 'c)', fontsize=16, weight = 'bold')
    
    plt.show()
    plt.savefig(model_directory + '/../../analysis/figures/' + model_directory[73:] + '.png')     

#%% Go through all files in directory and visualize for each run
cases = get_cases(model_directory1)       

labels = ['reference D$_{50}$', '240 $\mu$m','250 $\mu$m','reference D$_{50}$','480 $\mu$m', '490 $\mu$m', 'PSD 4', 'PSD 5', 'PSD 6', 'PSD 1', 'PSD 2', 'PSD 3']
colors = ['k', 'k', 'k', 'k', 'k', 'k', 'darkgreen', 'green', 'lightgreen', 'darkblue', 'blue', 'lightblue']
            
fluxes_overview_1day_shapesnorm_detailed(model_directory1, cases, labels, colors)

#%% Go through all files in directory and visualize for each run
cases = get_cases(model_directory2)       

labels = ['reference D$_{50}$', '240 $\mu$m','250 $\mu$m','reference D$_{50}$','480 $\mu$m', '490 $\mu$m', 'PSD 4', 'PSD 5', 'PSD 6', 'PSD 1', 'PSD 2', 'PSD 3']
colors = ['k', 'k', 'k', 'k', 'k', 'k', 'darkgreen', 'green', 'lightgreen', 'darkblue', 'blue', 'lightblue']
            
fluxes_overview_1year_shapesnorm_detailed(model_directory2, cases, labels, colors)

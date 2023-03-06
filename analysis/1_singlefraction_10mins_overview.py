# -*- coding: utf-8 -*-
"""
Created on Fri Jul  8 15:54:42 2022

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

 #%%% Define transport equations

def calc_z0(ks):
    z0    = ks / 30.
    return z0

def calc_ustar(uw, z0):
    kappa                   = 0.41
    z                       = 10.
    ustar = (kappa * uw) / np.log(z/z0)
    return ustar

def calc_ustar_th(d):
    Aa                      = 0.085
    g                       = 9.81               					# [m/s^2] Gravitational constant
    rhoa                    = 1.225              					# [kg/m^3] Air density
    rhog                    = 2650.0                            # [kg/m^3] Grain density
    ustar_th = Aa * np.sqrt(((rhog - rhoa) / rhoa) * g * d)
    return ustar_th
    
def calc_Q(uw, d, z0):
    C                       = 1.5
    g                       = 9.81               					# [m/s^2] Gravitational constant
    rhoa                    = 1.225              					# [kg/m^3] Air density

    ustar = calc_ustar(uw, z0)
    ustar_th = calc_ustar_th(d)
    
    Q = C * (rhoa / g) * (ustar - ustar_th)**3 * 600 # * np.sqrt(d / D) 
    if isinstance(Q, np.ndarray) is False:
        if Q < 0:
            Q = 0
    else:
        Q[Q<0] = 0
    return Q

#%%% Visualize transport for range of single fraction cases, including subplot with shear velocities 
    
model_directory = r"C:\Users\cijzendoornvan\Documents\DuneForce\AEOLIS\GrainSize\scenarios\A_singlefraction_10mins"
cases = get_cases(model_directory) 
    
fig1, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))  
fig1.subplots_adjust(left = 0.15, right=0.95, bottom = 0.15, wspace = 0.3)  

flux_250 = []
flux_375 = []
flux_500 = []
rhog = 2650
porosity = 0.4    

for i, case in enumerate(cases):
    # Load parameters resulting from model run
    with netCDF4.Dataset(os.path.join(model_directory+'/', case), 'r') as ds:
        t = ds.variables['time'][:]# / 3600.
        qs = ds.variables['qs'][:]
        uw = ds.variables['uw'][:,:]
        delta_t = t[1] - t[0]
              
        # Calculate cumulative flux over time at the end of the domain
        qs_cumsum = np.cumsum(qs[:, 0, -2, 0])*delta_t/(rhog*(1-porosity))

        if '250' in case:
            flux_250.append(qs_cumsum[-1])
        elif '375' in case:
            flux_375.append(qs_cumsum[-1])
        elif '500' in case:
            flux_500.append(qs_cumsum[-1])
            
uw = list(range(5,31))
ax1.plot(uw, sorted(flux_250), 'o-', color = 'silver', markersize=4, label = '250 $\mu$m')               
ax1.plot(uw, sorted(flux_375), 'o-', color = 'grey', markersize=4, label = '375 $\mu$m')    
ax1.plot(uw, sorted(flux_500), 'o-', color = 'k', markersize=4, label = '500 $\mu$m')   
# ax['bottom'].set_title('Transport at 30 mins with 375 $\mu$m and varying wind speeds', fontsize = 20) #, fontweight = 'bold'
ax1.legend()
ax1.set_xlabel('Wind speed (m/s)')
ax1.set_ylabel('Cumulative sediment\nflux ($m^3$/m)')
        
d_ini = [50, 75, 110, 150, 195,  235,  285,  345,  420,  505,  615,  745, 900, 1100, 1325, 1600, 2000]
d = [dx*1e-6 for dx in d_ini]
u = [10, 12.5 , 15]

colors2 = ['silver', 'darkgrey', 'grey', 'black']
labels2 = ['u$_w$ = 10 m/s', 'u$_w$ = 12.5 m/s', 'u$_w$ = 15 m/s']
# Calculate (u star - u star th) - gs dependent, varying d
ustar_ar = np.zeros((3,17))
ustar_th_ar = []
for j, gs in enumerate(d):
    for i, uw in enumerate(u):
        z0 = calc_z0(gs)    
        ustar = calc_ustar(uw, z0)
        ustar_ar[i,j] = ustar
for j, gs in enumerate(d):    
    ustar_th = calc_ustar_th(gs)
    ustar_th_ar.append(ustar_th)
    
ax2.plot(d_ini, sorted(ustar_th_ar), '-o', color = 'red', markersize=4, label = 'u$_{*, th}$')  
ax2.plot([0.2],  marker='None', linestyle='None', label='u$_*$ based on:')             
ax2.plot(d_ini, sorted(ustar_ar[0,:]), '-o', color = 'silver', markersize=4, label = 'u$_w$ = 10 m/s')    
ax2.plot(d_ini, sorted(ustar_ar[1,:]), '-o', color = 'grey', markersize=4, label = 'u$_w$ = 12.5 m/s')           
ax2.plot(d_ini, sorted(ustar_ar[2,:]), '-o', color = 'k', markersize=4, label = 'u$_w$ = 15 m/s')           
    
ax2.set_xlabel('Grain size ($\mu$m)')
ax2.set_ylabel('Shear velocity (m/s)')
ax2.legend()

h,l = ax2.get_legend_handles_labels()
by_label = dict(zip(l, h))
ax2.legend(by_label.values(), by_label.keys())

plt.gcf().text(0.16, 0.9, 'a)', fontsize=16, weight = 'bold')
plt.gcf().text(0.61, 0.9, 'b)', fontsize=16, weight = 'bold')
                    
plt.savefig(model_directory + '/../../analysis/figures/' + 'singlefraction_10mins.png')



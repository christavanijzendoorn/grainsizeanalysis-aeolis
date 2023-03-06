# -*- coding: utf-8 -*-
"""
Created on Wed Sep 21 11:21:49 2022

@author: cijzendoornvan
"""

import numpy as np

def calc_median_grain_size(fractions, mass, percent):
    '''Calculate grain size characteristics based on mass in each fraction
    Calculate grain size distribution for each cell based on weight 
    distribution over the fractions. Interpolates to the requested percentage 
    in the grain size distribution. For example, percent=50 will result 
    in calculation of the D50. Calculation is only executed for the top layer.
    '''
    from scipy.interpolate import interp1d
    # mass = mass[:,:,:,0,:] # only select upper, surface layer because that is the relevant layer for transport
    D = np.zeros((mass.shape[0], mass.shape[1], mass.shape[2], mass.shape[3]))
    for ti in range(mass.shape[0]):
        for yi in range(mass.shape[1]):
            for xi in range(mass.shape[2]):
                for li in range(mass.shape[3]):
                    diameters = np.insert(fractions, 0, 0)
                    cummulative_weights = np.cumsum(np.insert(mass[ti, yi, xi, li, :], 0, 0))
                    percentages = 100*cummulative_weights/np.max(cummulative_weights)
                    f_linear = interp1d(list(percentages), diameters, fill_value='extrapolate') # get interpolation function
                    
                    # Retrieve grain size characteristics based on interpolation
                    D[ti, yi, xi, li] = f_linear(percent)*1e6
    return D

def calc_mean_grain_size(fractions, mass):
    '''Calculate mean grain size based on mass in each fraction
    Calculate mean grain size for each cell based on weight distribution 
    over the fractions. Calculation is only executed for the top layer.
    '''
    D_mean = np.zeros((mass.shape[0], mass.shape[1], mass.shape[2], mass.shape[3]))
    for ti in range(mass.shape[0]):
        for yi in range(mass.shape[1]):
            for xi in range(mass.shape[2]):
                for li in range(mass.shape[3]):
                    diameters = fractions
                    weights = mass[ti, yi, xi, li, :]/ np.sum(mass[ti, yi, xi, li, :])
                    
                    # Retrieve mean grain size based on weight of mass
                    D_mean[ti, yi, xi, li] = np.sum(diameters*weights)*1e6
    return D_mean

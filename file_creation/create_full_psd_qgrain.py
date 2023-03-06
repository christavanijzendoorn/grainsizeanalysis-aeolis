# -*- coding: utf-8 -*-
"""
Created on Mon Oct 10 11:09:29 2022

@author: cijzendoornvan
"""

###################
# This script is run in a separate environment so no conflicts occur with the AeoLiS package
# Install Qgrain, numpy, scipy and matplotlib in this environment
###################

from QGrain.statistics import all_statistics
from QGrain.generate import random_sample, random_dataset, SIMPLE_PRESET
from QGrain.models import DistributionType

import numpy as np
import pandas as pd
import scipy
from scipy.stats import skewnorm, genhyperbolic, rayleigh
from scipy.interpolate import interp1d

import matplotlib.pyplot as plt
import math

savedir = "C:/Users/cijzendoornvan/Documents/DuneForce/AEOLIS/GrainSize/data/"

x = [0.000050, 0.000150, 0.000250, 0.000350, 0.000450, 0.000550, 0.000650, 0.000750, 0.000850, 0.000950, 0.001050, 0.001150, 0.001250, 0.001350, 0.001450, 0.001550, 0.001650, 0.001750, 0.001850, 0.001950]
diameters = x

diam_mm = [1000*diam for diam in x]
phi = [-math.log2(diam) for diam in diam_mm]

perc_cummul = [0.02, 0.05, 0.10, 0.16, 0.25, 0.36, 0.50, 0.64, 0.75, 0.84, 0.90, 0.95, 0.98]
perc_cummulative = [x*100 for x in perc_cummul]

#%% Get statistics from Noordwijk grain size data

def get_ratios(grainsize, percentage_cumm, diameters, loc):
    # add known values of 0 percent at 0 mm and 100% at 2000 micrometer.
    gs_zero = np.insert(grainsize,[0,len(grainsize)],[0,2000])
    percentage_cumm = np.insert(percentage_cumm,[0,len(percentage_cumm)],[0,1])
    
    f_linear = interp1d(gs_zero, percentage_cumm)
    
    ratios = np.zeros(len(diameters))
    for i, f in enumerate(diameters):
        if i == 0:
            ratios[i] = f_linear(f*1000000)
        else: 
            ratios[i] = f_linear(f*1000000) - f_linear(diameters[i-1]*1000000)

    ratios = ratios / sum(ratios) # Normalize so sum equals unity
    
    with open("grainsize_" + loc + ".txt", "w") as output:
        for item in ratios:
            output.write(str(item) + " ")
            
    return ratios

#%% Get statistics from Noordwijk grain size data
  
# Grain size for each DXX as based on the average of all samples per location.
grainsizes_NW = [125, 162, 188, 210, 237, 265, 298, 334, 369, 407, 470, 606, 1112]
ratios_NW = get_ratios(grainsizes_NW, perc_cummul, diameters, 'NW')
plt.plot(x, ratios_NW, 'grey')

statistics = all_statistics(np.array(diam_mm), np.array(phi), ratios_NW)

used_statistics =  statistics['logarithmic']
skewness = round(used_statistics['skewness'],2)
mean = round(used_statistics['mean'], 1)
std = round(used_statistics['std'], 2)
kurtosis = round(used_statistics['kurtosis'], 2)

#%% Fill parameters needed to create distributions with Qgrain

# a = skewness parameters
# loc = mean grain size
# scale = stddev

# Mean and stand deviation are varied to create 9 different shapes. Note that only 6 of these are used in the paper.

preset1 = dict(target=[
    [(skewness, 0.0), (mean, 0.0), (std-0.3, 0.0), (kurtosis, 0.0)]],
    distribution_type=DistributionType.SkewNormal)

preset2 = dict(target=[
    [(skewness, 0.0), (mean, 0.0), (std, 0.0), (kurtosis, 0.0)]],
    distribution_type=DistributionType.SkewNormal)

preset3 = dict(target=[
    [(skewness, 0.0), (mean, 0.0), (std+0.3, 0.0), (kurtosis, 0.0)]],
    distribution_type=DistributionType.SkewNormal)

preset4 = dict(target=[
    [(skewness, 0.0), (mean-0.5, 0.0), (std-0.3, 0.0), (kurtosis, 0.0)]],
    distribution_type=DistributionType.SkewNormal)

preset5 = dict(target=[
    [(skewness, 0.0), (mean-0.5, 0.0), (std, 0.0), (kurtosis, 0.0)]],
    distribution_type=DistributionType.SkewNormal)

preset6 = dict(target=[
    [(skewness, 0.0), (mean-0.5, 0.0), (std+0.3, 0.0), (kurtosis, 0.0)]],
    distribution_type=DistributionType.SkewNormal)

preset7 = dict(target=[
    [(skewness, 0.0), (mean+0.5, 0.0), (std-0.3, 0.0), (kurtosis, 0.0)]],
    distribution_type=DistributionType.SkewNormal)

preset8 = dict(target=[
    [(skewness, 0.0), (mean+0.5, 0.0), (std, 0.0), (kurtosis, 0.0)]],
    distribution_type=DistributionType.SkewNormal)

preset9 = dict(target=[
    [(skewness, 0.0), (mean+0.5, 0.0), (std+0.3, 0.0), (kurtosis, 0.0)]],
    distribution_type=DistributionType.SkewNormal)


#%% Create grain size distributions with Qgrain

def create_custom_graindist(preset, savedir, filename):
    dataset = random_dataset(**preset, n_samples=1,
                         min_size=50, max_size=1950.0, n_classes=20,
                         precision=2, noise=2)
    with open(savedir + filename + ".txt", "w") as output:
        for i, item in enumerate(dataset.components[0,0,:]):
            if len(dataset.components[0,0,:]) == i+1:
                output.write(str(item))
            else:
                output.write(str(item) + " ")
            
    percentages = np.cumsum(dataset.components[0,0,:]/np.sum(dataset.components[0,0,:]))*100
    f_linear = interp1d(list(percentages), dataset.classes, fill_value='extrapolate') # get interpolation function
    
    print(dataset.classes*1e-6)
    # Retrieve grain size characteristics based on interpolation
    d10 = f_linear(10)
    d50 = f_linear(50)
    d90 = f_linear(90)                
    return d10, d50, d90

d10_1, d50_1, d90_1 = create_custom_graindist(preset1, savedir, 'shape1')
d10_2, d50_2, d90_2 = create_custom_graindist(preset2, savedir, 'shape2')
d10_3, d50_3, d90_3 = create_custom_graindist(preset3, savedir, 'shape3')
d10_4, d50_4, d90_4 = create_custom_graindist(preset4, savedir, 'shape4')
d10_5, d50_5, d90_5 = create_custom_graindist(preset5, savedir, 'shape5')
d10_6, d50_6, d90_6 = create_custom_graindist(preset6, savedir, 'shape6')
d10_7, d50_7, d90_7 = create_custom_graindist(preset7, savedir, 'shape7')
d10_8, d50_8, d90_8 = create_custom_graindist(preset8, savedir, 'shape8')
d10_9, d50_9, d90_9 = create_custom_graindist(preset9, savedir, 'shape9')

#%% Visualize particle size distributions that are used in the paper 
# Note that shapes 1-3, which correspond to the mean grain size of the original Noordwijk distribution are not used.
# The shapes 4-6 correspond to psd 4-6 but shapes 7-9 correspond psd 1-3 in the figure.

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

model_directory = r"C:\Users\cijzendoornvan\Documents\DuneForce\AEOLIS\GrainSize\data"

matplotlib.rcParams['xtick.labelsize'] = 16
matplotlib.rcParams['ytick.labelsize'] = 16

colors = ['darkgreen', 'green', 'lightgreen', 'darkblue', 'blue', 'lightblue']
labels = ['PSD 4', 'PSD 5', 'PSD 6', 'PSD 1', 'PSD 2', 'PSD 3'] 
          
diameters = [50., 60.63316614, 73.52761673, 89.16424402,  108.12620843,
  131.12068721,  159.00524824,  192.81983268,  233.82553901,  283.55165509,
  343.85269226,  416.97754837,  505.65337935,  613.18730721,  743.58975748,
  901.72402614, 1093.48765382, 1326.03237177, 1608.03082214, 1950.]

txtfiles = os.listdir(model_directory) 

fig = plt.figure(figsize=(6,5))
fig.subplots_adjust(left = 0.17, right=0.95, bottom = 0.15)

i=0
j=0
for txt in txtfiles:
    if 'shape' in txt:
        print(txt)
        if i>2: # Skip shapes with original mean grain size.
            dist = np.loadtxt(model_directory + '\\' + txt, delimiter = ' ')
            plt.plot(diameters, dist, color = colors[j], label = labels[j], linewidth = 2)
            j+=1
        i+=1
plt.legend()

handles, labels = plt.gca().get_legend_handles_labels()
order = [3,4,5,0,1,2]
plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order], fontsize = 16)

plt.xlabel('Particle size ($\mu$m)', fontsize = 20)
plt.ylabel('Weight percentage (%)', fontsize = 20)


# -*- coding: utf-8 -*-
"""
Created on Thu Sep 22 16:57:01 2022

@author: cijzendoornvan
"""

# Load packages
import numpy as np
import matplotlib.pyplot as plt

import os
import sys
sys.path.append(os.getcwd())

# Initialize directories
data_directory = r"C:\Users\cijzendoornvan\Documents\DuneForce\AEOLIS\reference_cases\data"
aeolis_directory = r"C:\Users\cijzendoornvan\OneDrive - Delft University of Technology\Documents\DuneForce\AEOLIS\aeolis-python\aeolis"

sys.path.append(aeolis_directory)
from console_debug import *
from aeolis.model import WindGenerator

###############
#The wind files were created with the following commands. 
#Note that it is a random generator, so it will not give the same wind series as were used in our research.
###############

#mean_speed=9.0, max_speed=30.0, dt=60.0, n_states=30, shape=2.0, scale=2.0
# wind = WindGenerator(mean_speed = 10.).generate(duration=3600*24)
# wind.write_time_series(data_directory + '/' + 'wind_1day_.txt')
# wind.plot()

# wind = WindGenerator(mean_speed = 10., dt = 3600.).generate(duration=3600*24*365)
# #mean_speed=9.0, max_speed=30.0, dt=60.0, n_states=30, shape=2.0, scale=2.0
# wind.write_time_series(data_directory + '/' + 'wind_1year_.txt')
# wind.plot()


#%% Visualize

t_10min = range(0,600,1)
t_day = range(0,86400,60)
t_day = range(0,31536000,3600)

wind_10min_5ms = np.loadtxt(data_directory + '//' + 'wind_constant_5ms.txt')
wind_10min_10ms = np.loadtxt(data_directory + '//' + 'wind_constant_10ms.txt')
wind_10min_20ms = np.loadtxt(data_directory + '//' + 'wind_constant_20ms.txt')
wind_day = np.loadtxt(data_directory + '//' + 'wind_1day.txt')
wind_year = np.loadtxt(data_directory + '//' + 'wind_1year.txt')

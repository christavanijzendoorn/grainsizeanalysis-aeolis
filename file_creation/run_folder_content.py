# -*- coding: utf-8 -*-
"""
Created on Wed Sep 21 11:29:10 2022

@author: cijzendoornvan
"""

# Load packages
import os
import sys
sys.path.append(os.getcwd())

#%% Run all aeolis cases that are stored in a folder

# Initialize directories
model_directory = r"C:\Users\cijzendoornvan\Documents\DuneForce\AEOLIS\GrainSize\scenarios\A_singlefraction_10mins"
# model_directory = r"C:\Users\cijzendoornvan\Documents\DuneForce\AEOLIS\GrainSize\scenarios\B_singlefraction_5050mix_1day"
# model_directory = r"C:\Users\cijzendoornvan\Documents\DuneForce\AEOLIS\GrainSize\scenarios\C_8020mixes_1year"
# model_directory = r"C:\Users\cijzendoornvan\Documents\DuneForce\AEOLIS\GrainSize\scenarios\D_variedpercentagemixes_1year"
# model_directory = r"C:\Users\cijzendoornvan\Documents\DuneForce\AEOLIS\GrainSize\scenarios\E_fullPSDs_1day"
# model_directory = r"C:\Users\cijzendoornvan\Documents\DuneForce\AEOLIS\GrainSize\scenarios\E_fullPSDs_1year"
# model_directory = r"C:\Users\cijzendoornvan\Documents\DuneForce\AEOLIS\GrainSize\scenarios\F_crossshore_10mins"
# model_directory = r"C:\Users\cijzendoornvan\Documents\DuneForce\AEOLIS\GrainSize\scenarios\F_crossshore_1day"
# model_directory = r"C:\Users\cijzendoornvan\Documents\DuneForce\AEOLIS\GrainSize\scenarios\F_crossshore_1year"
# model_directory = r"C:\Users\cijzendoornvan\Documents\DuneForce\AEOLIS\GrainSize\scenarios\M_verticallayering_10mins"
# model_directory = r"C:\Users\cijzendoornvan\Documents\DuneForce\AEOLIS\GrainSize\scenarios\M_verticallayering_1day"
# model_directory = r"C:\Users\cijzendoornvan\Documents\DuneForce\AEOLIS\GrainSize\scenarios\M_verticallayering_1year"
aeolis_directory = r"C:\Users\cijzendoornvan\Documents\DuneForce\AEOLIS\aeolis-python\aeolis"

sys.path.append(aeolis_directory)
from console_debug import *

# Go through all files in directory and run if it is not a directory
for txtfile in os.listdir(model_directory):
    if os.path.isdir(model_directory + "/" + txtfile):
        pass
    elif os.path.getsize(model_directory + "/" + txtfile) > 9000:
        pass
    else: 
        case_dir = model_directory + "/" + txtfile
        aeolis_debug(case_dir)
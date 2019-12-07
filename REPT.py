# -*- coding: utf-8 -*-
"""
Created on Fri Nov 29 18:15:54 2019

@author: Aggelos
"""

import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
from spacepy import pycdf
from scipy.interpolate import interp1d
from scipy.interpolate import Rbf
startTime = datetime.now()                      #This exists for the purpose of calculating the time needed to run the script

Timestamps = 7731
TotalPitchAngles = 17
TotalChannels = 12

#Reading L3 Data File for REPT
data_rept = pycdf.CDF('rbspa_rel03_ect-rept-sci-L3_20170422_v5.4.0.cdf')

#Reading L2 Data File for REPT Quality Flags)
rept_flags = pycdf.CDF('rbspa_rel03_ect-rept-sci-L2_20170422_v5.4.0.cdf')

#Importing Data from L3 Data File
Epoch = data_rept['Epoch'][:]
pitchangle = data_rept['FEDU_Alpha'][:]
channel_energies = data_rept['FEDU_Energy'][:]
rept_flux = data_rept['FEDU'][:]


#Importing Magnetic Coordinates from L3 Data File:
#######################################
Lstar = data_rept['L_star'][:]
L = data_rept['L'][:]
MLT = data_rept['MLT'][:]
B_Calc = data_rept['B_Calc'][:]
B_eq = data_rept['B_Eq'][:]
#######################################

#Converting Epoch from Datetime to time
Time_rept = np.zeros(7731, dtype=float)                 #Initializing an array 7731x1 full of zeros
          

for i in range(len(data_rept['Epoch'])):
   Time_rept[i] = datetime.timestamp(Epoch[i])          #Every array element takes the value of its correspoding time
 
    
#Importing Quality Flags from L2 Data File:
rept_quality_flags = rept_flags['Data_Quality'][:]

rept_quality_flags = rept_quality_flags.astype(float)

##############################################################################
rept_flux_corr = np.zeros((7734, 17, 12))
                                                #Creating an identical copy of the fedu_corr matrix
rept_flux_corr = rept_flux.copy()

Error_Min = 0                                   #Desired min-max error value
Error_Max = 50
Quality_Flag_Value = 1                          #Desired quality flag value

for i in range(7731):
    for j in range(TotalPitchAngles):
        for k in range(TotalChannels):
            if np.abs(rept_quality_flags[i]) >= Quality_Flag_Value:                             #Clearing fedu_corr_clear matrix from values with lower than desired quality flag
                rept_flux_corr[i][j][k] = np.nan      
            if rept_flux_corr[i][j][k] <= 0:                                                    #Making sure every element of the 7731x17x12 matrix is either positive or NaN
                rept_flux_corr[i][j][k] = np.nan










rept_angles_copy = pitchangle.copy()
interpolate1 = np.zeros((17,1), dtype=float)
negative = 0
interpolatedpitchangles = np.zeros((7731,11,12), dtype=float)

for i in range(7731):
    rept_angles_copy = pitchangle.copy()
    interpolate1 = np.zeros((17,1), dtype=float)
    negative = 0
    
    for k in range(12):
        rept_angles_copy = pitchangle.copy()
        interpolate1 = np.zeros((17,1), dtype=float)
        negative = 0
        
        for x in range(17):
             interpolate1[x] = rept_flux_corr[i][x][k]
        for y in range(17):
            if np.isnan(interpolate1[y]): rept_angles_copy[y] = np.nan
        for z in range(17):
            if np.isnan(interpolate1[z]): negative = negative + 1
                
                      
        if any([np.isnan(interpolate1[0]), np.isnan(interpolate1[16]), negative >= 9]):
            for w in range(11): interpolatedpitchangles[i][w][k] = np.nan           
               
        else:
            interpolate1 = interpolate1[np.logical_not(np.isnan(interpolate1))]
            rept_angles_copy = rept_angles_copy[np.logical_not(np.isnan(rept_angles_copy))]
            interpolate1 = np.asarray(interpolate1).squeeze()
            rept_angles_copy = np.asarray(rept_angles_copy).squeeze()
            func = interp1d(rept_angles_copy, interpolate1, kind='quadratic', fill_value='extrapolate') 
           
            for j in range(11):
                interpolatedpitchangles[i][j][k] = func(pitchangle_mageis[j])
            
print(datetime.now() - startTime)                                                               #Printing the time needed for the script to run







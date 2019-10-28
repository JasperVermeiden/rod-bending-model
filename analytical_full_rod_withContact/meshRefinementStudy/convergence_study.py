# -*- coding: utf-8 -*-
"""
Created on Mon Sep 30 16:16:38 2019

@author: Jasper Vermeiden
Convergence study based on:
    P_opt       = 28	                # [W]        Power of optical module, light energy input [W]
    Eff_opt     = 0.33                  # [-]        Optical efficiency, percentage of light input converted to output
    width       = 1.6*mm	            # [m]        Rod width
    height      = 1.05*mm               # [m]        Rod height (thickness)
    length      = 59*mm                 # [m]        Rod length
    k_rod 	    = 4.5	                # [W/mK]     Internal heat conductivity of rod
    k_air	    = 0.0265                # [W/mK]     Conductivity of air
    T_RH        = 70	                # [°C]       Rod holder temp, assumed uniform
    
    avgRayAngle = 23                    # [deg]      Average angle rays entering the rod 
    refl_eff    = .9                    # [-]        Reflection efficiency, assumed to be constant for all wavelengths
    lightsource = 'padbar'              # [-]        Type of LED lightsource used, can be 'CSP' or 'padbar'
    frac_body   = .11                   # [-]        Fraction of heat load assigned as uniform body heat load
"""

import numpy as np
import matplotlib.pyplot as plt

meshes  = np.arange(10,110,10)
for size in meshes:
    temperatures = np.load('temperatures_'+str(size)+'.npy')
    plt.figure(1, figsize = (12,6))
    plt.plot(np.arange(0,1.6, 1.6/size),temperatures[:,1], label = 'N elements  ' + str(size))

plt.legend()
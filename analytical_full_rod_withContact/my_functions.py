# -*- coding: utf-8 -*-
"""
Created on Tue Sep  3 11:42:49 2019

@author: 70238599
"""

import numpy as np
import matplotlib.pyplot as plt
from my_datastructures import *

def discretise(rod):
    # create discretised version of the rod, first element located on left edge
    dx = rod.length/rod.num_el_x
    elements = []
    for el_numb in range (rod.num_el_x+1): 
        x_pos =  el_numb*dx
        elements.append(Element(el_numb, x_pos ))
    rod.set_elements(elements)
    

def thermalSolution(rod, Q_heat, T_RH, k_air, roughnessGap):
    
    # Temperature difference over the height of the rod, required to transport all incoming heat
    dT_rod =    Q_heat* rod.height/(rod.k* rod.length* rod.width)
    x_positions, conj_moment, heatDeflection, totalDisplacement = rod.displacements()
    for i, element in enumerate (rod.elements) :
        dT_gap_i = Q_heat*(totalDisplacement[i] + roughnessGap)/(k_air* rod.length* rod.width)
        element.dT_gap = dT_gap_i
        element.T_bottom =  T_RH + dT_gap_i
        element.T_top    = element.T_bottom + dT_rod
    plotTempDistribution(rod, x_positions)    

def plotDisplacements(rod, force_posF_posC, x_spring):
    x_positions, conj_moment, heatDeflection, totalDisplacement = rod.displacements()
    maxDeflection, loc_maxDeflection, minDeflection, loc_minDeflection  = rod.max_min_deflection()
    nodeIDreaction   = int(force_posF_posC[1])
    x_reaction       = x_positions[nodeIDreaction]
   
    plt.figure(2, figsize = (12,4))
    plt.grid('true')
    plt.plot(x_positions,conj_moment)
    plt.plot(x_positions,heatDeflection)
    plt.plot(x_positions,totalDisplacement)
    plt.plot(loc_maxDeflection,maxDeflection,"v")
    plt.plot(loc_minDeflection,minDeflection,"v")
    plt.axvline(x= rod.length/2, linestyle = '--', color = 'black', linewidth = 1)   
    plt.title('Displacement [m]')
    plt.legend(['Spring deflection', 'Heat deflection', 'total deflection', 'min gap', 'max gap', 'Sym axis'])
        
    plt.figure(5, figsize = (12,3))
    plt.grid('true')
    plt.plot(x_positions*1000,totalDisplacement*1E6, label = 'Resultant displacement')
    plt.plot(loc_maxDeflection*1000,maxDeflection*1E6,"v",label = 'maximum gap')
    plt.plot(loc_minDeflection*1000,minDeflection*1E6,"v",label = 'minimum gap')
    plt.plot(x_reaction*1000,totalDisplacement[nodeIDreaction]*1E6,"v",label = 'application point reaction force')
    plt.axvline(x = x_spring*1000, linestyle = '--', color = 'black', linewidth = 1, label = 'spring location')   

    plt.title('Displacement [$\mu$m]')
    plt.xlabel('Rod position [mm]')
    plt.ylabel('Displacement [$\mu$m]')
    
    plt.legend()

    print('maximum gap: %.4f [$\mu$m]' %(maxDeflection*1e6))
    print('minimum gap: %.4f [$\mu$m]' %(minDeflection*1e6))

def plotTempDistribution(rod, x_coordinates, T_2D):
    y_coordinates = np.arange(rod.height-rod.height/rod.num_el_y,-rod.height/rod.num_el_y,-rod.height/rod.num_el_y)
    ## Find locations of min and max temperature
    min_y, min_x = np.where(T_2D <= np.min(T_2D)+ 0.000001)
    max_y, max_x = np.where(T_2D >= np.max(T_2D)- 0.000001)
    min_x, min_y = x_coordinates[min_x], y_coordinates[min_y]
    max_x, max_y = x_coordinates[max_x], y_coordinates[max_y]
    
    plt.figure(4,figsize = (20,6))
    plt.imshow(T_2D, extent = [0,np.max(x_coordinates)*1000,0,np.max(y_coordinates)*1000])
    plt.colorbar(orientation = 'vertical')
    
    '''
    for i, y in enumerate (min_y):
        plt.plot(min_x[i]*1000,y*1000,'^', color = 'y', label = 'T max = ' + str(np.min(T_2D)))
    for i, y in enumerate (max_y):
        plt.plot(max_x[i]*1000,y*1000,'^', color = 'r', label = 'T max = ' + str(np.max(T_2D)))
    '''
  
    plt.title('Temperature distribution')
    plt.ylabel('[mm]')
    plt.xlabel('[mm]')
    print('Max temperature = %.2f' % (np.max(T_2D)), '[C]')
    print('Min temperature = %.2f' % (np.min(T_2D)), '[C]')

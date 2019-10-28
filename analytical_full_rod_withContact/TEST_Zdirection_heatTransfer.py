# -*- coding: utf-8 -*-
"""
Created on Mon Sep 30 12:55:50 2019

@author: Jasper Vermeiden 
Rod temperature computation using thermal lumped mass model (discretised rod with each element a lumped mass) 
2D discretisation of rod geometry over width of rod. Bottom contacts RH, side walls proximity variable left and right. 
state space based solving routine 
"""


from my_datastructures   import Rod
from my_functions        import *
from my_absorbtion_functions import absorbedPower
import matplotlib.pyplot as plt
import numpy as np

mm = 1/1000

# Parameters, dimensions
P_opt       = 28	                # [W]        Power of optical module, light energy input [W]
Eff_opt     = 0.33                  # [-]        Optical efficiency, percentage of light input converted to output
width       = 1.05*mm	            # [m]        Rod width
height      = 1.05*mm               # [m]        Rod height (thickness)
length      = 59*mm                 # [m]        Rod length
k_rod 	    = 4.5	                # [W/mK]     Internal heat conductivity of rod
k_air	    = 0.0265                # [W/mK]     Conductivity of air
CTE_rod     = 8.00E-06	            # [1/K]      Expansion coefficient
T_RH        = 70	                # [°C]       Rod holder temp, assumed uniform
gap_unconstr= 5e-6                  # [m]        Assumed air gap thickness, unconstrained, undeformed rod
F_spring    = -6                    # [N]        Spring force 
#X_spring    = 14*mm                 # {m}        Distance of springs from edges
#Q_heat      = P_opt*(1-Eff_opt)     # [w]        Heat load, due to conversion lossses etc. 
Q_heat      = 22.5
E           = 3E+11                 # [Pa]       Modulus
#roughnessGap= 3e-6                  # [m]        Roughness gap between rod and rod holder

avgRayAngle = 23                    # [deg]      Average angle rays entering the rod 
refl_eff    = .9                    # [-]        Reflection efficiency, assumed to be constant for all wavelengths
lightsource = 'padbar'              # [-]        Type of LED lightsource used, can be 'CSP' or 'padbar'
frac_body   = .11                   # [-]        Fraction of heat load assigned as uniform body heat load
#R_corr      = 1.24                  # [-]        Radius correction factor  1 corresponds to toponly heat lod

## Analysis settings
Horizontal_transfer = True          #            Perform thermal analysis with or without lateral heat transfer
num_el_x    = 200                   #            Discretisation in x direction
num_el_y    = 30                    #            Discretisation in y direction
num_el_z    = 30                    #            Discretisation in z direction                                                               
heatload    = 'distribution'        #            Can be 'distribution' or 'toponly', with toponly the fraction body load is also applied to the top 


Q_body      = Q_heat*frac_body      # [W]        Heat load applied as body load
Q_distr     = Q_heat*(1-frac_body)  # [W]        Heat load applied distributed according to absorbtion spectrum

assert heatload == 'distribution' or heatload == 'toponly'

if heatload == 'toponly':    # only apply curvature correction factor if heatload is distributed over rod
    R_corr = 1

# Create data structure containing a finite element representation of the rod
rod = Rod(length, height, width, k_rod, CTE_rod,E, num_el_x, num_el_y, num_el_z)

discretise(rod)

# Optical calculations: stokes losses as function of position given absorbtion and LED emittance spectra
absorbedPower = absorbedPower(rod,lightsource, avgRayAngle, refl_eff, plot = True)

'''
Thermal solution over width of rod 
'''

# Extra variables: 
gap_bottom  = 3E-6
gap_left    = 5E-6
gap_right   = 5E-6
#Discretisation 
num_el_z = rod.num_el_z             # number of elements over the width of the rod
num_el_y = rod.num_el_y             # nuber of elements over the thickess (height) of the rod

   
## Thermal model
element_height= rod.height/num_el_y
element_width = rod.width/num_el_z

g_y     = element_width*rod.length*rod.k/element_height             # [W/K]      Conductivity of rod material per element thickness, for area of element 
g_z     = element_height*rod.length*rod.k/element_width             # [W/K]      Conductivity of rod material per element thickness, for area of element 
g_bottom= element_width*rod.length*k_air/gap_bottom                 # [W/k]      Conductivity of air gap with known thickness
g_left  = element_height*rod.length*k_air/gap_left                  # [W/k]      Conductivity of air gap with known thickness
g_right = element_height*rod.length*k_air/gap_right                 # [W/k]      Conductivity of air gap with known thickness

## Heat capacities matrix (E) 
E = np.zeros((num_el_y*num_el_z,num_el_y*num_el_z), float)
np.fill_diagonal(E, 1)      # set mc_p to 1 for all elements, value is irrelevant for steady state result, validated 26-09-2019
                            # Could also use np.identity(num_elements, float)

## Heat transfer matrix (K)
K = np.zeros((num_el_y*num_el_z,num_el_y*num_el_z), float)

# Add vertical heat transfer to K matrix
for j in range (num_el_z):
    # heat transfer to rod holder
    K[(j+1)*num_el_y-1, (j+1)*num_el_y-1] -= g_bottom
    
    # loop over all elments in vertical direction to add vertical heat transfer for each element    
    for i in range (num_el_y):
        index = j*num_el_y + i
        
        if (i-1) >= 0:
            K[index, index-1] += g_y
            K[index, index] -= g_y
            
        if (i +1) <= (num_el_y-1):    
            K[index, index+1] += g_y
            K[index, index] -= g_y

        # Add heat transfer to sides of rod holder for elements on the vertical edges
        if j == 0:
            K[index, index] -= g_left
        if j == num_el_z-1:
            K[index, index] -= g_right
        
if Horizontal_transfer:
    # Add Z direction heat transfer to K matrix
    for j in range (num_el_z):
        # loop over all elments in vertical direction to add horizontal heat transfer for each element    
        for i in range (num_el_y):
            index = j*num_el_y + i
            
            if (index-num_el_y) >= 0:
                K[index, index-num_el_y] += g_z
                K[index, index] -= g_z
                
            if (index + num_el_y) <= (num_el_y*num_el_z-1):    
                K[index, index+num_el_y] += g_z
                K[index, index] -= g_z

## Heat load input vector
if heatload == 'toponly': 
    Q_distr += Q_body
if heatload == 'distribution':
    Q_distr = Q_distr/sum(absorbedPower)    # normalise distributed Q so that the actual distributed heat load corresponds to this value

u = np.array(([T_RH], [Q_distr/num_el_z],[Q_body/(num_el_y*num_el_z)]))    
print(u)

## Heat load matrix (L)
L = np.zeros((num_el_z*num_el_y,np.size(u)), float)     

# Heat flux to rod holder bottom + sides 
for j in range (num_el_z):
    index = (j+1)*num_el_y-1 
    L[index,0] += g_bottom

L[[np.arange(0,num_el_y)],0] += g_left
L[[np.arange((num_el_z-1)*num_el_y, num_el_y*num_el_z)],0] += g_right

# Optical heat loads
if heatload == 'toponly':           # Uniform heat influx on rod top surface
    for j in range (num_el_z):
        index = j*num_el_y 
        L[index,1] = 1
    
if heatload == 'distribution':
    for j in range (num_el_z):
        index = list(np.arange(0,num_el_y)+j*num_el_y) 
        L[index,1] = [absorbedPower] 
    # Uniform heat load (reabsorbtion etc.)
    L[:,2] = 1




## state space solution, steady state
A = -np.linalg.inv(E)@K
B = -np.linalg.inv(E)@L
C = np.identity(num_el_z*num_el_y, float) # Interested in T values for all elements
D = np.zeros(np.shape(L),float)
x = -np.linalg.inv(A)@B@u
T_static = C@x +D@u

## Unpack 1D temperatures vector into a matrix representing the rod geometry
T_staticMatrix = np.zeros((num_el_y, num_el_z))
for i in range (num_el_z):
    T_staticMatrix[:,i] = list(T_static[i*num_el_y:(i+1)*num_el_y ]) 

q_in  = sum(L[:,1:]@u[1:])


plotTempDistribution(rod,np.arange(element_width,rod.width+ element_width, element_width),T_staticMatrix)



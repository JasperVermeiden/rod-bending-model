"""
Analytical rod bending, original by Pedro Favera
Implementation in python by Jasper Vermeiden, 
started 03-09-2019
"""

from my_datastructures   import Rod
from my_functions        import *
from my_structural_functions      import *
from my_2D_lumpedMassFEmodel import thermalSolution2D
from my_absorbtion_functions import absorbedPower
import matplotlib.pyplot as plt
import numpy as np
import time
        
mm = 1/1000

# Parameters, dimensions
P_opt       = 21.06	                # [W]        Power of optical module, light energy input [W]
Eff_opt     = 0.33                  # [-]        Optical efficiency, percentage of light input converted to output
width       = 1.05*mm	            # [m]        Rod width
height      = 1.05*mm               # [m]        Rod height (thickness)
length      = 59*mm                 # [m]        Rod length
k_rod 	    = 4.5	                # [W/mK]     Internal heat conductivity of rod
k_air	    = 0.0265                # [W/mK]     Conductivity of air
CTE_rod     = 8.00E-06	            # [1/K]      Expansion coefficient
T_RH        = 70	                # [°C]       Rod holder temp, assumed uniform
F_spring    = -6                    # [N]        Spring force 
X_spring    = 13.74*mm              # {m}        Distance of springs from edges
Q_heat      = P_opt*(1-Eff_opt)     # [w]        Heat load, due to conversion lossses etc. 
Q_heat      = 16                    # [w]        Heat load, due to conversion lossses etc. 

E           = 3E+11                 # [Pa]       Modulus
roughnessGap= 3e-6                  # [m]        Roughness gap between rod and rod holder

avgRayAngle = 23                    # [deg]      Average angle rays entering the rod 
refl_eff    = .9                    # [-]        Reflection efficiency, assumed to be constant for all wavelengths
lightsource = 'padbar'              # [-]        Type of LED lightsource used, can be 'CSP' or 'padbar'
frac_body   = .11                   # [-]        Fraction of heat load assigned as uniform body heat load
R_corr      = 1.27                  # [-]        Radius correction factor  1 corresponds to toponly heat lod

## Analysis settings
Horizontal_transfer = True          #            Perform thermal analysis with or without lateral heat transfer
num_el_x    = 200                   #            Discretisation in x direction
num_el_y    = 10                    #            Discretisation in y direction
num_el_z    = 10                    #            Discretisation in z direction 
heatload    = 'toponly'        #            Can be 'distribution' or 'toponly', with toponly the fraction body load is also applied to the top 


Q_body      = Q_heat*frac_body      # [W]        Heat load applied as body load
Q_distr     = Q_heat*(1-frac_body)  # [W]        Heat load applied distributed according to absorbtion spectrum

if heatload == 'toponly':    # only apply curvature correction factor if heatload is distributed over rod
    R_corr = 1

# Create data structure containing a finite element representation of the rod
rod = Rod(length, height, width, k_rod, CTE_rod,E, num_el_x, num_el_y, num_el_z)

discretise(rod)

# Displacement assuming the rod is free to bend under the first assumption thermal load 
heat_only_disp = free_rod(rod, Q_heat, R_corr)

F_min_contact = minContactForce(heat_only_disp, X_spring, rod)

# Effects of springs applied to the rod on its resulting shape, assuming contact with the rod hoder at each tip and possibly in the middle
# force_posF_posC stores development of reaction force location + magnitude over iteration proedure
#force_posF_posC = findReactionForceV2(rod, X_spring, F_spring, F_min_contact, heat_only_disp)

force_posF_posC = findReactionForceV2(rod, X_spring, F_spring, F_min_contact, heat_only_disp)

x_positions, conj_moment, heatDeflection, totalDisplacement = rod.displacements()

# Optical calculations: stokes losses as function of position given absorbtion and LED emittance spectra
absorbed_power = absorbedPower(rod,lightsource, avgRayAngle, refl_eff, plot = False)

# Thermal solution 
T_2D, q_out,L,u = thermalSolution2D(rod, Q_body, Q_distr, absorbed_power, T_RH, k_air, roughnessGap, totalDisplacement, heatload, Horizontal_transfer)


#output and results processing
rod.analyse_beam_results()
#rod.plotDisplacements(force_posF_posC[-1])
plotDisplacements(rod,force_posF_posC[-1], X_spring)

plotTempDistribution(rod,x_positions,T_2D)

'''
plt.figure(6,figsize = (20,6))
plt.plot(x_positions, T_1D[-1,:], label = 'Vertical heat transfer only')
plt.plot(x_positions, T_2D[-1,:], label = 'Vertical + horizontal heat transfer')
plt.legend()
'''

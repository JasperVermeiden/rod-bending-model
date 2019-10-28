# -*- coding: utf-8 -*-
"""
Created on Tue Sep 10 10:29:27 2019

@author: J.C.Vermeiden
All functions required to analyse the sturctural loads and deformations of the loaded rod
"""

import numpy as np
from my_functions        import *

# Computes temperatures on top and bottom surfaces of unconstrained rod based on pure circular bending
# Based on this simplified temperature profile deformations, moments etc are determined 
def free_rod(rod, Q_heat, R_corr):
    # Temperature difference over the height of the rod, required to transport all incoming heat
    dT_rod =    Q_heat* rod.height/(rod.k* rod.length* rod.width)
        
    # Temperature induced rod bending
    alpha = rod.length* rod.CTE*dT_rod/rod.height   # [rad] Angle over which the entire rod is bent due to thermal expansion 
    radius= rod.length/alpha                        # [m]   Resulting bend radius due to thermal expansion
    radius= radius*R_corr                           # Correct radius of curvature for non_uniform thermal expansion given heat load case
    
    h     = radius - radius*np.cos(np.arcsin((rod.length/2)/radius))
    
    # local deflection for each element location, assuming beam tips are fixed (jointed) i.e. central section is bowing upwards
    beta_0          = np.arcsin(0.5*rod.length/radius)                          # local rotation anle with respect to symmetry axis, point on curve with thermally induced radius
    
    for element in (rod.elements):
        beta        = np.arcsin((element.get_xpos()-0.5*rod.length)/radius)     # local rotation anle with respect to symmetry axis, point on curve with thermally induced radius
        deflection  = radius*(np.cos(beta) - np.cos(beta_0))                   # Last term is compensation for the fixaction points, deflection is upwards
        element.heatDeflection = deflection
       
    print("radius of curvature: \t \t", radius, "[m]")
    print("central air gap, without spring : \t", h*1000, "[mm]")
    return h
    
# Structural effects of spring force on displacement etc. as function of spring force and location
# assuming 2 springs placed at distance X spring from either end of the rod. 
def spring_load(rod, X_spring, F_spring, F_reaction, elementID_reaction):
    rod.clearElements()
    cell_width = rod.cellWidth()
    elementID_spring = round(X_spring/cell_width)
    reactionForce = F_spring + F_reaction
    
    # Check wether the reaction force acts on the centre or on two distinct points, if the force is central its magnitude will be doubled 
    if elementID_reaction == (len(rod.elements)-elementID_reaction+1):
        F_reaction = 2*F_reaction

    # assign forces to correct elements
    rod.elements[0].force, rod.elements[-1].force = -reactionForce , - reactionForce 
    rod.elements[elementID_reaction-1].force, rod.elements[-elementID_reaction].force = F_reaction , F_reaction
    rod.elements[elementID_spring].force, rod.elements[-elementID_spring-1].force = F_spring , F_spring
    # Shear force and moment in the rod due to the spring loading
    
    rod.elements[-1].shear = rod.elements[-1].force*-1
    for i, element in enumerate (rod.elements) : 
        for prevElement in rod.elements[0:i+1]:
            element.shear += prevElement.force
        for adjElement in rod.elements[i:]:
            arm = adjElement.get_xpos()- element.get_xpos()
            element.moment += adjElement.force*arm
        element.conj_force = (element.moment/(rod.E*rod.I_x()))
    
    conj_reaction_force  = 0
    ## conjugate beam analysis
        # total conjugate reaction force, combined over the two joints (left and right contact point)
    for element in (rod.elements):        
        conj_reaction_force += element.conj_force*cell_width      
    
    # Conjugate beam, loaded with moment distribution of normal beam as force distribution. 
    # Calculate conjugate shear at each position of the beam 
    for i, element in enumerate (rod.elements[0:]) : 
        for adj_element in (rod.elements[0:i+1]):
            element.conj_shear += adj_element.conj_force*cell_width     
        element.conj_shear += -conj_reaction_force/2

    # Conjugate Moment = conjugate shear* arm to pivot integrated form location of interest to pivot location  
    for i, element in enumerate (rod.elements) : 
        element.conj_moment = 0
        # integrate over all subsequent elements
        for adj_element in (rod.elements[i:]):
            arm = (adj_element.get_xpos()- element.get_xpos())
            element.conj_moment += arm*adj_element.conj_force*cell_width
        element.conj_moment += (rod.length - element.get_xpos())*.5*-conj_reaction_force 


# Determine minimum required spring force to obtain contact of the rod with the rod holder, anywhere between the two outer contact points
# 2 step approach, first a coarse estimation, iterate untill the deflection is below 0, then iterate back up with smaller steps
def minContactForce(heat_only_disp, X_spring, rod):
    dF = 0.1
    F_spring = 0
    while rod.max_min_deflection()[2] > 0:
        F_spring -= dF
        spring_load(rod, X_spring, F_spring, 0, round(len(rod.elements)/2))

    dF = 0.001
    while rod.max_min_deflection()[2] < 0:
        F_spring += dF
        spring_load(rod, X_spring, F_spring, 0, round(len(rod.elements)/2))
    print('Force required to make rod contact RH: \t', F_spring , '[N]' )
    return F_spring


# Determine resulting contact force given spring force and thermal deformation
# Assumed is that the contact force always acts on two points equidistant from the centre (possibly in the centre)
# If the spring force is high enough to ensure contact the  alogrithm applies a reaction force to the centre of the rod. 
# While the lowest point lies further form the centre than the applied reaction force the reaction force is shifted away from the centre in 1 element increments

def findReactionForceV2(rod, X_spring, F_spring, F_min_contact, heat_only_disp):
      
    # Get baseline deformation
    spring_load(rod, X_spring, F_spring, 0, round(len(rod.elements)/2))
    
    force_posF_posC = np.array([]).reshape((0,3))
    
    if abs(F_min_contact) >= abs(F_spring):     
        print('Applied spring force not high enough to create contact between rod and rod holder')

    else: 
        print('Center of rod contacts rod holder, contact force and position are being calculated')
        # Force acts on the lowest point, from init
        F_contact_guess = 0
        elementID_reaction = round(len(rod.elements)/2) 
        F_contact = ReactionForceMagnitude(rod, X_spring, F_spring, F_contact_guess, elementID_reaction)
        force_posF_posC = np.vstack((force_posF_posC,[F_contact,elementID_reaction, rod.max_min_deflection()[3]/rod.cellWidth() ])) 

        # While the contact point is not the point where load is applied, move the assumed contact point to the left by one element
        # Use previous contact force as initial contact force guess
        while force_posF_posC[-1,1] >= force_posF_posC[-1,2]:
            elementID_reaction -= 1
            if elementID_reaction == 3: break 
            F_contact_guess = F_contact
            F_contact = ReactionForceMagnitude(rod, X_spring, F_spring, F_contact_guess, elementID_reaction)
            force_posF_posC = np.vstack((force_posF_posC,[F_contact,elementID_reaction, rod.max_min_deflection()[3]/rod.cellWidth() ]))  
                
    print('Surface reaction force [N]: \t', F_contact,'\n','Location of reaction force (from left corner) [mm]:', rod.max_min_deflection()[3]*1000 )
    return force_posF_posC

# Determine resulting contact force given spring force and thermal deformation
# Assumed is that the contact force always acts on two points equidistant from the centre (possibly in the centre)
# If the spring force is high enough to ensure contact the  alogrithm applies a reaction force to the centre of the rod. 
# Resulting average displacement between the springs is evaluated, case where resulting average displacement is lowest is selected. 
def findReactionForceV3(rod, X_spring, F_spring, F_min_contact, heat_only_disp):
      
    # Get baseline deformation
    spring_load(rod, X_spring, F_spring, 0, round(len(rod.elements)/2))
    
    force_posF_posC = np.array([]).reshape((0,3))
    
    if abs(F_min_contact) >= abs(F_spring):     
        print('Applied spring force not high enough to create contact between rod and rod holder')

    else: 
        print('Center of rod contacts rod holder, contact force and position are being calculated')
        # Force acts on the lowest point, from init
        F_contact_guess = 0
        
        # Create a range of 5 equally spaced points beteen spring and rod centre to get a rough idea of optimal spring location
        elID_spring    = round(X_spring/rod.cellWidth())
        elID_rodCentre = round(rod.num_elements()/2)

        meanDisp = []                   # Mean displacement between springs
        meanDispFlat = []               # Mean displacement between lowest points on rod
        elID = elID_spring + 2          # ELement where the RH reaction force is applied
        F_contact = 0                   # Initial guess for contact force
        while elID < elID_rodCentre:    # Start evaluating at the spring location, moving towards the rod centre
            F_contact = ReactionForceMagnitude(rod, X_spring, F_spring, F_contact_guess, elID)
            force_posF_posC = np.vstack((force_posF_posC,[F_contact,elID, rod.max_min_deflection()[3]/rod.cellWidth() ])) 
            springSection = rod.displacements()[-1][elID_spring:-elID_spring]
            centralSection= rod.displacements()[-1][elID:-elID]
            meanDisp.append(np.sum(abs(springSection))/len(springSection)) 
            meanDispFlat.append(np.sum(abs(centralSection))/len(centralSection)) 
            F_contact_guess = F_contact # Start next computation with contact force of previous iteration
            elID +=1
        elID = np.argmin(meanDisp)+2
        F_contact = force_posF_posC[elID-2,0]
        elID += elID_spring    
        print('Last computation\n element ID:',elID, 'force', F_contact)
        spring_load(rod, X_spring, F_spring, F_contact, elID)
        force_posF_posC = np.vstack((force_posF_posC,[F_contact,elID, rod.max_min_deflection()[3]/rod.cellWidth() ])) 

    return force_posF_posC, meanDisp, meanDispFlat
        

def ReactionForceMagnitude(rod, X_spring, F_spring,F_contact, elementID_reaction):
        while rod.max_min_deflection()[2] < 0:
            F_contact += 0.1        
            spring_load(rod, X_spring, F_spring, F_contact, elementID_reaction)
            #print(' element ID reaction force: \t', elementID_reaction, '\n contact force: \t \t', F_contact, '\n min deflection: \t \t', rod.max_min_deflection()[2])

        # Iterate back upwards, gives extra digit of accuracy
        while rod.max_min_deflection()[2] > 0:
            F_contact -= 0.01   
            spring_load(rod, X_spring, F_spring, F_contact, elementID_reaction)
            #print(' element ID reaction force: \t', elementID_reaction, '\n contact force: \t \t', F_contact, '\n min deflection: \t \t', rod.max_min_deflection()[2])

        # And back down again, gives another extra digit of accuracy     
        while rod.max_min_deflection()[2] < 0:
            F_contact += 0.001       
            spring_load(rod, X_spring, F_spring, F_contact, elementID_reaction)

        return F_contact

    # Same as reaction force Magnitude, only with an extra digit of precision for reaction force magnitude
    # Fixes negative displacements
def ReactionForceMagnitudeXP(rod, X_spring, F_spring,F_contact, elementID_reaction):
        while rod.max_min_deflection()[2] < 0:
            F_contact += 0.1        
            spring_load(rod, X_spring, F_spring, F_contact, elementID_reaction)
            #print(' element ID reaction force: \t', elementID_reaction, '\n contact force: \t \t', F_contact, '\n min deflection: \t \t', rod.max_min_deflection()[2])

        # Iterate back upwards, gives extra digit of accuracy
        while rod.max_min_deflection()[2] > 0:
            F_contact -= 0.0005   
            spring_load(rod, X_spring, F_spring, F_contact, elementID_reaction)
            #print(' element ID reaction force: \t', elementID_reaction, '\n contact force: \t \t', F_contact, '\n min deflection: \t \t', rod.max_min_deflection()[2])
        return F_contact

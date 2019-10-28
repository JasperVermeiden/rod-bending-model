# -*- coding: utf-8 -*-
"""
Created on Fri Sep 13 15:56:34 2019

@author: 70238599
Finite element lumped mass model using state space solution format
Still assuming uniform heat input at top surface
"""
import numpy as np

def thermalSolution2D(rod,Q_body, Q_distr, absorbedPower, T_RH, k_air, roughnessGap, gap, heatload, Horizontal_transfer):

    #Discretisation 
    num_el_x = rod.num_el_x+1             # number of elements over the lenght of the rod
    num_el_y = rod.num_el_y             # nuber of elements over the thickess (height) of the rod
    
    # Specify thermal distance gap as number of times the specified roughness gap
    gap = gap + roughnessGap
    gap_factor = np.divide(gap,roughnessGap)
   
    ## Thermal model
    element_height= rod.height/num_el_y
    element_length= rod.length/num_el_x
    
    g_y     = (element_length*rod.width)*rod.k/element_height    # [W/K]      Conductivity of rod material per element thickness, for area of element 
    g_x     = (element_height*rod.width)*rod.k/element_length
    g_gap   = (element_length*rod.width)*k_air/roughnessGap      # [W/k]      Conductivity of air gap with known thickness
    
    ## Heat capacities matrix (E) 
    E = np.zeros((num_el_x*num_el_y,num_el_x*num_el_y), float)
    np.fill_diagonal(E, 1)      # set mc_p to 1 for all elements, value is irrelevant for steady state result, validated 26-09-2012
                                # Could also use np.identity(num_elements, float)
    
    ## Heat transfer matrix (K)
    K = np.zeros((num_el_x*num_el_y,num_el_x*num_el_y), float)
    
    # Add vertical heat transfer to K matrix
    for j in range (num_el_x):
        # heat transfer to rod holder
        K[(j+1)*num_el_y-1, (j+1)*num_el_y-1] -= g_gap/gap_factor[j]
        
        # loop over all elments in vertical direction to add vertical heat transfer for each element    
        for i in range (num_el_y):
            index = j*num_el_y + i
            
            if (i-1) >= 0:
                K[index, index-1] += g_y
                K[index, index] -= g_y
                
            if (i +1) <= (num_el_y-1):    
                K[index, index+1] += g_y
                K[index, index] -= g_y
    
    if Horizontal_transfer:
        # Add horizontal heat transfer to K matrix
        for j in range (num_el_x):
            # loop over all elments in vertical direction to add horizontal heat transfer for each element    
            for i in range (num_el_y):
                index = j*num_el_y + i
                
                if (index-num_el_y) >= 0:
                    K[index, index-num_el_y] += g_x
                    K[index, index] -= g_x
                    
                if (index + num_el_y) <= (num_el_y*num_el_x-1):    
                    K[index, index+num_el_y] += g_x
                    K[index, index] -= g_x
    
    ## Heat load input vector
    if heatload == 'toponly': 
        Q_distr += Q_body
    if heatload == 'distribution':
        Q_distr = Q_distr/sum(absorbedPower)    # normalise distributed Q so that the actual distributed heat load corresponds to this value
    
    u = np.array(([T_RH], [Q_distr/num_el_x],[Q_body/(num_el_x*num_el_y)]))     # Case assuming uniform heat loading on the top surface of the rod and uniform rod hoder temperature
    print(u)
    ## Heat load matrix (L)
    L = np.zeros((num_el_x*num_el_y,np.size(u)), float)     
    
    # Heat flux to rod holder
    for j in range (num_el_x):
        index = (j+1)*num_el_y-1 
        L[index,0] = g_gap/gap_factor[j]

    if heatload == 'toponly':           # Uniform heat influx on rod top surface
        for j in range (num_el_x):
            index = j*num_el_y 
            L[index,1] = 1
        
    if heatload == 'distribution':
        for j in range (num_el_x):
            index = list(np.arange(0,num_el_y)+j*num_el_y) 
            L[index,1] = [absorbedPower] 
        # Uniform heat load (reabsorbtion etc.)
        L[:,2] = 1



    
    ## state space solution, steady state
    A = -np.linalg.inv(E)@K
    B = -np.linalg.inv(E)@L
    C = np.identity(num_el_x*num_el_y, float) # Interested in T values for all elements
    D = np.zeros(np.shape(L),float)
    x = -np.linalg.inv(A)@B@u
    T_static = C@x +D@u
    
    ## Unpack 1D temperatures vector into a matrix representing the rod geometry
    T_staticMatrix = np.zeros((num_el_y, num_el_x))
    for i in range (num_el_x):
        T_staticMatrix[:,i] = list(T_static[i*num_el_y:(i+1)*num_el_y ]) 
    q_out = np.subtract(T_staticMatrix[-1,:],T_RH)*(np.divide(g_gap,gap_factor))
    q_in  = sum(L[:,1:]@u[1:])
    return T_staticMatrix, q_out,L,u
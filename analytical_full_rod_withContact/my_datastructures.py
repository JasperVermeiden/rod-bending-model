# -*- coding: utf-8 -*-
"""
Created on Tue Sep  3 10:26:02 2019

@author: 70238599
"""

import numpy 
import math
import matplotlib.pyplot as plt

class Rod:   
    
    def __init__ (self,length, height, width, k_rod, CTE, E, num_el_x, num_el_y, num_el_z):
        self.length   = length
        self.height   = height
        self.width    = width
        self.k        = k_rod
        self.CTE      = CTE
        self.E        = E
        self.num_el_x = num_el_x
        self.num_el_y = num_el_y
        self.num_el_z = num_el_z
    
    def set_elements (self, elements):
        self.elements = elements 
    
    def cellWidth(self):
        return (self.length/(len(self.elements)-1))

    def clearElements(self):
        for element in self.elements:
            element.clear()
    
    def num_elements (self):
        return (len(self.elements) -1)
    
    # Area moment of inertia about horizontal axis
    def I_x (self):
        return self.width*math.pow(self.height,3)/12
    
    def analyse_beam_results(self):
        shear = numpy.array([])
        moment =numpy.array([])
        conj_moment = numpy.array([])
        conj_shear = numpy.array([])
        conj_force = numpy.array([])
        heatDeflection = numpy.array([])
        for element in self.elements:
            shear = numpy.append(shear,element.shear)
            conj_shear = numpy.append(conj_shear, element.conj_shear)
            moment = numpy.append(moment, element.moment)
            conj_moment = numpy.append(conj_moment, element.conj_moment)
            conj_force = numpy.append(conj_force, element.conj_force)
            heatDeflection = numpy.append(heatDeflection, element.heatDeflection)
        plt.figure(1, figsize = (12,12))
        plt.subplot(221)
        plt.plot(shear)
        plt.title('shear')
        
        plt.subplot(222)
        plt.plot(moment)
        plt.title('moment')

        plt.subplot(223)
        plt.plot(conj_shear)
        plt.title('Theta [rad] (conjugate shear)')
      
        plt.subplot(224)
        plt.plot(conj_force)
        plt.title('conjugate force')      

               
        # Reads displacement values for all the nodes and returns them in numpy array form.
        # returns:  spring resultant displacement,          [m]
        #           deflection due to thermal expansion     [m]
        #           total deflection due to combined effect [m]
        #           maximum delfection in central section,  [m]  excluding the outer 20% of each half of the rod where penetration can occur
        #           loc_max deflection [m], distance from edge to location of max displacement 
    
    def displacements(self):
        x_positions = numpy.array([])
        conj_moment = numpy.array([])
        heatDeflection = numpy.array([])
        for element in self.elements:
            element.displacement = element.conj_moment + element.heatDeflection
            conj_moment = numpy.append(conj_moment, element.conj_moment)
            heatDeflection = numpy.append(heatDeflection, element.heatDeflection)
            x_positions = numpy.append(x_positions, element.get_xpos())
        totalDisplacement = numpy.add(conj_moment, heatDeflection)
        return x_positions, conj_moment, heatDeflection, totalDisplacement

        # Find positions of maximum displacement and maximum air gap, looking at the left half of the rod only
    def max_min_deflection(self):
        x_positions, conj_moment, heatDeflection, totalDisplacement = self.displacements() 
        maxDeflection = max(totalDisplacement[round(len(self.elements)/10) :round(len(self.elements)/2)])     
        minDeflection = min(totalDisplacement[round(len(self.elements)/10) :round(len(self.elements)/2)])
        cellID_max  = numpy.where(totalDisplacement == maxDeflection)
        cellID_min  = numpy.where(totalDisplacement == minDeflection)
        loc_maxDeflection = self.elements[cellID_max[0][0]].get_xpos()
        loc_minDeflection = self.elements[cellID_min[0][0]].get_xpos()
        return maxDeflection, loc_maxDeflection, minDeflection, loc_minDeflection
              
        

        
class Element:
    
    # Set element number and x position of element 
    def __init__(self, el_numb, x_pos):
        self.__x_pos    = x_pos
        self.__el_numb  = el_numb
        self.conj_shear = 0
        self.conj_force = 0
        self.conj_moment = 0
        self.moment = 0
        self.shear = 0
        self.displacement = 0
        self.heatDeflection = 0
        self.force = 0
        self.conj_reaction_force = 0
        
    def clear(self):
        self.conj_shear = 0
        self.conj_force = 0
        self.conj_moment = 0
        self.conj_reaction_force = 0
        self.moment = 0
        self.shear = 0
        self.displacement = 0
        self.force = 0

        
        
    def get_xpos(self):
        return self.__x_pos
    
        
        
        
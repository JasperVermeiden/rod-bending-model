# -*- coding: utf-8 -*-
"""
Created on Mon Sep 16 16:29:01 2019

@author: 70238599
Absorbtion functions, calculate amount of light energy absorbed by a species as function of absorbtion spectrum, light spectrum and position in the sample
"""
import numpy as np
import matplotlib.pyplot as plt


# Fitted curve on emittance spectrum, data from 'Copy of LuAG and LuYAG material data for LT simulation 2019-09-04' accessed 16-09-2019
def LEDspectrum(lightsource):
        
    if lightsource == 'CSP':
        mu_inv  = 449.0                 # [nm]  Intensity peak location
        sigma   = 9.99e-5              
        amp	    = 0.00188               # Amplitude scaler
        c	    = 2.265                 
        alpha	= 0.142                 # Skewness factor
    
    if lightsource == 'padbar':
        mu_inv  = 446.5                 # [nm]  Intensity peak location
        sigma   = 9.89e-5              
        amp	    = 0.00199               # Amplitude scaler
        c	    = 2.615
        alpha	= 0.046                 # Skewness factor
    
    # normalised spectral intensity as function of wavelength, fitted curve
    intensities = []
    wavelength = []       # wavelength in nm
    for wl in range (200,801):
        x = (1/wl-1/mu_inv)
        intensity = amp/(np.sqrt(2*np.pi*sigma))*np.exp(-0.5*np.power(x/sigma,2))/(1+c*np.power(x/sigma,2))*0.5*(1+np.math.erf(alpha*(x)/(np.sqrt(2)*sigma)))
        intensities.append(intensity)
        wavelength.append(wl)
    return(wavelength, intensities)

def Absorbtion(wavelengths, lightSpectrum, absorbtionSpectrum, avgRayAngle, rod):
    Transmission_perc = np.array([wavelengths])
    Transmission_watt = np.array([wavelengths])
    dy = rod.height*1000/rod.num_el_y
    y_coordinates  = np.arange(0, rod.height*1000 + dy, dy)  # Y coordinate of element, list is 1 element longer than list of elements

    for depth  in y_coordinates:      
        mean_pathLength = depth*np.cos(avgRayAngle/(180*np.pi))
        Trans_perc = np.array([])
        Trans_watt = np.array([])

        for i, wl in enumerate (wavelengths): 
            # percentage of light that is transmitted, as function of wavelength
            Trans_perc = np.append(Trans_perc, np.exp(-absorbtionSpectrum[i,2]*mean_pathLength))
            Trans_watt = np.append(Trans_watt, np.exp(-absorbtionSpectrum[i,2]*mean_pathLength)*lightSpectrum[i])
        Transmission_perc = np.vstack((Transmission_perc, Trans_perc))        
        Transmission_watt = np.vstack((Transmission_watt, Trans_watt))    
    
    # Transmitted light flux [w] at each cell face (normalised to 1 watt at the entrance surface of the rod)
    P_light = sum(Transmission_watt.T[:,1:])
    
    absorbed_power = []
    # Absorbed light power per cell
    for i in range (0, len(P_light)-1):
        absorbed_power.append(P_light[i]- P_light[i+1])

    transmitted_power = P_light     
    
    return absorbed_power, transmitted_power, Transmission_watt[-1,:]


# Returns absorbed power per cell, discretised along the height of the rod with number of elements specified in rod.num_el_y
def absorbedPower(rod,lightsource, avgRayAngle, refl_eff, plot):
    # Load rod absorbance spectra data for the 7 tested rods.  
    absorbtionSpectrum = np.load('spectral_data/absorbancesV2.npy')
    outputSpectrumHLD  = np.load('spectral_data/00098_no_reflector_18mm_PROCESSED.npy')
    outputSpectrumHLD  = np.vstack((outputSpectrumHLD[0,1:],np.mean(outputSpectrumHLD[1:,1:], axis = 0)))
    wavelengths, lightSpectrum = LEDspectrum(lightsource)    
    
    # Calculate: absorbed power for each cell, transmitted power at each cell interface & exiting light spectrum
    absorbed_power, transmitted_power, exitSpectrum = Absorbtion(wavelengths, lightSpectrum, absorbtionSpectrum, avgRayAngle, rod)
    
    lightSpectrum2 = np.multiply(exitSpectrum, refl_eff)
    
    # same for the first reflection
    absorbed_power2, transmitted_power2, exitSpectrum2= Absorbtion(wavelengths, lightSpectrum2, absorbtionSpectrum, avgRayAngle, rod)
    
    # Filp the datasets since the reflection passes though the rod in opposite direction
    absorbed_power2     = np.flipud(absorbed_power2)
    transmitted_power2  = np.flipud(transmitted_power2)
    
    # Spectrum averaged wavelength, input and output spectra -> stokes efficiency
    WLmeanOutput    = np.sum(np.multiply(outputSpectrumHLD[0,:], outputSpectrumHLD[1,:]))/sum(outputSpectrumHLD[1,:])
    WLmeanInput     = np.sum(np.multiply(wavelengths, lightSpectrum))/sum(lightSpectrum)
    Stokes_eff      = (1/WLmeanOutput)/(1/WLmeanInput)
    
    if plot:
        plt.figure(10, figsize = (10,6))
        plt.plot(wavelengths, lightSpectrum, label = 'Source LED spectrum')
        plt.plot(wavelengths, exitSpectrum, label = 'Exit spectrum after first pass through rod')
        plt.plot(wavelengths, lightSpectrum2, label = 'Reflected LED spectrum')
        plt.plot(wavelengths, exitSpectrum2, label = 'Exit spectrum after second pass through rod')
        plt.legend()
        plt.title('Light intensities, entrance spectrum normalised to 1W integrated')
        
        y_values = np.arange(rod.height/(2*rod.num_el_y),rod.height +rod.height/(2*rod.num_el_y), rod.height/rod.num_el_y)
        plt.figure(11,figsize = (10,6))
        plt.plot(y_values, absorbed_power, '*-', label ='First pass')
        plt.plot(y_values, absorbed_power2, '*-', label ='Second pass')
        plt.plot(y_values, np.add(absorbed_power,absorbed_power2), '*-', label = 'Summed first & second pass')
        plt.plot(rod.height, sum(absorbed_power), '*', label ='cumulative power absorbed on first pass' )
        plt.plot(0, sum(absorbed_power2), '*', label = 'Cumulative power absorbed on second pass')
        plt.title('Absorbed power per cell [W], spectral power normalised to 1W')
        plt.grid('True')
        plt.legend()
    
    print('Absorbed power first pass:\t \t \t', sum(absorbed_power))
    print('Absorbed power first + second pass:\t \t', sum(absorbed_power)+ sum(absorbed_power2))
    print('Percentage extra absorbtion with reflector:\t', (100*sum(absorbed_power2))/sum(absorbed_power))
    print('stokes efficiency for light spectrum with absorbtion spectrum:\t', Stokes_eff)
    
    return np.add(absorbed_power,absorbed_power2)
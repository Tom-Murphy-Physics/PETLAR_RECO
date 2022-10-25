"""
This code simulates electron drift and readout. 
The pixel plane geometry is using the so called square formation.
It is characterized by a 8x8 grid of pixels.
"""

import numpy as np
from numpy import random
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import math
from matplotlib.collections import PolyCollection
import Drift_Axial
import Generate

#   Function that returns the drift time of ionized electrons
#   The 100 assumes a 100 mm drift length
#   100 was determined by looking at the results from the GEANT4 code
def t(pendepth):
    return (100-pendepth)/V

#   Noise is the number of electrons below which the signal is swamped by noise
#   Another collaboration found this number so we are using it as an estimate
noise = 87

#   Number of ionized electrons (PETLAr GEANT4)
#   The denominator is used to scale the energy of the photon 
#   Use 600 to test 100 keV
#   Use 600 to test 511 keV
#   I was looking at sensitivity as it pretained to CT scans
#N0 = int(1500)

#   This function will convert a number of pixels into an ADC count
#   If below noise, say it is zero
#   Assumes an 8 bit ADC that can go upto 4096 electrons
def ADC(x):
    if math.isnan(x)==False:
        if x < noise:
            return 0
        else:
            x_ = np.floor(x)
            nume_diff = 16
            ADC_value = 0
            while x_>0:
                x_ = x_ - nume_diff
                ADC_value += 1
            return ADC_value

#   Actually simulates the pixels 
def run(pixel_size,target_r,plot,N0):
       #   Defines a vector containing all the verticies 
    verticies = np.empty(shape=(64,4,2),dtype='float')
    
    #   Defines the plots if required
    if plot == True:
        fig2, ax2 = plt.subplots()
        ax2.set_facecolor('0.5')
        ax2.set_aspect('equal', adjustable='box')
    index = 0
    
    #   Defines the verticies and draws the plots if applicable
    #   The for loops and if statments help define the packed geometry
    for i in range(8):
        for j in range(8):
            if plot == True:
                ax2.add_patch(Rectangle((-i*pixel_size+3*pixel_size,-j*pixel_size+3*pixel_size), pixel_size, pixel_size, color = "white", ec="black"))
            verticies[index] = np.array([[-i*pixel_size+3.*pixel_size,-j*pixel_size+3.*pixel_size],[-i*pixel_size+3.*pixel_size,-j*pixel_size+4.*pixel_size],[-i*pixel_size+4.*pixel_size,-j*pixel_size+4.*pixel_size],[-i*pixel_size+4.*pixel_size,-j*pixel_size+3.*pixel_size]])
            index += 1
                 
    #   r0 is where, in the transverse direction the photon Compton scatters
    #   Uniform circular distribution with maximum radius r0
    r0 = random.uniform(0,target_r)
    #r0 = 0
    
    #   Picks an angle in the x-y plane
    theta0 = random.uniform(0,2*np.pi)
    
    #   Converts to x-y plane 
    x0 = r0*np.cos(theta0)
    y0 = r0*np.sin(theta0)
    
    #   Penetration depth of photon
    pd = random.normal(40, 5)
        
    #   Creates a distribution of electrons
    #   Distribution is governed by diffusion
    generated_z = (Generate.generate()[1])/1000
    generated_r = Generate.generate()[0]/1000
    xdata = random.normal(x0,Drift_Axial.stdev(generated_z),N0)
    ydata = random.normal(y0,Drift_Axial.stdev(generated_z),N0)
        
    #   Used in comparison to Scintillation photons
    #xdata = np.ones(N0)*x0
    #ydata = np.ones(N0)*y0
    
    #   e_on_p = Number of electrons on Pixel
    #   each index represents a different pixel
    #   If no electrons you would get 64 zeros
    e_on_p = np.zeros(64)
    
    i = 0
    counter = 0
    
    #   Counts the number of electrons on each pixel and stores it in the e_on_p array
    while i < N0:
        j = 0
        while j < 64:
            if xdata[i]>verticies[j][0][0]:               
                if xdata[i]<verticies[j][2][0]: 
                    if ydata[i]<verticies[j][2][1]: 
                        if ydata[i]>verticies[j][0][1]:
                            e_on_p[j] = e_on_p[j]+1
            j += 1
        i += 1
        
    #   num_adc = Total number of ADC counts registered
    num_adc = 0
    
    #   List of values from reconstruction
    xreco = []
    yreco = []
    
    #   max_adc = ADC counts on the pixel with the most electrons
    max_adc = 0
        
    #   Reconstructs the x-y position of Compton scattering
    #   Adds the number of ADC counts to the plane 
    #   Adds a color coding to the pixel based on ADC counts
    for i in range(64):
        num_adc = num_adc + ADC(e_on_p[i])
        #   Finds center of pixel
        xcent = (float(verticies[i][2][0]) + float(verticies[i][0][0]))/2
        ycent = (float(verticies[i][1][1]) + float(verticies[i][3][1]))/2
        #   Plots things if required
        if plot == True:
            ax2.add_patch(Rectangle((xcent-.5*pixel_size,ycent-.5*pixel_size), pixel_size, pixel_size, color = str(1-ADC(e_on_p[i])/128), ec="black"))
            ax2.text(xcent, ycent, ADC(e_on_p[i]), ha="center", va="center", color = "black", fontweight = "bold")
            #if e_on_p[i] > 0:
             #   ax2.text(xcent, ycent, ADC(e_on_p[i]), ha="center", va="center", color = "black", fontweight = "bold")
            
            Truth = plt.Circle((x0,y0), 0.05*pixel_size, color = "red")
            ax2.add_patch(Truth)
                           
             
        #   Used in the weighted average (reconstruction)
        xreco+=[ADC(e_on_p[i])*xcent]
        yreco+=[ADC(e_on_p[i])*ycent]
        #   Finds the x-y position by just looking at the maximum
        if ADC(e_on_p[i])>max_adc:
            max = ADC(e_on_p[i])
            x1 = xcent
            y1 = ycent
    
    if plot == True:
        k = 0
        while k < N0:
            
            data_point = plt.Circle((xdata[k],ydata[k]), 0.01 * pixel_size, color = "grey", alpha = 0.25)
            ax2.add_patch(data_point)
            k = k + 1
    
    # If counts were recorded do weighted average
    if num_adc > 0:
        x2 = sum(xreco)/num_adc
        y2 = sum(yreco)/num_adc
    # Assigns a value, taken care of at the end (returns "None")
    else:
        x2 = x0
        y2 = y0
        x1 = x0
        y1 = y0
        
    if plot == True:
        #   values is used if you want to produce a plot with the pixel's index printed on it
        values=np.arange(64)
        #   Create collection of rectangles.
        pc = PolyCollection(verticies, closed=True, edgecolors="k", linewidth=0.72, cmap="inferno")
        pc.set_array(values)
        ax2.add_collection(pc)
        Reco = plt.Circle((x2,y2), 0.05*pixel_size, color = "black")
        ax2.add_patch(Reco)
        plt.xlim(-3,3)
        plt.ylim(-3,3)
        plt.savefig("square_sample.png", dpi = 2000)
        plt.show()
    
    #   If ADC counts were recorded return values otherwise return "None"
    if num_adc>0:
        return x0-x1, y0-y1, x0-x2, y0-y2, generated_r, generated_z
    else:
        return None, None, None, None

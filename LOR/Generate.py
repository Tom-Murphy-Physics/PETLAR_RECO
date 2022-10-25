"""
Written by Tom Murphy
This code uses data generated from GEANT4 simulation
It produces the position of a compton scatter 
Note that this is an approixmation
More accurate infomration can be found using the GEANT code
"""

import numpy as np
import matplotlib.pyplot as plt

#   Z = height of detector
#   R0 = Inner cylinder's radius
#   R1 = Outter cylinder's radius
Z = 200
R0 = 50
R1 = 150

#   Used in 2D histogram 
#   N = number of iteractions
N = 100000

#   This is the maximum angle that the photon could travel
#   Anything more and it would not be in the fiducial volume
Theta_lim = np.arctan(Z/R0)

#   Generates the interaction point
def generate():
    valid = False
    while valid == False:
        #   Produces an angle
        theta = np.random.uniform(-Theta_lim, Theta_lim)
        #   Produces a penetration depth
        PD = np.random.normal(40,6)
        #   Converts to cylindrical coordinates
        R = 50 + PD*np.cos(theta)
        z = R*np.tan(theta)
        #   Checks to see if the interactions is inside the detector
        if R <= R1:
            if abs(z) <= Z:
                valid = True
    return R,z

#   This code makes the 2D histogram of interaction points
"""
data = [generate() for i in range(N)]

x = []
y = []

for i in range(N):
    #plt.scatter(data[i][0],data[i][1],color = "black", alpha = 0.2)
    x += [data[i][0]]
    x += [-data[i][0]]
    y += [data[i][1]]
    y += [-data[i][1]]

x_min = -R1
x_max = R1
  
y_min = -Z
y_max = Z
  
x_bins = np.linspace(x_min, x_max, 100)
y_bins = np.linspace(y_min, y_max, 100)
  
fig, ax = plt.subplots(figsize =(10, 7))
# Creating plot
plt.hist2d(x, y, bins =[x_bins, y_bins])
plt.title("Distribution of 1st Compton Scatter")
  
ax.set_xlabel('R-axis') 
ax.set_ylabel('Z-axis') 
  
# show plot
plt.colorbar()
plt.tight_layout() 
plt.show()
"""
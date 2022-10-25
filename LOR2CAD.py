"""
Written by Tom Murphy
This code turns LOR to a CAD model by:
    Converting the LOR to a sinogram 
    Preforming an inverse radon transform (Reconstruction)
    Creating 2D projections of the reconstruction
    Transforming the projections to a 3D map
    Filtering the data
    Writing the filtered data to an stl file
"""
from skimage.transform import iradon
import numpy as np
import matplotlib.pyplot as plt
from skimage.data import shepp_logan_phantom
from skimage.transform import radon, rescale
import skimage.io
from PIL import Image
import matplotlib.patheffects as path_effects
from mpl_toolkits.axes_grid1 import AxesGrid
import matplotlib as m
import to_stl
import sys

#   Allows you to use code from the "LOR" directory
#   That code will create the "Lines of Response" produced by PET data
sys.path.insert(0, './../LOR')

#   Imports the "LOR.py" File from ./../LOR/
import LOR

#   Generates N lines of response 
N = 100
LOR_data = LOR.Generate_LOR(N)

#   Defines the number of pixels in the geometry 
#   I still need to verify this is the right number of pixels
pixel_density = 101

#   Defines the detector length
L_det = 100

#   xs and ys are arrays containing 2 numbers
#   This function finds the equation of the line that passes through both points
#   In cartesian coordinates
def find_line(xs, ys):
    m = (ys[1]-ys[0])/(xs[1]-xs[0])
    b = ys[0] - m*xs[0]
    return m, b

#   This converts the line into sinogram coordinates 
def convert(m,b,x_p,y_p,x0):
    if m != 0:
        m_ = -1/m
    else:
        m_ = 1000000
    b_ = y_p - m_ * (x_p - x0)
    x_i = (b_ - b)/(m - m_) + x0
    y_i = m_ * (x_i - x0) + b_
    distance = -y_i/abs(y_i)*np.sqrt((x_i - x0 - x_p)**2 + (y_i - y_p)**2)
    angle = np.arctan((x_i - x0)/y_i)*180/np.pi
    return distance, angle

#   Decleration of data arrays 
data_xy = []
x_xy = []
y_xy = []

data_yz = []
y_yz = []
z_yz = []

data_xz = []
x_xz = []
z_xz = []

#   These help define a cube and a sphere
L = 50
R = 20



#   Uses the LOR raw data to generate a sinogram
for i in range(N):
    line_xy = find_line(LOR_data[0][i],LOR_data[1][i])
    data_xy += [convert(line_xy[0],line_xy[1],0,0,0)]
    x_xy += [data_xy[i][0]]
    y_xy += [data_xy[i][1]]

    line_yz = find_line(LOR_data[2][i],LOR_data[1][i])
    data_yz += [convert(line_yz[0],line_yz[1],0,0,0)]
    y_yz += [data_yz[i][0]]
    z_yz += [data_yz[i][1]]
    
    line_xz = find_line(LOR_data[2][i],LOR_data[0][i])
    data_xz += [convert(line_xz[0],line_xz[1],0,0,0)]
    x_xz += [data_xz[i][0]]
    z_xz += [data_xz[i][1]]
        
#   Generates the Sinogram (represented as a 2D histogram)
#   Three projections onto the xy, yz, xz planes
#   Bins are the pixels
#   Horizontal axis is the s coordinates
#   Vertical Coordinate is the angle
hist_xy, xbins_xy, ybins_xy = np.histogram2d(x=x_xy, y=y_xy, bins=(np.arange(-L_det,L_det,2*L_det/pixel_density), np.arange(-90, 90, 180/pixel_density)))
hist_yz, ybins_yz, zbins_yz = np.histogram2d(x=y_yz, y=z_yz, bins=(np.arange(-L_det,L_det,2*L_det/pixel_density), np.arange(-90, 90, 180/pixel_density)))
hist_xz, xbins_xz, zbins_xz = np.histogram2d(x=x_xz, y=z_xz, bins=(np.arange(-L_det,L_det,2*L_det/pixel_density), np.arange(-90, 90, 180/pixel_density)))

#   Defines figure 
fig, (ax1, ax2, ax3) = plt.subplots(1,3)

#   Adds Signograms to figure
ax1.imshow(hist_xy)
ax2.imshow(hist_yz)
ax3.imshow(hist_xz)

#   Adding some labels
ax1.set_xlabel("XY Plane")
ax2.set_xlabel("YZ Plane")
ax3.set_xlabel("XZ Plane")

ax2.set_title("Sinogram Projections onto Planes")


####################
##                ##
##      RECO      ##
##                ##
####################

plt.show()

#   Defines figure for the reconstruction projections
fig, (ax1, ax2, ax3) = plt.subplots(1,3)

#   iradon is a function that preforms an inverse radon transform
#   iradon converts s, theta into cartesian coordinates
recon_xy = np.flip(iradon(hist_xy,filter_name='ramp'))
recon_xy = [np.flip(i) for i in recon_xy]
ax1.imshow((recon_xy), cmap=plt.cm.Greys_r)

recon_yz = np.flip(iradon(hist_yz,filter_name='ramp'))
recon_yz = [np.flip(i) for i in recon_yz]
ax2.imshow((recon_yz), cmap=plt.cm.Greys_r)

recon_xz = np.flip(iradon(hist_xz,filter_name='ramp'))
recon_xz = [np.flip(i) for i in recon_xz]
ax3.imshow((recon_xz), cmap=plt.cm.Greys_r)

#   Adding some labels
ax1.set_xlabel("XY Plane")
ax2.set_xlabel("YZ Plane")
ax3.set_xlabel("XZ Plane")

ax2.set_title("Reconstruction Projections")

plt.show()

####################
##                ##
##    3D-itize    ##
##                ##
####################

#   Defines a data format
#   Basically tells you the value at (x[i],y[j],z[k])
data = np.zeros(shape = (pixel_density,pixel_density,pixel_density))

#   Takes the Reconstructed 2D histograms and just copies them 
#   For the xy plot it makes the number of copies into the z direction = number of pixels
#   Same data structure as "data"
#   Simillarly for the other directions
dataxyz = [recon_xy for i in range(pixel_density)]/np.amax(recon_xy)
datayzx = [recon_yz for i in range(pixel_density)]/np.amax(recon_yz)
dataxzy = [recon_xz for i in range(pixel_density)]/np.amax(recon_xz)

#   Adds up the different data_ _ _
for i in range(pixel_density-1):
    j = 0
    for j in range(pixel_density-1):
        k = 0
        for k in range(pixel_density-1):
            data[i][j][k] =  dataxzy[pixel_density -2 - j][pixel_density -2 - i][k]   +   dataxyz[pixel_density - 2 - k][pixel_density -2 - j][i]  +  datayzx[i][pixel_density - 2 -j][k]

    #   FILTERING OF DATA   #

#   Sets up the data structure for filtered data       
c = []
x_pos = []
y_pos = []
z_pos = []

#   Defines th parameter space for the pixels
x = np.arange(-100,100,200/pixel_density)
y = np.arange(-100,100,200/pixel_density)
z = np.arange(-100,100,200/pixel_density)

#   Filters the data 
#   Still need to work on this a bit_length
#   RN basically says if it isn't one of the highest values throw it out 
#   Could screw up reconstruction of multiple objects 
for j in range(pixel_density-1):
    #thresh = np.quantile(abs(data), .9999)
    thresh = 1.185
    for i in range(pixel_density-1):
        for k in range(pixel_density-1):
            if data[i][j][k] >= thresh:
                if data[i][j][k] > recon_xy[i][j]+0.8 and data[i][j][k] > recon_yz[j][k]+0.8 and data[i][j][k] > recon_xz[i][k]+0.8:
                    data[i][j][k] = data[i][j][k]
                    c += [data[i][j][k]]
                    x_pos += [x[i]]
                    y_pos += [y[j]]
                    z_pos += [z[k]]
            else: 
                data[i][j][k] = 0

#   Plots different slices of the 3D reconstruciton
"""
fig, axs = plt.subplots(2,3)   
vals = [data[0],data[10],data[20],data[40],data[50],data[60],data[70],data[80]]
for i in vals:
    max = 0
    if np.amax(i) > max:
        max = np.amax(i)
        im = axs[1][2].imshow(i, interpolation='nearest')

axs[0][0].imshow(data[0], interpolation='nearest')
axs[0][1].imshow(data[0], interpolation='nearest')
axs[0][2].imshow(data[5], interpolation='nearest')
axs[1][0].imshow(data[10], interpolation='nearest')
axs[1][1].imshow(data[15], interpolation='nearest')
axs[1][2].imshow(data[20], interpolation='nearest')
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
fig.colorbar(im, cax=cbar_ax)
plt.show()
print(im)
"""

#   Makes a 3D scatter plot of the filtered data
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
img = ax.scatter(x_pos, y_pos, z_pos, c=c, cmap=plt.hot())
ax.axes.set_xlim3d(left=-100, right=100) 
ax.axes.set_ylim3d(bottom=-100, top=100) 
ax.axes.set_zlim3d(bottom=-100, top=100) 
fig.colorbar(img)
ax.set_box_aspect((1,1,1))
plt.show()

####################
##                ##
## STL CONVERSION ##
##                ##
####################

#   Converts the 3D scatter plot into a bunch of boxes
#   Then gives the verticies of the boxes 
#   This is used to make the stl file
cube_verticies = []
for i in range(len(x_pos)):
    cube_verticies += [[[x_pos[i]],[y_pos[i]],[z_pos[i]],[x_pos[i]+2],[y_pos[i]+2],[z_pos[i]+2]]]

#   Makes a file called "Test.stl"
#   See to_stl.py for more information about how the file is made
F = open("test.stl", "w")
to_stl.intro(F)
for i in range(len(cube_verticies)-1):
    to_stl.make_cube(F,cube_verticies[i][0][0],cube_verticies[i][1][0],cube_verticies[i][2][0],cube_verticies[i][3][0],cube_verticies[i][4][0],cube_verticies[i][5][0])
to_stl.end(F)
F.close()

#   Takes an image, generates a sinograph from it and then reconstructs image
#   Need to "3D-itize" this part of the code
"""
image = Image.open('Untitled.png')
image = image.convert('1')
image = skimage.img_as_float(image)

image = rescale(image, scale=1, mode='reflect', channel_axis=None)
print(image)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 4.5))

ax1.set_title("Original")
ax1.imshow(image, cmap=plt.cm.Greys_r)

theta = np.linspace(0., 180., max(image.shape), endpoint=False)
sinogram = radon(image, theta=theta)

print(np.amax(sinogram))

dx, dy = 0.5 * 180.0 / max(image.shape), 0.5 / sinogram.shape[0]
ax2.set_title("Radon transform\n(Sinogram)")
ax2.set_xlabel("Projection angle (deg)")
ax2.set_ylabel("Projection position (pixels)")
ax2.imshow(sinogram, cmap=plt.cm.Greys_r,
           extent=(-dx, 180.0 + dx, -dy, sinogram.shape[0] + dy),
           aspect='auto')

fig.tight_layout()
plt.show()

reconstruction_fbp = iradon(sinogram, theta=theta, filter_name='ramp')

error = reconstruction_fbp - image
print(f'FBP rms reconstruction error: {np.sqrt(np.mean(error**2)):.3g}')

imkwargs = dict(vmin=-0.2, vmax=0.2)
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 4.5),
                               sharex=True, sharey=True)
ax1.set_title("Reconstruction\nFiltered back projection")
ax1.imshow(reconstruction_fbp, cmap=plt.cm.Greys_r)
ax2.set_title("Reconstruction error\nFiltered back projection")
ax2.imshow(reconstruction_fbp - image, cmap=plt.cm.Greys_r, **imkwargs)
plt.show()
"""
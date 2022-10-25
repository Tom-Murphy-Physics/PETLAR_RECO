# PETLAR_RECO
This code is intented to simulate the reconstruction process for the PETLAr experiment.

# Dependencies
numpy
matplotlib
sklearn
PIL

## What it does
In PET scans, two back to back photons are detected with a cylindrical detector.
A line is then drawn between the two of detected photons.
This line called a "line of response" (LOR)

A LOR is then converted into a new coordinate system to store its information. 
A point is taken as the origin, then the shortest distance from the origin to the line is taken.
This distance forms the first coordinate. 
The second coordinate is the angle of the LOR from the horizantal axis. 

This is done for each set of LORs. To reduce any ambiguties, the data from each LOR is grouped into projections.
XY, YZ, XZ are used. 

Once the LOR are converted to the new coordinates, we need to reconstruct the image.
This is done using an "inverse radon transform". This code uses a pre-defined library to accomplish this. 

This will give us 3, 2D projections of the reconstructed images in the XY, YZ, and XZ planes. 
From here, a simple algorithim converts these projections into a 3D image. 
This is done by extending each 2D projection into 3D space. 
Then a filter is applied to try and extract the real positions. 
Currently this has shown mixed success and seems to be affected by the resolution of the reconstruction. 
With high statistics, it usually returns a fairly accurate result.
These data points are then converted into an .stl file.
This .stl file can then be viewed in your favorite CAD software.

## How it does it
To run the program simply go into the Sinogram directory and execute the LOR2CAD.py script. 

### LOR.py
LOR2CAD will first execute a file called LOR.py (in the LOR directory).
LOR.py will use a user defined distribution of the source.
It will use this, and data collected from the PETLAr GEANT4 simulation to determine where the compton scatter will occur.
Then, it will execute either Test_Axial.py or Test_Radial.py
This will return the error of the reconstructed photon position. 
LOR.py will return the reconstrcuted position (with errors)

### Test_Axial.py
This file will simulate the drifting of charge along the axial direction.
This field is uniform and has a strength of 1kV/cm.
The electrons behave as they would in simillar LArTPCs.
The functions for the electron transport are found in Drift_Axial.py

### Test_Radial.py
This file will simulate the drifting of charge along the radial direction.
The key difference between this and the axial case is, in this case, the field is not uniform. 
There is a radial dependence that will effect the electron transport.
To deal with these complications, all of the kinematics are calculated within the program at various time steps. 
This is also done in the axial case (for consistentcy).
The functions for the electron transport are found in Drift_Radial.py

### LOR2CAD.py
Once the LORs are made, they are passed back to LOR2CAD.py
Here, they will be transformed into a set of sinograms. 
These sinograms will be converted to the image.
The reconstructed image will then be converted to an stl file with the help of to_stl.py

### to_stl.py
This file contains the functions that write an stl file
Each scatter plot point is converted to a box with side length equal to the pixel size in the reconstructed image 2D histogram.
These boxes are written out to an stl file. 


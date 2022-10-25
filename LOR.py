"""
Written by Tom Murphy
In Positron Emission Tomography:    
    A radionuclide decays via beta +
    The resulting positron annhilates with a neighboring electron
    This produces 2, 511 keV photons back to back
    We detect these photons with a ring shaped detector
    A "Line of Response" (LOR) is drawn connecting the dots
    
This file simulates that and makes the lines of response 
"""
import numpy as np
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
import Test_Axial
import Test_Radial
import Drift_Radial
import Drift_Axial
import Generate

#   This class draws a cylinder with a radius on the same scale as the axes 
class data_linewidth_plot():
    def __init__(self, x, y, z, **kwargs):
        self.ax = kwargs.pop("ax", plt.gca())
        self.fig = self.ax.get_figure()
        self.lw_data = kwargs.pop("linewidth", 1)
        self.lw = 1
        self.fig.canvas.draw()

        self.ppd = 72./self.fig.dpi
        self.trans = self.ax.transData.transform
        self.linehandle, = self.ax.plot([],[],**kwargs)
        if "label" in kwargs: kwargs.pop("label")
        self.line, = self.ax.plot(x, y, z, **kwargs)
        self.line.set_color(self.linehandle.get_color())
        self._resize()
        self.cid = self.fig.canvas.mpl_connect('draw_event', self._resize)

    def _resize(self, event=None):
        lw =  ((self.trans((1, self.lw_data))-self.trans((0, 0)))*self.ppd)[1]
        if lw != self.lw:
            self.line.set_linewidth(lw)
            self.lw = lw
            #self._redraw_later()

    def _redraw_later(self):
        self.timer = self.fig.canvas.new_timer(interval=1)
        self.timer.single_shot = True
        self.timer.add_callback(lambda : self.fig.canvas.draw_idle())
        self.timer.start()

#   Function that makes "N" LORs
def Generate_LOR(N):
    #   Defines data structure
    #   One coordinate for axial drift data 
    #   One coorindate for radial drift data
    LOR_xa_data = []
    LOR_xr_data = []
    LOR_ya_data = []
    LOR_yr_data = []
    LOR_za_data = []
    LOR_zr_data = []
    
    #   The way the charge readout works is it will trigger within 0.2 microseconds 
    #   So the uncertainty in the drift direction (z) is the drift velocity times the velocity
    #   Axial drifts have constant electric field
    #   So axial drifts have a constant speed (and therefore uncertainty)
    #   In the radial direction, the field varies so its not as simple
    delta_t = 0.2   
    v = 2.175
    delta_drift_axial = delta_t*v
    
    #   This finds the drift speed in the radial direction
    E = Drift_Radial.efield(Generate.R1)
    delta_drift_radial = Drift_Radial.mu(E)*E*delta_t/10/1000000
    
    #   Here we will plot the lines of response
    """
    #   Declares figure
    fig = plt.figure()
    ax = fig.add_subplot(1, 2, 1, projection='3d')
    ax2 = fig.add_subplot(1, 2, 2, projection='3d')

    #   Sets limits on plots 
    ax.axes.set_xlim3d(left=-100, right=100) 
    ax.axes.set_ylim3d(bottom=-100, top=100) 
    ax.axes.set_zlim3d(bottom=-100, top=100) 
    
    ax2.axes.set_xlim3d(left=-100, right=100) 
    ax2.axes.set_ylim3d(bottom=-100, top=100) 
    ax2.axes.set_zlim3d(bottom=-100, top=100) 
    
    #   Used to make visualizations in CAD geometry
    ax.set_facecolor((167/255, 223/255, 237/255))
    ax.set_axis_off()
    ax2.set_axis_off()
    """

    #   Iterates N times
    for i in range(N):
        
        #   Generates a random angle
        phi_axial = np.random.uniform(-np.pi,np.pi)
        phi_radial = np.random.uniform(-np.pi,np.pi)
        
        #   First interaction   #
        #   Preforms the Axial and radial drifts
        data_axial = Test_Axial.run(1,1,False,1500)
        data_radial = Test_Radial.run(1,1,False,Drift_Radial.N1(1500,E))
        
        #   Grabs the errors on the axial and radial drifts
        er_x1 = data_radial[2]
        er_y1 = data_radial[3]
        
        #   Generates a drift error
        er_z1 = np.random.uniform(0,delta_drift_radial)

        

        #   Second interaction  #
        data_axial = Test_Axial.run(1,1,False,1500)
        data_radial = Test_Radial.run(1,1,False,Drift_Radial.N1(1500,E)) 
        er_x2 = data_radial[2]
        er_y2 = data_radial[3]
        er_z2 = np.random.uniform(0,delta_drift_radial)
        
        #   Axial data for plotting
        """
        Ra = (((abs(er_x1)+abs(er_x2))/2)**2+((abs(er_y1)+abs(er_y2))/2)**2+((abs(er_z1)+abs(er_z2))/2)**2)**0.5/400
        x0a = (data_axial[4])*np.cos(phi_axial)*1000 
        xoff = np.random.uniform(-10,10)
        yoff = np.random.uniform(-10,10)
        zoff = np.random.uniform(-10,10)
        np.random.uniform(-5,5)
        xa = [-x0a+xoff+er_x1, x0a+xoff+er_x2]
        LOR_xa_data += [xa]
        y0a = data_axial[4]*np.sin(phi_axial)*1000
        ya = [-y0a+yoff+er_y1, y0a+yoff+er_y2]
        LOR_ya_data += [ya]
        z0a = data_axial[5]*1000
        za = [-z0a+zoff+er_z1, z0a+zoff+er_z2]
        LOR_za_data += [za]
        """
        
        #   Radial error data for plotting
        Rr = (((abs(er_x1)+abs(er_x2))/2)**2+((abs(er_y1)+abs(er_y2))/2)**2+((abs(er_z1)+abs(er_z2))/2)**2)**0.5/400
        
        #   NEMA nu 2 Standards
        #   Still need to do 
        theta = np.pi/3
        #R  = 2.5*3.7
        R = 20
        sphere_r = 5
        while theta < 2*np.pi:
            r_x = np.random.uniform(-sphere_r,sphere_r)
            r_y = np.random.uniform(-sphere_r,sphere_r)
            r_z = np.random.uniform(-sphere_r,sphere_r)
            r_test = np.sqrt(r_x**2+r_y**2+r_z**2)
            while r_test > sphere_r:
                r_x = np.random.uniform(-sphere_r,sphere_r)
                r_y = np.random.uniform(-sphere_r,sphere_r)
                r_z = np.random.uniform(-sphere_r,sphere_r)
                r_test = np.sqrt(r_x**2+r_y**2+r_z**2)
            xoff = R*np.cos(theta) + r_x
            yoff = R*np.sin(theta) + r_y
            zoff = 0 + r_z
            x0r = (data_radial[4])*np.cos(phi_radial)*1000 
            xr = [-x0r+xoff+er_x1, x0r+xoff+er_x2]
            LOR_xr_data += [xr]
            y0r = data_radial[4]*np.sin(phi_radial)*1000
            yr = [-y0r+yoff+er_y1, y0r+yoff+er_y2]
            LOR_yr_data += [yr]
            z0r = data_radial[5]*1000
            zr = [-z0r+zoff+er_z1, z0r+zoff+er_z2]
            LOR_zr_data += [zr]
            theta += np.pi/3
        """
        
        #   This section forms the objects that emit the radiation
        #   This code is broken up into 2 different sections
        #   The first of which creates a cube of side length 10 (R)
        #   To do this we make a random distribution of points between -R and R
        #   We can add some offset by adding a constant
        #   Here the object is moved 10 units on the x axis
                #   25 units on the y axis  
                #   50 units on the z axis 
        
        R = 10
        xoff =  np.random.uniform(-R,R)+10
        yoff =  np.random.uniform(-R,R)+25
        zoff =  np.random.uniform(-R,R)+50
        
        x0r = (data_radial[4])*np.cos(phi_radial)*1000 
        xr = [-x0r+xoff+er_x1, x0r+xoff+er_x2]
        LOR_xr_data += [xr]
        y0r = data_radial[4]*np.sin(phi_radial)*1000
        yr = [-y0r+yoff+er_y1, y0r+yoff+er_y2]
        LOR_yr_data += [yr]
        z0r = data_radial[5]*1000
        zr = [-z0r+zoff+er_z1, z0r+zoff+er_z2]
        LOR_zr_data += [zr]
        
        #   Second cube creation
        
        if i%2 == 0:
            xoff =  np.random.uniform(-R/2,R/2)-10
            yoff =  np.random.uniform(-R/2,R/2)-25
            zoff =  np.random.uniform(-R/2,R/2)-50

        x0r = (data_radial[4])*np.cos(phi_radial)*1000 
        xr = [-x0r+xoff+er_x1, x0r+xoff+er_x2]
        LOR_xr_data += [xr]
        y0r = data_radial[4]*np.sin(phi_radial)*1000
        yr = [-y0r+yoff+er_y1, y0r+yoff+er_y2]
        LOR_yr_data += [yr]
        z0r = data_radial[5]*1000
        zr = [-z0r+zoff+er_z1, z0r+zoff+er_z2]
        LOR_zr_data += [zr]
        """       
        
        #   Creates line from data
        """
        l = data_linewidth_plot( xa, ya, za, ax=ax, solid_capstyle="round", 
                            linewidth = 2*Ra, color = "black", alpha = 0.05)
              
        l = data_linewidth_plot( xr, yr, zr, ax=ax, solid_capstyle="round", 
                            linewidth = 2*Rr, color = "red", alpha = 0.05)
        
        l = data_linewidth_plot( xr, yr, zr, ax=ax2, solid_capstyle="round", 
                            linewidth = 2*Rr, color = "red", alpha = 0.05)
        """
    #   Shows plot
    plt.show()
    
    #   Returns data for a drift type
    return LOR_xr_data,LOR_yr_data,LOR_zr_data#, LOR_xa_data,LOR_ya_data,LOR_za_data
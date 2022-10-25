"""
Written by Tom Murphy
This code simulates axial drift of electrons in a cylindrical TPC 
"""

import matplotlib.pyplot as plt
import numpy as np
import random
import Generate

#   Useful constants 

sigmaC = 5.0638886*10**(-7)
epsilon0 = 8.854187817*10**(-12)
mu_ = 0.32
drift_length = Generate.R1/1000
edep_l = 0.0003
q_e = 1.602*10**(-19)
m_e = 9.109*10**(-31)

a0 = 551.6
a1 = 7158.3
a2 = 4440.43
a3 = 4.29
a4 = 43.63
a5 = 0.2053

b0 = 0.0075
b1 = 742.9
b2 = 3269.6
b3 = 31678.2


N = 1000
tmax = 0.00002
divisions = 1000

t = np.linspace(0,tmax,divisions + 1)
dt = tmax/divisions
r = np.zeros(divisions+1)

phi0 = random.uniform(0,2*np.pi)
rho0 = Generate.R0/1000

#   Electron mobility calculation
def mu(E):      # cm^2/Vs
    E = E/100000
    return (a0+a1*E+a2*E**(3/2)+a3*E**(5/2))/(1+(a1/a0)*E+a4*E**2+a5*E**3)*(87/89)**(-3/2)

#   Effective longitudinal electron energy
def eL(E):      # eV
    E = E/100000
    return (b0+b1*E+b2*E**2)/(1+(b1/b0)*E+b3*E**2)

#   Finds the transeverse diffusion coefficent given an Electric field
def Dt(E, mu, epsilon_L):       # cm^2/s
    E = E/100000
    dmudE = (((a1+3/2*a2*E**(1/2)+5/2*a3*E**(3/2))*(1+(a1/a0)*E+a4*E**2+a5*E**3))-((a1/a0+2*a4+3*a5)*(a0+a1*E+a2*E**(3/2)+a3*E**(5/2))))/((1+(a1/a0)*E+a4*E**2+a5*E**3)**2)
    DL = mu*epsilon_L*(87/89)**(-3/2)
    return DL/(1+E/mu*dmudE)

#   Finds the electric field
def efield():   # V/m
    return 100000

#   Finds the diffusion 
def diffusion(rho_t, e_id):     # mm
    return np.sqrt(2*Dt*t)

#Takes meters as units
#   Finds the standard deviation of the electron cloud for a given starting point
def stdev(z0):      # mm
    E_ = efield()
    mu_ = mu(E_)
    eL_ = eL(E_)
    sigma = np.sqrt((Generate.Z/1000-abs(z0))/(mu(E_)*E_)*Dt(E_,mu_,eL_)*2)*1000
    return sigma


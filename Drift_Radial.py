"""
Written by Tom Murphy
This code simulates radial drift of electrons in a cylindrical TPC 
"""

import matplotlib.pyplot as plt
import numpy as np
import random
import Generate

#   Useful constants 

sigmaC = 2.78162514*10**(-4) # delta V = 5,500 Volts Efield is 1kV/cm at R = 0.05 m
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
rho0 = Generate.R0

#   Electron mobility calculation
def mu(E):              # cm^2/Vs
    E = E/100000
    return (a0+a1*E+a2*E**(3/2)+a3*E**(5/2))/(1+(a1/a0)*E+a4*E**2+a5*E**3)*(87/89)**(-3/2)

#   Effective longitudinal electron energy
def eL(E):              # eV
    E = E/100000
    return (b0+b1*E+b2*E**2)/(1+(b1/b0)*E+b3*E**2)

#   Recombination Coefficent
def Reco(E):            #
    E = E/100000
    A = 0.8
    k = 0.0486
    dedx = 2.1
    rho = 1.396
    return A/(1+k/(E*rho)*dedx)

#   Converts from number of electrons (N0) produced by the 1kV/cm to whatever 
def N1(N0,E):
    E = E
    return int(N0 * Reco(E)/Reco(100000))

#   Finds the transeverse diffusion coefficent given an Electric field
def Dt(E, mu, epsilon_L):   # cm^2/s
    E = E/100000
    dmudE = a0*(a3*E**(3/2)+3*np.sqrt(E)*(a3*E+a2)/2+a1)/(E*(a0*E*(a5*E+a4)+a1)+a0) - a0*(E**(3/2)*(a3*E+a2)+a1*E+a0)*(E*(a0*(a5*E+a4)+a0*a5*E)+a0*E*(a5*E+a4)+a1)/(E*(a0*E*(a5*E+a4)+a1)+a0)**2
    DL = mu*epsilon_L
    return DL/(1+E/mu*dmudE)

#   Finds the electric field at a radius of r
def efield(r):  # V/m
    return sigmaC/(2*np.pi*r*epsilon0)

#   Calculates diffusion over a distance 
def diffusion(rho_t, e_id): 
    return np.sqrt(2*Dt*t)

#   Calcualtes the stdev of electron cloud spread from a given starting point. Takes meters as units
def stdev(r0): # mm
    N = 1000
    tmax = 0.0002
    divisions = 1000

    t = np.linspace(0,tmax,divisions + 1)
    dt = tmax/divisions
    r = np.zeros(divisions+1)

    phi0 = random.uniform(0,2*np.pi)
    rho0 = Generate.R0/1000
    sigmas = []
    for i in range(len(t)):
        if r[i-1]<drift_length:
            if i == 0:
                r[i] = r0/1000
            else:
                r[i] = r[i-1]+mu_*efield(r[i-1]*1000)*dt/10000
        else:
            r[i] = drift_length
        
        mu_ = mu(efield(r[i]*1000))
        eL_ = eL(efield(r[i]*1000))
        Dt_ = Dt(efield(r[i]*1000), mu_, eL_)
        if r[i-1]<drift_length:
            sigmas += [np.sqrt(2*Dt_*dt)]

    sigma = 1/np.sqrt(len(sigmas))*sum(sigmas) * 10
    return sigma
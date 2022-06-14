#Mike Lafky
#USPAS June 2022
#HW2.6

from cmath import sqrt
import numpy as np
import matplotlib.pyplot as plt

#constants
c = 299792458               #c, m/s
h_bar = 6.58E-16            #Planck constant, ev*s
Ebeam = 1.5                 #Energy of e-beam, GeV
m0 = .511                   #Mass of e-, MeV
m0_gev = m0/1000            #Mass of e-, GeV
gamma = (Ebeam*1000+m0)/m0  #Lorentz factor
L = 5                       #Undulator length, m
e = 1.6E-19                 #electron charge, coulombs
alpha = 1/137               #fine structure constant, unitless

#free parameters
#ONLY CHANGE THESE
#gap < lam_u
lam_u = .05                 #undulator period length, m
gap = .015                  #gap, m


#Given undulator parameters
Nu = 100            #Undulator period, unitless
lam_u = .05         #Undulator wavelength, m
k = 1               #Deflection parameter, unitless

#Given beam parameters
emit_x = 4E-9       #Horizontal emittance, m-rad
emit_y = 4E-11      #Vertical emittance, m-rad
beta_x = 5          #Horizontal beta (twiss parameter), m
beta_y = 5          #Vertical beta (twiss parameter), m
curr = 100          #e- beam current, A

#Calculate beam size from given parameters
sigx = sqrt(emit_x*beta_x)
sigx = sqrt(emit_y*beta_y)
#establish a logspace for indpendent variable, photon energy
Eph = np.logspace(10,10000,1000)    #energy from 10eV to 10keV, log scale
flux = np.zeros(len(Eph))           #array for the flux calculations

#re-write equation 2.75 as a function of photon energy by combining eqs 2.52 & 2.53, and solving those for gamma^2 and inserting into 2.75
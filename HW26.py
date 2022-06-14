#Mike Lafky
#USPAS June 2022
#HW2.6

from math import sqrt, pi
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy import special
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
Eph = np.logspace(1,4,1000)    #energy from 10eV to 10keV, log scale
flux = np.zeros(len(Eph))           #array for the flux, photons/(s*mrad*.1%BW)
brightness = np.zeros(len(Eph))     #array for brightness, photons/(s*mm*2*mrad^2*.1%BW)
m_rel = np.logspace(0,1,1000)      #Just a test, e- energy in GeV

#emittance-r = lambda (photons/4*pi, 1.59)
#50 eV photon lambda = 2.48E-8 m
#10 keV photon lambda = 1.24E-10 m
#lam_r = lam_wavelength/4*math.pi
#apparently I should be able to solve 2.109 and incorporate 2.106
#Eph=freq*h_bar=h_bar*


#re-write equation 2.75 as a function of photon energy by combining eqs 2.52 & 2.53, and solving those for gamma^2 and inserting into 2.75
#Calculated fundamental parameters
for x in range(len(Eph)):
    x_1 = (1*k**2)/(4+2*(k**2))                                 #input for bessell function
    jj_1 = ((-1)**0)*(special.jv(0,(x_1))-special.jv(1,(x_1)))  #bessell function
    flux[x]=(1.74E14*(Nu**2)*(m_rel[x]**2)*curr*(k**2)*(jj_1**2))/(1/(1+(k**2)/2))
    #flux[x] = (1.74E14*Nu**2*((lam_u*Eph[x])/(h_bar*4*pi*c)))*(((m0_gev**2)*(k**2)*(jj_1**2)*curr)/(1+k**2/2))
    #brightness[x]=flux[x]/(emit_x*emit_y)
    #print("Eph = " + str(round(Eph[x],3)))
    flux[x]=(math.pi*alpha*Nu*(curr/e)*.1*k**2*jj_1**2)/(1/(1+(k**2)/2))*((h_bar*2*math.pi*c)/(Eph[x]**2))

#plt.plot(Eph,brightness)
print(Eph)
plt.plot(Eph,flux)
plt.yscale('log')
plt.xscale('log')
plt.xlabel('Photon energy (eV)')
plt.ylabel('Spectral Flux (photons/s/mm*2/mrad*2 .1 %BW')
plt.show()
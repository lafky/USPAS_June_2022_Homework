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
e = 1.6E-19                 #electron charge, coulombs
alpha = 1/137               #fine structure constant, unitless

#Given undulator parameters
Nu = 100            #Undulator period, unitless
lam_u = .05         #Undulator wavelength, m
l = 5               #Undulator length, m
k = 1               #Deflection parameter, unitless

#Given beam parameters
emit_x = 4E-9       #Horizontal emittance, m-rad
emit_y = 4E-11      #Vertical emittance, m-rad
beta_x = 5          #Horizontal beta (twiss parameter), m
beta_y = 5          #Vertical beta (twiss parameter), m
curr = .100          #e- beam current, A

#Calculate beam size from given parameters
sigx = sqrt(emit_x*beta_x)  #horizontal e- beam size, m
sigxp = emit_x/sigx         #horizontal beam angular divergnce, rad
sigy = sqrt(emit_y*beta_y)  #vertical e- beam size, m
sigyp = emit_y/sigy         #vertical beam angular divergnce, rad

#establish a logspace for indpendent variable, photon energy
Eph = np.logspace(1,4,1000)         #energy from 10eV to 10keV, log scale
brightness = np.zeros(len(Eph))     #array for brightness, photons/(s*mm*2*mrad^2*.1%BW)
coh = np.zeros(len(Eph))            #array for coherence
x_1 = (1*k**2)/(4+2*(k**2))                                 #input for bessell function
jj_1 = ((-1)**0)*(special.jv(0,(x_1))-special.jv(1,(x_1)))  #bessell function

for x in range(len(Eph)):
    lam_e = (h_bar*2*pi*c)/Eph[x]                               #photon wavelength at this energy, m
    sigr = sqrt(2*lam_e*l)/(4*pi)                               #photon beam size, m
    sigrp = sqrt(lam_e/(2*l))                                   #photon beam divergence at this energy, rad
    bigsigx = sqrt(sigx**2+sigr**2)                             #Inputs for brilliance calc
    bigsigxp = sqrt(sigxp**2+sigrp**2)
    bigsigy = sqrt(sigy**2+sigr**2)
    bigsigyp = sqrt(sigyp**2+sigrp**2)
    brightness[x]=(pi*alpha*Nu*(curr/e)*(.1/100)*k**2*jj_1**2)/((4*pi**2)*(bigsigx*bigsigxp*bigsigy*bigsigyp)*(1/(1+(k**2)/2)))
    coh[x] = (bigsigx*bigsigy*bigsigy*bigsigyp-((lam_e/(4*pi)))**2)/((lam_e/(4*pi))**2)
    #print(bigsigx,bigsigxp,bigsigy,bigsigyp)
    #print(brightness[x]/1E12)
    #print(bigsigx*bigsigxp, bigsigy*bigsigyp,curr/e, alpha*Nu)

#Convert to mm-mrad
brightness=brightness/1E12
plt.plot(Eph,brightness)
plt.yscale('log')
plt.xscale('log')
plt.xlabel('Photon energy (eV)')
plt.ylabel('Spectral Flux (photons/s/mm*2/mrad*2 .1 %BW')
plt.show()
plt.plot(Eph,coh)
plt.yscale('log')
plt.xscale('log')
plt.ylabel('Transverse Coherence (Sigx*Sigx\'*Sigy*Sigy\'/(lam_1/4*pi)^2)')
plt.xlabel('Photon energy (eV)')
plt.show()

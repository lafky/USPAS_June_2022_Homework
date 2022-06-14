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
sigx = sqrt(emit_x*beta_x)*1000  #horizontal e- beam size, mm
sigxp = emit_x/sigx*1000         #horizontal beam angular divergnce, mrad
sigy = sqrt(emit_y*beta_y)*1000  #vertical e- beam size, mm
sigyp = emit_y/sigy*1000         #vertical beam angular divergnce, mrad

#establish a logspace for indpendent variable, photon energy
Eph = np.logspace(1,4,1000)         #energy from 10eV to 10keV, log scale
brightness = np.zeros(len(Eph))     #array for brightness, photons/(s*mm*2*mrad^2*.1%BW)

for x in range(len(Eph)):
    x_1 = (1*k**2)/(4+2*(k**2))                                 #input for bessell function
    jj_1 = ((-1)**0)*(special.jv(0,(x_1))-special.jv(1,(x_1)))  #bessell function
    lam_e = (h_bar*2*pi*c)/Eph[x]                               #photon wavelength at this energy, m
    sigr = sqrt(2*lam_e*l)/(4*pi)*1000                          #photon beam size, mm
    sigrp = sqrt(lam_e/(2*l))*1000                              #photon beam divergence at this energy, mrad
    bigsigx = sqrt(sigx**2+sigr**2)                             #Inputs for brilliance calc
    bigsigxp = sqrt(sigxp**2+sigrp**2)
    bigsigy = sqrt(sigy**2+sigr**2)
    bigsigyp = sqrt(sigyp**2+sigrp**2)
    brightness[x]=(pi*alpha*Nu*(curr/e)*.1*k**2*jj_1**2)/(4*pi**2*bigsigx*bigsigxp*bigsigy*bigsigyp*(1/(1+(k**2)/2)))
    #print(bigsigx*bigsigxp, bigsigy*bigsigyp,curr/e, alpha*Nu)

plt.plot(Eph,brightness)
plt.yscale('log')
plt.xscale('log')
plt.xlabel('Photon energy (eV)')
plt.ylabel('Spectral Flux (photons/s/mm*2/mrad*2 .1 %BW')
plt.show()
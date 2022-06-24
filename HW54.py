#Mike Lafky
#USPAS June 2022
#HW5.4

from math import sqrt, pi
import math
from selectors import EVENT_READ
import numpy as np
import matplotlib.pyplot as plt
from scipy import special

#constants
c = 299792458               #c, m/s
h_bar_ev = 6.58E-16         #Planck constant, ev*s
h_bar_j = 1.054571817E34    #Planck constant, J⋅s
m0 = .511                   #Mass of e-, MeV
m0_gev = m0/1000            #Mass of e-, GeV
e = 1.6E-19                 #electron charge, C
alpha = 1/137               #fine structure constant, unitless
IA = 17045                  #Alfven constant, A

#Fel parameters
eBeam = 13.6                   #Energy of e-beam, GeV
eBeam_j = eBeam*1.60218E-10     #Energy of e-beam, J
gamma = (eBeam*1000+m0)/m0      #Lorentz factor
emit_norm =  4E-7               #normalized emittance, m rad
emit = 4E-7/gamma               #unnormalized emittance, m rad
curr = 3000                     #beam current, A
lam_u = .03                     #period length, m
beta = 140/(2*np.pi)            #average beta function, m
sig_g = 0.0001                  #energy spread 
lambda_1 = 1.5E-10              #fundamental wavelength, m
k = 3.5                         #FEL focusing parameter, unitless

#function parameters
sig_x = sqrt(beta*emit)         #beam size, m
sig_xp = sqrt(emit/sig_x)       #beam divergence, rad
kb = 1/beta                     #average focusing parameter, 1/m
a = np.array([.45, .57, .55, 1.6, 3, 2, .35, 2.9, 2.4, 51, .95, 3, 5.4, .7, 1.9, 1140, 2.2, 2.9, 3.2])
eta_d = .5*(1/sqrt(3))*(sig_x**2)       #diffraction parameter
eta_e = (2/sqrt(3))*(kb**2)*(sig_x**2)  #angular spread parameter
eta_g = sig_g/sqrt(3)                   #energy spread parameter

#Calculations
x_1 = (1*k**2)/(4+2*(k**2))                                 #input for bessell function
jj_1 = ((-1)**0)*(special.jv(0,(x_1))-special.jv(1,(x_1)))  #bessell function
#ro = ((1/(8*pi))*(curr/IA)*(((k*jj_1)/(1+((k**2)/2)))**2)*((gamma*(lambda_1**2))/(2*pi*sig_x**2)))**(1/3) #pierce parameter

term1 = ((1/(8*pi))*(curr/IA))
term2 = (((k*jj_1)/(1+((k**2)/2)))**2)
term3 = ((gamma*lambda_1**2)/(2*pi*(sig_x**2)))
ro = (term1*term2*term3)**(1/3)
print(term1, term2, term3, ro)

big_lamb_1 = a[0]*eta_d**a[1]+a[2]*eta_e**a[3]+a[4]*eta_g**a[6]+a[6]*(eta_e**a[7])*(eta_g**a[8])
big_lamb_2 = a[9]*(eta_d**a[10])*(eta_g**a[12])+a[12]*(eta_d**a[13])*(eta_e**a[14])
big_lamb_3 = a[15]*(eta_d**a[16])*(eta_e**a[17])*(eta_g**a[18])
big_lamb = big_lamb_1 + big_lamb_2 + big_lamb_3 #parameter that defines gain length 

#print(ro)

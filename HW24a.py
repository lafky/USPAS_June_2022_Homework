#Mike Lafky
#USPAS June 2022
#HW2.4a

from cmath import nan, pi
import numpy as np
import matplotlib.pyplot as plt
import cmath
import math

#constants
c = 299792458   #c, m/s
h_bar = 6.58E-16    #Planck constant, ev*s
Ebeam = 1.4     #Energy of e-beam, GeV
m0 = .511       #Mass of e-, MeV
gamma = (Ebeam*1000+m0)/m0 #Lorentz factor
g_min = .015    #gap minimum, m

#free parameters
#gap must be > lam_u
lam_u = .02     #undulator period length, m
gap = .03       #gap, m

#Calculated values
B0 = 3.44*math.exp(-(gap/lam_u)*(5.08-1.54*(gap/lam_u)))  #Field strength, T
k = .934*lam_u*B0 #deflection parameter, m*T

lam_1 = (lam_u/(2*gamma*gamma))*(1+k*k/2)
omg_1 = (2*pi*c)/lam_1
E_1 = omg_1*h_bar

print('e- Energy: ' + str(Ebeam) + ' GeV')
print('Lorentz Factor: ' + str(round(gamma)))
print('For undulator period length ' + str(lam_u*100) + ' (cm) and gap ' + str(round(gap*100,3)) + ' (cm) ')
print('B0 = ' + str(round(B0,3)) + ' T and K = ' + str(round(k,4)) + ' (T * cm)')
print('λ1 = ' + str(round(lam_1,12)) + ' (m)' + ' ω1 = ' + str(round(omg_1,12)))
print('E1 = ' + str(round(E_1,1)) + ' eV')
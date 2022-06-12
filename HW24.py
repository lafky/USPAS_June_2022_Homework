#Mike Lafky
#USPAS June 2022
#HW2.4a

from cmath import pi
import math
from scipy import special

#constants
c = 299792458               #c, m/s
h_bar = 6.58E-16            #Planck constant, ev*s
Ebeam = 1.4                 #Energy of e-beam, GeV
m0 = .511                   #Mass of e-, MeV
m0_gev = m0/1000            #Mass of e-, GeV
curr = .4                   #beam current, Amps
gamma = (Ebeam*1000+m0)/m0  #Lorentz factor
L = 5                       #Undulator length, m

#free parameters
#ONLY CHANGE THESE
#gap < lam_u
lam_u = .05                 #undulator period length, m
gap = .015                  #gap, m

#Calculated undulator parameters
B0 = 3.44*math.exp(-(gap/lam_u)*(5.08-1.54*(gap/lam_u)))    #Field strength, T
k = 93.4*lam_u*B0                                           #deflection parameter, unitless
Nu = int(L/lam_u)                                           #number of periods


#Calculated fundamental parameters
lam_1 = (lam_u/(2*gamma*gamma))*(1+k*k/2)                   #fundamental wave length, m
omg_1 = (2*pi*c)/lam_1                                      #fundamental frequency, s^-1
E_1 = omg_1*h_bar                                           #fundamental energy, eV
sig_1 = (lam_1/(2*L))**.5                                   #RMS Angular Divergence
x_1 = (1*k**2)/(4+2*(k**2))                                 #input for bessell function
jj_1 = ((-1)**0)*(special.jv(0,(x_1))-special.jv(1,(x_1)))  #bessell function
flux_1 = (1.74E14)*(Nu**2)*((gamma*m0_gev)**2)*curr*(1**2)*(k**2)*(jj_1)**2/(1+(k**2)/2)    #flux, photons/(s*mrad*.1% BW)

#Calculated third harmonic parameters
gam_r = gamma                                               #resonance gamma of the third harmonic, unitless
lam_3 = ((1+k**2/2)/(2*gam_r**2))*(lam_u/3)                 #3rd harmonic wavelength, m
omg_3 = omg_1 = (2*pi*c)/lam_3                              #3rd harmonic frequency, s^-1
E_3 = omg_3*h_bar                                           #3rd harmonic energy, eV
sig_3 = (lam_3/(2*L))**.5                                   #RMS Angular Divergence
x_3 = (1*k**2)/(4+2*(k**2))                                 #input for bessell function
jj_3 = ((-1)**1)*(special.jv(1,(x_3))-special.jv(2,(x_3)))  #bessell function
flux_3 = (1.74E14)*(Nu**2)*((gamma*m0_gev)**2)*curr*(3**2)*(k**2)*(jj_3)**2/(1+(k**2)/2)    #flux, photons/(s*mrad*.1% BW)

print('e- Energy: ' + str(Ebeam) + ' GeV')
print('Lorentz Factor: ' + str(round(gamma)))
print('For undulator period length ' + str(lam_u*100) + ' (cm), gap ' + str(round(gap*100,3)) + ' (cm), ' + 'Nu = ' + str(Nu))
print('B0 = ' + str(format(B0,'.3E')) + ' T and K = ' + str(format(k,'.3E')))
print('FUNDAMENTAL:')
print('λ1 = ' + str(format(lam_1,'.3E')) + ' (m)' + ' ω1 = ' + str(format(omg_1,'.3E')) + ' (s^-1)')
print('E1 = ' + str(round(E_1,1)) + ' eV.  RMS Angular Divergence = ' + str(format(sig_1,'.3E')) + ', flux = ' + str(format(flux_1,'.3E')) + ' photons/(s*mrad*.1% BW)')
print('THIRD HARMONIC:')
print('λ3 = ' + str(format(lam_3,'.3E')) + ' (m)' + ' ω3 = ' + str(format(omg_3,'.3E')) + ' (s^-1)')
print('E3 = ' + str(round(E_3,1)) + ' eV.  RMS Angular Divergence = ' + str(format(sig_3,'.3E')) + ', flux = ' + str(format(flux_3,'.3E')) + ' photons/(s*mrad*.1% BW)')
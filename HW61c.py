#Mike Lafky
#USPAS June 2022
#HW6.1c

from math import sqrt, pi
import cmath
import numpy as np
import matplotlib.pyplot as plt
from scipy import special

#constants
c = 299792458               #c, m/s
h_bar_ev = 6.58E-16         #Planck constant, ev*s
h_bar_j = 1.054571817E34    #Planck constant, J⋅s
m0 = .511                   #Mass of e-, MeV
m0_gev = m0/1000            #Mass of e-, GeV
mo_kg = 9.1093837014E-31    #Mass of 2-, kg
e = 2.718281828459          
alpha = 1/137               #fine structure constant, unitless
IA = 17045                  #Alfven constant, A

compl = 1j
mu = (-1+compl*sqrt(3))/2

#k1 = 3.7                       #LCLS Focusing parameter, unitless
z = 70                          #LCLS undulator length, m
lambda_u = .03                  #LCLS undulator wavelength, m
ku = 2*pi/lambda_u              #LCLS Ku
rho = 6.23E-4                   #LCLS Pierce Parameter
B = 2*rho*ku*z*mu*compl         #Exponent

k = np.logspace(.1,10,100)  #Array of K values
h = np.array([3])  #Array of harmonics
p_h = np.zeros((len(h),len(k)))      #Power calc array
rho = np.zeros((len(h),len(k)))      #bandwidth array

for x in range(len(h)):
    for y in range(len(k)):
        
        x_1 = (1*(k[y]**2))/(4+2*(k[y]**2))    #input for fundamental bessel function function
        jj_1 = ((-1)**1)*(special.jv(0,(x_1))-special.jv(1,(x_1)))  #fundamental bessel function function
        
        x_h = h[x]*(h[x]*(1+k[y]**2)-2) / (4+2*(h[x]*(1+k[y]**2)-2))                                           #input for harmonic bessell function function
        jj_h = ((-1)**h[x])*(special.jv((h[x]-1)/2,(x_h))-special.jv((h[x]+1)/2,(x_h)))  #bessell function function at each harmonic

        rho_p_h = ((sqrt(h[x]*(1+(k[y]**2))-2))*jj_h/ (1+((h[x]*(1+k[y]**2)-2)**2)/2))**2 #pierce parameter at the harmonic
        rho_f_h = ((((k[y]*jj_1)/(1+(k[y]**2)/2))**(2)))**2  #pierce parameter at the fundamental
        rho[x][y] = (rho_p_h/rho_f_h)**(1/3)
        #print(jj_h, rho_p_h,rho_f_h, rho[x][y])


#print(rho)
fig, ax1 = plt.subplots()
ax1.plot(k, rho[0], 'b', label = '3rd')
ax1.set_title('Homework 6.1 a ii')
ax1.legend()
ax1.set_yscale('linear')
ax1.set_xscale('linear')
ax1.set_xlim([0,10])
ax1.set_xlabel('k')
ax1.set_ylabel('Ratio of rho(3) to rho(1)')

plt.show()

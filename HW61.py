#Mike Lafky
#USPAS June 2022
#HW6.1

from math import sqrt, pi
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


k = 1                       #Arbitrarily chosen undulator k parameter
A = .1                      #Arbitrarily chosen 6.12 coefficient
B = 1                       #Arbitrarily chosen for z*mu*i in 6.12 exponential

z = 70                      #LCLS undulator length, m
lambda_u = .03              #LCLS undulator wavelength, m
ku = 2*pi/lambda_u
rho = 6.23E-4               #LCLS Pierce Parameter
B = -rho*ku*z               #Closer to reality B
print(B)

h = np.array([3, 5, 7, 9])  #Array of harmonics
p_h = np.zeros(len(h))     #Power calc array

x_1 = (1*(k**2))/(4+2*(k**2))    #input for fundamental bessel function function
jj_1 = ((-1)**0)*(special.jv(0,(x_1))-special.jv(1,(x_1)))  #fundamental bessel function function

x_h = (h*(k**2))/(4+2*(k**2))                                 #input for harmonic bessell function function
jj_h = ((-1)**h)*(special.jv((h-1)/2,(x_h))-special.jv((h+1)/2,(x_h)))  #bessell function function at each harmonic

for x in range(len(h)):
    p_h[x] = h[x]* ((jj_h[x]**2)/(jj_1**2)) * A**(h[x]-1) *e**(B*(h[x]+1))

#xax = []
#xax = h
fig, ax1 = plt.subplots()
ax1.bar(h, p_h)
ax1.set_yscale('log')
ax1.set_xlabel('harmonic #)')
ax1.set_ylabel('% Power of fundamental')
plt.show()
#print(x_h,jj_h, jj_1, p_h)
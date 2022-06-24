#Mike Lafky
#USPAS June 2022
#HW5.4

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
e = 1.6E-19                 #electron charge, C
alpha = 1/137               #fine structure constant, unitless
IA = 17045                  #Alfven constant, A
a = np.array([.45, .57, .55, 1.6, 3, 2, .35, 2.9, 2.4, 51, .95, 3, 5.4, .7, 1.9, 1140, 2.2, 2.9, 3.2])  #fitting parameter array

#Fel parameters
eBeam = 13.6                    #Energy of e-beam, GeV
eBeam_j = eBeam*1.60218E-10     #Energy of e-beam, J
gamma = (eBeam*1000+m0)/m0      #Lorentz factor
emit_norm =  4E-7               #normalized emittance, m rad
emit = 4E-7/gamma               #unnormalized emittance, m rad
curr = 3000                     #beam current, A
lam_u = .03                     #period length, m
sig_g = 0.0001                  #energy spread 
lambda_1 = 1.5E-10              #fundamental wavelength, m
k = 3.5                         #FEL focusing parameter, unitless
#pBeam = eBeam_j*(curr/e)              #beam power, W
pBeam = gamma*mo_kg*(c**2)*(curr/e)    #beam power, W
#pBeam  = eBeam*(curr/(e*1000))        #beam power in book units

#parameters for plotting
beta = np.linspace(1,25,1000)
p = np.zeros(len(beta))
LG = np.zeros(len(beta))

#loop to calculate P at each z
for x in range(len(beta)):   

    #Power Calculations
    sig_x = sqrt(beta[x]*emit)              #beam size, m
    sig_xp = sqrt(emit/sig_x)               #beam divergence, rad    
    kb = 1/beta[x]                          #average focusing parameter, 1/m

    #Pierce parameter Calculations
    x_1 = (1*k**2)/(4+2*(k**2))                                 #input for bessell function
    jj_1 = ((-1)**0)*(special.jv(0,(x_1))-special.jv(1,(x_1)))  #bessell function
    ro_term1 = ((1/(8*pi))*(curr/IA))
    ro_term2 = (((k*jj_1)/(1+((k**2)/2)))**2)
    ro_term3 = ((gamma*lambda_1**2)/(2*pi*(sig_x**2)))
    ro = (ro_term1*ro_term2*ro_term3)**(1/3)                    #Pierce parameter

    sig_d = .5*(1/sqrt(3))*(sig_x**2)                           #diffraction parameter
    sig_e = (2/sqrt(3))*(kb**2)*(sig_x**2)                      #angular spread parameter
    sig_g = sig_g/sqrt(3)                                       #energy spread parameter

    big_lamb_1 = a[0]*sig_d**a[1]+a[2]*sig_e**a[3]+a[4]*sig_g**a[6]+a[6]*(sig_e**a[7])*(sig_g**a[8])
    big_lamb_2 = a[9]*(sig_d**a[10])*(sig_g**a[12])+a[12]*(sig_d**a[13])*(sig_e**a[14])
    big_lamb_3 = a[15]*(sig_d**a[16])*(sig_e**a[17])*(sig_g**a[18])
    big_lamb = big_lamb_1 + big_lamb_2 + big_lamb_3 #parameter that determines power growth 
    
    LG0 = lam_u/(4*pi*sqrt(3)*ro)                   #1-D gain length, m
    LG[x]  = LG0*(1+big_lamb)                       #gain length, m
    #p[x]=1.6*((LG0/LG[x])**2)*ro*pBeam             #Saturation Power, alternate calc
    p[x] = (1.6/((1+big_lamb)**2))*ro*pBeam         #Saturation Power, units depend on input power units


print(ro, LG0, pBeam)

#plot
fig, ax1 = plt.subplots()
ax2 = ax1.twinx()
ax1.plot(beta, p, 'r-')
ax2.plot(beta, LG, 'b-')
ax1.set_xlabel('Avg. Beta (m)')
ax1.set_ylabel('Power (w)', color='r')
ax2.set_ylabel('LG (m)', color='b')
plt.show()
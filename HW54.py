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

#Input FEL Parameters 
eBeam = 13.6                    #Energy of e-beam, GeV
eBeam_j = eBeam*1.60218E-10     #Energy of e-beam, J
gamma = (eBeam*1000+m0)/m0      #Lorentz factor
curr = 3000                     #Beam current, A
lam_u = .03                     #Period length, m
sig_e  = 0.0001                 #Energy spread 
lambda_1 = 1.5E-10              #Fundamental wavelength, m
k = 3.5                         #FEL focusing parameter, unitless
emit_norm =  4E-7               #Normalized emittance, m rad

#Calculated FEL Parameters
emit = 4E-7/gamma                       #unnormalized emittance, m rad
#pBeam = eBeam_j*(curr/e)               #beam power, W
pBeam = gamma*mo_kg*(c**2)*(curr/e)     #beam power, W
#pBeam  = eBeam*(curr/(e*1000))         #beam power in book units
x_1 = (1*k**2)/(4+2*(k**2))                                 #input for bessell function
jj_1 = ((-1)**0)*(special.jv(0,(x_1))-special.jv(1,(x_1)))  #bessell function

#parameters for plotting
beta = np.linspace(1,100,1000)
p = np.zeros(len(beta))
LG = np.zeros(len(beta))

#loop to calculate P at each beta(z)
for x in range(len(beta)):   

    #Power Calculations
    sig_x = sqrt(beta[x]*emit)              #beam size, m
    kb = 1/beta[x]                          #average focusing parameter, 1/m

    #Pierce parameter Calculations
    ro_term1 = ((1/(8*pi))*(curr/IA))
    ro_term2 = (((k*jj_1)/(1+((k**2)/2)))**2)
    ro_term3 = ((gamma*lambda_1**2)/(2*pi*(sig_x**2)))
    ro = (ro_term1*ro_term2*ro_term3)**(1/3)                    #Pierce parameter

    sig_x_hat = 2*pi*sig_x*sqrt((2*ro)/(lambda_1*lam_u))        #eq 5.124
    sig_eta_hat = sig_e/ro                                      #eq 5.121
    k_beta_hat = (emit*2*pi)/(lambda_1*(sig_x_hat**2))          #eq 5.122

    eta_d = (1/(2*sqrt(3)*(sig_x_hat**2)))                      #eq 5.158, diffraction parameter
    eta_e = (k_beta_hat**2)*(sig_x_hat**2)*(2/sqrt(3))          #eq 5.159, angular spread parameter
    eta_g = sig_eta_hat/sqrt(3)                                 #eq 5.160, energy spread parameter

    big_lamb_1 = a[0]*eta_d**a[1] + a[2]*eta_e**a[3] + a[4]*eta_g**a[5] + a[6]*(eta_e**a[7])*(eta_g**a[8])
    big_lamb_2 = a[9]*(eta_d**a[10])*(eta_g**a[11]) + a[12]*(eta_d**a[13])*(eta_e**a[14])
    big_lamb_3 = a[15]*(eta_d**a[16])*(eta_e**a[17])*(eta_g**a[18])
    big_lamb = big_lamb_1 + big_lamb_2 + big_lamb_3 #Lambda, parameter that determines power growth 
    
    LG0 = lam_u/(4*pi*sqrt(3)*ro)                   #1-D gain length, m

    #Actual plot variables:
    LG[x]  = LG0*(1+big_lamb)                       #gain length, m
    #p[x]=1.6*((LG0/LG[x])**2)*ro*pBeam             #Saturation Power, alternate calc
    p[x] = (1.6/((1+big_lamb)**2))*ro*pBeam         #Saturation Power, units depend on input power units

#plot
fig, ax1 = plt.subplots()
ax2 = ax1.twinx()
ax1.plot(beta, p, 'r-')
ax2.plot(beta, LG, 'b-')
ax1.set_xlabel('Avg. Beta (m)')
ax1.set_ylabel('Power (w)', color='r')
ax2.set_ylabel('LG (m)', color='b')
plt.show()
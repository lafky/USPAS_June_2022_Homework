#Mike Lafky
#USPAS June 2022
#HW2.3

from cmath import nan, pi
import numpy as np
import matplotlib.pyplot as plt
import cmath
import math
from scipy.fftpack import fft, rfft, fftfreq, rfftfreq

#Parameters
l = 5               #undulator length, m
c = 299792458       #c, m/s
Nu = 50             #Number of undulator periods, unitless
lamdaFund = l/Nu    #fundamental wavelength, I think, m
#omegaFund = (2*pi*c)/Nu    #Fundamental frequency, 1/s
omegaFund = 1       #Fundamental frequency, 1/s
h1 = 1              #1st Harmonic, unitless
h3 = 3              #3rd harmonic, unitless
freq = np.arange(.9*omegaFund,1.1*omegaFund,omegaFund/1000)/omegaFund  #independent variable
print(len(freq))
s1 = np.zeros((len(freq))) #dependent variable, h = 1
s3 = np.zeros((len(freq))) #dependent variable, h = 3
#print(s)
for x in range(len(freq)):
    s1[x] = abs(math.sin(pi*Nu*(freq[x]-h1*omegaFund)/omegaFund)/(pi*Nu*(freq[x]-h1*omegaFund)/omegaFund))**2
    s3[x] = abs(math.sin(pi*Nu*(freq[x]-h3*omegaFund)/omegaFund)/(pi*Nu*(freq[x]-h3*omegaFund)/omegaFund))**2
fig1, axs = plt.subplots(1,2,figsize=(16,9))
fig1.suptitle('')
axs[0].plot(freq,s1)
axs[0].set_title('Fundamental')
axs[0].set_xlabel('ω/ω1')
axs[0].set_ylabel('S(ω,phi')
axs[1].plot(freq,s3)
axs[1].set_xlabel('ω/ω1')
axs[1].set_ylabel('S(ω,phi)')
axs[1].set_title('3rd Harmonic')
plt.show()
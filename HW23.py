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
Nu = 100             #Number of undulator periods, unitless
lamdaFund = l/Nu    #fundamental wavelength, I think, m
numPoints = 10000   #number of points to plot
#omegaFund = (2*pi*c)/Nu    #Fundamental frequency, 1/s
omegaFund = 1       #Fundamental frequency, 1/s
h1 = 1              #1st Harmonic, unitless
h3 = 3              #3rd harmonic, unitless
h5 = 5
freq1 = np.arange(.9*omegaFund,1.1*omegaFund,omegaFund/numPoints)/omegaFund  #independent variable
freq3 = np.arange(2.9*omegaFund,3.1*omegaFund,omegaFund/numPoints)/omegaFund  #independent variable
freq5 = np.arange(4.9*omegaFund,5.1*omegaFund,omegaFund/numPoints)/omegaFund  #independent variable
s1 = np.zeros((len(freq1))) #dependent variable, h = 1
s3 = np.zeros((len(freq3))) #dependent variable, h = 3
s5 = np.zeros((len(freq5))) #dependent variable, h = 5

#Calculate "S" at each frequency for each harmonic
for x in range(len(freq1)-1):

    s1[x] = abs(math.sin(pi*Nu*(freq1[x]-h1*omegaFund)/omegaFund)/(pi*Nu*(freq1[x]-h1*omegaFund)/omegaFund))**2
    s3[x] = abs(math.sin(pi*Nu*(freq3[x]-h3*omegaFund)/omegaFund)/(pi*Nu*(freq3[x]-h3*omegaFund)/omegaFund))**2
    s5[x] = abs(math.sin(pi*Nu*(freq5[x]-h5*omegaFund)/omegaFund)/(pi*Nu*(freq5[x]-h5*omegaFund)/omegaFund))**2

#Calculate FWHM from given formula
FWHMCalcH1 = np.round(.9/(h1*Nu),4)
FWHMCalcH3 = np.round(.9/(h3*Nu),4)
FWHMCalcH5 = np.round(.9/(h5*Nu),4)

#Measure it from plot
s3_max=np.argmax(s3) #return max index from array
print(s3_max)



fig1, axs = plt.subplots(1,3,figsize=(16,9))
fig1.suptitle('')
axs[0].plot(freq1,s1)
axs[0].set_title('Fundamental')
axs[0].set_xlabel('ω/ω1')
axs[0].set_ylabel('S(ω,phi)')
axs[0].text(1,.1,'Calc FWHM: ' + str(FWHMCalcH1))
axs[1].plot(freq3,s3)
axs[1].set_xlabel('ω/ω1')
axs[1].set_ylabel('S(ω,phi)')
axs[1].set_title('3rd Harmonic')
axs[1].text(3,.1,'Calc FWHM: ' + str(FWHMCalcH3))
axs[2].plot(freq5,s5)
axs[2].set_xlabel('ω/ω1')
axs[2].set_ylabel('S(ω,phi)')
axs[2].set_title('5th Harmonic')
axs[2].text(5,.1,'Calc FWHM: ' + str(FWHMCalcH5))
plt.show()
#Mike Lafky
#USPAS June 2022
#HW1.3

from cmath import nan, pi
import numpy as np
import matplotlib.pyplot as plt
import cmath
import math
from scipy.fftpack import fft, rfft, fftfreq, rfftfreq


numE = 100                      #Number of Electrons
bunchLength = 100               #bunch Length, arbitrary time units
time = np.arange(-75,75,.1)      #Independent Axis, arbitrary time units
eField = np.zeros(len(time),dtype=complex)    #Dependent Axis, Calculated E Field at each point of time, not SI
#timej = (np.random.rand(numE)*bunchLength) #Generate a random distribution of times for temporal incoherence, between 0 and the bunchLength
timej = (np.random.rand(numE)*bunchLength)-bunchLength*.5 #Generate a random distribution of times for temporal incoherence, corrected +/- half the bunch length
#timej = np.zeros(numE)
#np.round(timej,2)              #Doesn't work, don't know why
sigmaT = 5
modes = bunchLength/(sigmaT*4)
omegaOne = 2*pi

for t in range(len(time)):  #Calculates the E field at each unit of time
    for y in range(numE):
        eField[t] = eField[t] + cmath.exp( (-((time[t]-timej[y])**2)/(4*(sigmaT**2))) - 1j*omegaOne*(time[t]-timej[y]) )
    #print(t, y, time[t], eField[t])

#Calculate the Fourier Transform
#efft = np.fft.fft(eField)[:len(eField)//2]
efft = np.fft.fft(eField)
efft = np.abs(efft)
#tf = fftfreq(len(eField),150)[:len(eField)//2]
tf = fftfreq(len(eField))

#tnum = len(eField)
#tnum2 = int(tnum/2)
#freq = time[:tnum2]/(len(time)*.1)

fig1, axs = plt.subplots(1,2,figsize=(16,9))
fig1.suptitle('')
axs[0].plot(time,eField)
axs[0].set_title('E(t) vs t')
axs[0].text(3,4,'sigmaT = ' + str(sigmaT) + ' # modes = ' + str(int(modes)))
axs[1].plot(tf,efft)
#axs[1].plot(freq,efft)
axs[1].set_title('E(t) vs s')
#plt.plot(time,eField)
#plt.plot(tf, efft)
plt.grid()
plt.show()
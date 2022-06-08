#Mike Lafky
#USPAS June 2022
#HW1.3

from cmath import nan, pi
from os import times_result
import numpy as np
import matplotlib.pyplot as plt
import random
import cmath
from scipy.optimize import curve_fit
import math

numE = 100                      #Number of Electrons
bunchLength = 100               #bunch Length, arbitrary time units
time = np.arange(-75,75,1)      #Independent Axis, arbitrary time units
eField = np.zeros(len(time))    #Dependent Axis, Calculated E Field at each point of time, not SI
#timej = (np.random.rand(numE)*bunchLength) #Generate a random distribution of times for temporal incoherence, between 0 and the bunchLength
timej = (np.random.rand(numE)*bunchLength)-bunchLength*.5 #Generate a random distribution of times for temporal incoherence, corrected +/- half the bunch length
#timej = np.zeros(numE)
#np.round(timej,2)              #Doesn't work, don't know why
sigmaT = 2
omegaOne = 2*pi

for t in range(len(time)):  #Calculates the E field at each unit of time
    for y in range(numE):
        eField[t] = eField[t] + cmath.exp( (-((time[t]-timej[y])**2)/(4*(sigmaT**2))) - 1j*omegaOne*(time[t]-timej[y]) )

#Calculate the Fourier Transform

plt.plot(time,eField)
plt.show()
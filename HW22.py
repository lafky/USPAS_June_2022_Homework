#Mike Lafky
#USPAS June 2022
#HW2.2

import math
import numpy as np
import matplotlib.pyplot as plt

x = np.arange(-5,5,.01)     #gamma*phi
dPxdPhi = np.zeros(len(x))  #arrays for power 
dPydPhi = np.zeros(len(x))

for a in range(len(x)):
    dPxdPhi = 1/(1+x**2)**(5/2)
    dPydPhi = (5*x**2)/(1+x**2)**(7/2)

fig1, axs = plt.subplots(1,2,figsize=(16,9))
fig1.suptitle('')
axs[0].plot(x,dPxdPhi)
axs[0].set_xlabel('λ*γ')
axs[0].set_ylabel('dPxdPhi')
axs[1].plot(x,dPydPhi)
axs[1].set_xlabel('λ*γ')
axs[1].set_ylabel('dPydPhi')
plt.grid()
plt.show()
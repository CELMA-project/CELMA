from boutdata import collect
import numpy as np
from scipy.signal import periodogram

import os, sys
sys.path.append("/home/mmag/CELMA-dev/celma/common/")
from CELMAPython.statsAndSignals.averages import polAvg

phi = collect("phi", yind=[16,16])
t = collect("t_array")
fs = 1/(t[1]-t[0])

phiAvg = polAvg(phi)
phiFluct = phi - phiAvg

# Actually selects one point here
# TODO: Consider polAvg <------- NO!!!! polAvg woul make these 0
#       What about the PSD...take polAvg of that?
phiFluctSlice = phiFluct[:,:,0,1]

phiPSD = []
for i in range(phiFluctSlice.shape[-1]):
    phiPSD.append(periodogram(phiFluctSlice[:,i], fs=fs, window=None, scaling="density"))

PSDFreq = np.array([amp[0] for amp in phiPSD])
PSDAmp  = np.array([amp[1] for amp in phiPSD])

# Transpose in order to have higher frequencies upwards
logPSDAmp = np.log10(PSDAmp.transpose())

import pickle
with open("/home/mmag/PSD", "wb") as f:
    pickle.dump(logPSDAmp, f)

#%%
import matplotlib.pylab as plt
import numpy as np

import pickle
with open("B0.02.pickle", "rb") as f:
    logPSDAmp = pickle.load(f)
    
#logPSDAmp = np.log10(np.exp(logPSDAmp))

# Just tried this...don't know if -20 makes sense
vmin = -20
fig, ax = plt.subplots()
ax.contourf(logPSDAmp, 100, cmap="inferno", vmax=-1, vmin=-10)
plt.show()

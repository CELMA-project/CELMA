from boutdata import collect
import numpy as np
import sys
import pickle
sys.path.append("/home/mmag/CELMA-dev/celma/common/")
from CELMAPython.statsAndSignals.averages import timeAvg, polAvg

phi = collect("phi")

phiAvgTZ = polAvg(timeAvg(phi))

tLen, xLen, yLen, zLen = phi.shape

phiTZFluct = np.zeros(phi.shape)

for t in range(tLen):
    phiTZFluct[t,:,:,:] = phi[t,:,:,:] - phiAvgTZ[0,:,:,:]

stdDev = np.sqrt(polAvg(timeAvg(phiTZFluct**2.0)))
with open("/home/mmag/phiStdDev", "wb") as f:
    pickle.dump(stdDev, f)

#%%
import matplotlib.pylab as plt

with open("phiStdDev", "rb") as f:
    phiStdDev = pickle.load(f)

# z should be the same for all, as we are taking polAvg
phiStdDev = phiStdDev[0,:,0,0]


fig, ax = plt.subplots()

ax.plot(phiStdDev)

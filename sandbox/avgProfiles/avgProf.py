# Collect stuff
#%%
%load_ext autoreload
%autoreload 2

from boutdata import collect
import numpy as np

import os, sys
sys.path.append("/home/mmag/CELMA-dev/celma/common/")
from CELMAPython.statsAndSignals.averages import timeAvg

import pickle

# Must include all zind as we are taking a pol average
lnN = collect("lnN", yind=[16,16])
t = collect("t_array")
n = np.exp(lnN)
polSliceN = n[:,16,0,:]

# nt = [n,t]
# with open("/home/mmag/nt", "wb") as f:
#     pickle.dump(nt, f)

from CELMAPython.statsAndSignals.averages import timeAvg, polAvg

tAvgN, tAvgT = timeAvg(n,t)
nAvgTZ = polAvg(tAvgN)

tLen, xLen, yLen, zLen = n.shape

nTZFluct = np.zeros(n.shape)

for t in range(tLen):
    nTZFluct[t,:,:,:] = n[t,:,:,:] - nAvgTZ[0,:,:,:]

stdDev = np.sqrt(polAvg(timeAvg(nTZFluct**2.0)))

avgStd = [nAvgTZ, stdDev]
with open("/home/mmag/avgStd.pickle", "wb") as f:
    pickle.dump(avgStd, f)

#%%
%load_ext autoreload
%autoreload 2

# Find steady state bg
from boutdata import collect
import numpy as np

import pickle

# Must include all zind as we are taking a pol average
lnN = collect("lnN", yind=[16,16])
t = collect("t_array")
n = np.exp(lnN)
steadyStateN = n[0,:,0,0]

with open("steadyStateN.pickle", "wb") as f:
    pickle.dump(steadyStateN, f)
#%%





# Plot stuff
#%%
import numpy as np
import matplotlib.pylab as plt

import pickle

with open("avgStd.pickle", "rb") as f:
    avg, std = pickle.load(f)
with open("steadyStateN.pickle", "rb") as f:
    steadyStateN = pickle.load(f)

avg = avg[0,:,0,0]
std = std[0,:,0,0]

fig, (ax1, ax2) = plt.subplots(nrows=2, sharex=True)

# Dummy
rho = range(len(avg))
ax1.plot(rho, avg)
ax1.plot(rho, steadyStateN)

ax1.fill_between(rho,\
                avg+std, avg-std,\
                facecolor="blue", edgecolor="none", alpha=0.5)

# Find the gradient
# FIXME: Need to add the gridspacing there as well, this is given as a
#        positional vararg
gradN       = np.gradient(avg, edge_order=2)
gradNSteady = np.gradient(steadyStateN, edge_order=2)
# From error propagation
gradNStd = np.sqrt((np.abs(gradN)**2)*(std**2))

ax2.plot(rho, gradN)
ax2.plot(rho, gradNSteady)
ax2.fill_between(rho,\
                gradN+gradNStd, gradN-gradNStd,\
                facecolor="blue", edgecolor="none", alpha=0.5)

plt.show()
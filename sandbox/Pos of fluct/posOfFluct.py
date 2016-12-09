from boutdata import collect
import numpy as np

import os, sys
sys.path.append("/home/mmag/CELMA-dev/celma/common/")
from CELMAPython.statsAndSignals.polAvg import polAvg

# Must include all zind as we are taking a pol average
phi = collect("phi", yind=[16,16])

phiAvg = polAvg(phi)
phiFluct = phi - phiAvg

phiRatio = phiFluct/phi
phiRatio = phiRatio[:,:,0,0]
import pickle
with open("/home/mmag/ratio", "wb") as f:
    pickle.dump(phiRatio, f)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
from matplotlib import animation
import matplotlib.pylab as plt
import numpy as np

import pickle
with open("ratio", "rb") as f:
    phiRatio = pickle.load(f)

phiMax = np.max(phiRatio[:600,:])
phiMin = np.min(phiRatio[:600,:])

fig, ax = plt.subplots()

# Also have real position here
def animRatio(t):
    ax.cla()
    ax.plot(phiRatio[t,:])
    ax.set_ylim([phiMin, phiMax])
    ax.set_xlabel(r"rho")
    ax.set_ylabel(r"level")
    ax.text(0.98,0.98,r"tind = {}".format(t),
            ha="right",va="top",transform=ax.transAxes)

anim = animation.FuncAnimation(fig, animRatio, frames = phiRatio.shape[0], interval = 1)
plt.show()
# anim.save("anim.avi",fps=5)

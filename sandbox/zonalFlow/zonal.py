from boutdata import collect
import numpy as np
import pickle

phi = collect("phi", yind=[16,16])
dx  = collect("dx")
dx = dx[0,0]
# Background is the zeroth mode: Should polAvg, but as modes propto y => background mode is the only one not having exp(y) <- integrated out
phiFFT = np.fft.fft(phi, axis=-1)

#Filter
phiFFT[:,:,:,1:] = 0.0 + 0.0j
# Should be no variation in theta as we have filtered the signal (i.e.
# all z points will on a radius will be identical)
phiFiltered = np.fft.ifft(phiFFT)[:,:,0,0]

radGradPhi = np.gradient(phi, dx, edge_order=2)
vPol       = np.zeros(phiFiltered.shape)
for t in range(phiFiltered.shape[0]):
    # remember the B field in non-normalized
    vPol[t, :] = np.gradient(phiFiltered[t,:].real, dx, edge_order=2)

# om x rho = v => om = v/rho, should divide by rho, but index number suffice for now (in the end will just have wrong scaling)
angVPol = vPol/np.array(range(1,vPol.shape[-1]+1))

with open("/home/mmag/angVPol", "wb") as f:
    pickle.dump(angVPol, f)

#%%
import pickle
import matplotlib.pylab as plt
import numpy as np
from matplotlib import animation

with open("angVPol", "rb") as f:
    angVPol = pickle.load(f)
        
maxAngVPol = np.max(angVPol)
minAngVPol = np.min(angVPol)

def animateVPol(t):
    ax.cla()
    ax.plot(angVPol[t,:])
    ax.set_ylim((minAngVPol, maxAngVPol))

    ax.set_xlabel(r"rho")
    ax.set_ylabel(r"ang freq")
    
    ax.text(0.98,0.98,r"t={}".format(t),\
            ha="right",va="top",transform=ax.transAxes)

fig, ax = plt.subplots()
    
anim = animation.FuncAnimation(\
            fig, animateVPol, frames = angVPol.shape[0], interval = 1)
anim.save("anim.mp4",fps=5)
#plt.show()

#%%
import pickle
import matplotlib.pylab as plt
import numpy as np
from matplotlib import animation

with open("angVPol", "rb") as f:
    angVPol = pickle.load(f)

vPol = angVPol*np.array(range(1,angVPol.shape[-1]+1))

maxVPol = np.max(vPol)
minVPol = np.min(vPol)

def animateVPol(t):
    ax.cla()
    ax.plot(vPol[t,:])
    ax.set_ylim((minVPol, maxVPol))

    ax.set_xlabel(r"rho")
    ax.set_ylabel(r"ang freq")
    
    ax.text(0.98,0.98,r"t={}".format(t),\
            ha="right",va="top",transform=ax.transAxes)

fig, ax = plt.subplots()
    
anim = animation.FuncAnimation(\
            fig, animateVPol, frames = angVPol.shape[0], interval = 1)
#anim.save("anim.mp4",fps=5)
plt.show()
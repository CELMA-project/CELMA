from boutdata import collect
from boututils.datafile import DataFile
f = DataFile("BOUT.dmp.0.nc")
f.list()

colMe=[
 'parKinEE',
 'parKinEI',
 'perpKinEE',
 'perpKinEI',
 'polAvgParKinEE',
 'polAvgParKinEI',
 'polAvgPerpKinEE',
 'polAvgPerpKinEI',
 'polAvgSumKinEE',
 'polAvgSumKinEI',
 'sumKinEE',
 'sumKinEI',
 'polAvgPotEE',
 'potEE',
 'particleNumber',
]
d = {}
for col in colMe:
    d[col] = collect(col)
import pickle
with open("/home/mmag/ens", "wb") as f:
    pickle.dump(d, f)

#%%

# !!!!!!!! THIS SHOWS PREDATOR PREY BEHAVIOUR: 
#          SHOULD LOOK FOR ZONAL FLOWS
#          Pot E rise => kin E rise => Pot E decline = kin E decline => pot E rises

import numpy as np
import matplotlib.pylab as plt
import pickle
with open("ens", "rb") as f:
        enDict = pickle.load(f)

def renormSignalToZeroAndMax(var):
    """
    Subtracts the min value, so that the signal min is at 0.
    Normalizes the resulting signal to 1
    """
    
    putToZero = var - np.min(var)
    return putToZero/np.max(putToZero)
    
normKinEE = renormSignalToZeroAndMax(enDict["sumKinEE"])
normKinEI = renormSignalToZeroAndMax(enDict["sumKinEI"])
# normal pot E
normParticleNumber = renormSignalToZeroAndMax(enDict["particleNumber"])
# jmads pot E
normPotEE = renormSignalToZeroAndMax(enDict["potEE"])

# Plot
fig, ax = plt.subplots()

ax.plot(normKinEE)
ax.plot(normKinEI)
ax.plot(normParticleNumber)
ax.plot(normPotEE)

# Grid
ax.grid()

# Show
plt.show()
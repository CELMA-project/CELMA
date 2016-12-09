from boutdata import collect
import numpy as np

lnN = collect("lnN", yind=[16,16])

n = np.exp(lnN)
polSliceN = n[:,16,16,:]
polSliceN = n[:,16,0,:]
import pickle
with open("/home/mmag/polSliceN", "wb") as f:
    pickle.dump(polSliceN, f)


#%%
import matplotlib.pylab as plt
import numpy as np

import pickle
with open("polSliceN", "rb") as f:
    polSliceN = pickle.load(f)

fig, ax = plt.subplots()
ax.contourf(polSliceN.transpose(), 100, cmap="inferno")
plt.show()
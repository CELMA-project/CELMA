%load_ext autoreload
%autoreload 2
commonDir = "/home/mmag/CELMA-dev/celmaRC/common/"

# Check that the max is working
from CELMAPython.calcHelpers import collectRadialProfileTime
nRad = collectRadialProfileTime(("../../nout_2_timestep_25/nz_256/geom_Lx_7.8633_geom_Ly_275.2144_input_B0_0.1_ownFilters_type_none_switch_useHyperViscAzVortD_False_tag_CSDXMagFieldScanAr-1-expand_0/",), "n", 16, 0, tInd=[2,2])

# FIXME: Should fix this
from boutdata import collect
dx = collect("dx", path="../../nout_2_timestep_25/nz_256/geom_Lx_7.8633_geom_Ly_275.2144_input_B0_0.1_ownFilters_type_none_switch_useHyperViscAzVortD_False_tag_CSDXMagFieldScanAr-1-expand_0/")

from CELMAPython.calcHelpers import DDX
ddxN = DDX(nProf,dx[0,0]).flatten()

from CELMAPython.calcHelpers import findLargestRadialGrad
largest, ind  = findLargestRadialGrad(nProf, dx[0,0])

n = nProf.flatten()

import pickle
with open("/home/mmag/nDdxNIndLargest.pickle", "wb") as f:
    pickle.dump((n,ddxN,ind,largest), f)

#%%
import pickle
with open("nDdxNIndLargest.pickle","rb") as f:
    n, ddxN, ind, largest = pickle.load(f)
import matplotlib.pylab as plt
fig, (ax1, ax2) = plt.subplots(nrows=2, sharex=True)
ax1.plot(n)
ax2.plot(ddxN)
ax1.grid()
ax2.grid()
plt.show()
print(ind, largest)

# Seen that above works :D






#%%
%load_ext autoreload
%autoreload 2
import sys
commonDir = "/home/mmag/CELMA-dev/celmaRC/common/"
sys.path.append(commonDir)
from CELMAPython.superClasses import PointsSuperClass
steadyStatePath = "../../nout_2_timestep_25/nz_256/geom_Lx_7.8633_geom_Ly_275.2144_input_B0_0.1_ownFilters_type_none_switch_useHyperViscAzVortD_False_tag_CSDXMagFieldScanAr-1-expand_0/"
xInd = None
yInd = 16
zInd = 0
nPoints = 5
pc = PointsSuperClass("."                 ,
                 xInd                     ,
                 yInd                     ,
                 zInd                     ,
                 varName         = "n"    ,
                 mode            = "fluct",
                 equallySpace    = "x"    ,
                 nPoints         = None   ,
                 steadyStatePath = steadyStatePath,
                 )












COMMON CLASS ABOVE IS SUPERCLASS FOR THE POSTBOUT DRIVERS!

















# Get the time trace in the end
from CELMAPython.timeTrace import calcTimeTrace, calcTimeTrace4d
# 14 from max grad, 16 for fun, 0 standard
calcTimeTrace((".",), "n", xInd=14, yInd=16, zInd=0, convertToPhysical=True, mode='fluct')

%load_ext autoreload
%autoreload 2
import sys
commonDir = "/home/mmag/CELMA-dev/celmaRC/common/"
sys.path.append(commonDir)
from CELMAPython.fields1D import CollectAndCalcFields1D

paths = (".",)
convertToPhysical = True
mode = "parallel"
processing = None
ccf1D = CollectAndCalcFields1D(paths,\
                               mode = mode,\
                               processing = processing,\
                               convertToPhysical = convertToPhysical)

xSlice = 16
ySlice = None
zSlice = 0
tSlice = None
ccf1D.setSlice(xSlice,\
               ySlice,\
               zSlice,\
               tSlice)
ccf1D.setVarName("lnN")
par1D = ccf1D.executeCollectAndCalc()

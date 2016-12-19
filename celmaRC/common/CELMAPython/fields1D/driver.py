%load_ext autoreload
%autoreload 2
import sys
commonDir = "/home/mmag/CELMA-dev/celmaRC/common/"
sys.path.append(commonDir)
from CELMAPython.fields1D import CollectAndCalcFields1D
from CELMAPython.collectAndCalcHelpers import calcN, calcUIPar, calcUE

mainFields  = ("lnN"       ,\
               "vort"      ,\
               "jPar"      ,\
               "phi"       ,\
               "vort"      ,\
               "vortD"     ,\
               "momDensPar",\
               "S"         ,\
              )

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

par1D = {}

for field in mainFields:
    ccf1D.setVarName(field)
    par1D.update(ccf1D.executeCollectAndCalc())

# Non-collects
par1D.update({"n"    : calcN(par1D["lnN"])})
par1D.update({"uIPar": calcUIPar(par1D["momDensPar"], par1D["n"])})
par1D.update({"uEPar": calcUEPar(par1D["uIPar"], par1D["jPar"], par1D["n"], not(ccf1D.convertToPhysical))})

# FIXME: Utested!
plotOrder = ("lnN"  , "phi"       ,\
             "n"    , "vortD"     ,\
             "jPar" , "momDensPar",\
             "uIPar", "uIPar"     ,\
             "uEPar", "uEPar"     ,\
             "S")

p1DPar = Plot1DPar(".", ccf1D.convertToPhysical)
p1DPar.setParData(par1D)
p1DPar.setPlotOrder(plotOrder)
p1DPar.plotAndSavePar(plotOrder)

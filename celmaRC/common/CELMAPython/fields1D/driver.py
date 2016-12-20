%load_ext autoreload
%autoreload 2
import sys
commonDir = "/home/mmag/CELMA-dev/celmaRC/common/"
sys.path.append(commonDir)
from CELMAPython.fields1D import CollectAndCalcFields1D, PlotAnim1DRadial
from CELMAPython.collectAndCalcHelpers import calcN, calcUIPar, calcUEPar

mainFields  = ("lnN"       ,\
               "jPar"      ,\
               "phi"       ,\
               "vort"      ,\
               "vortD"     ,\
               "momDensPar",\
               "S"         ,\
              )

paths = (".",)
convertToPhysical = True
mode = "radial"
processing = None
ccf1D = CollectAndCalcFields1D(paths,\
                               mode = mode,\
                               processing = processing,\
                               convertToPhysical = convertToPhysical)

xSlice = None
ySlice = 16
zSlice = 0
tSlice = None
ccf1D.setSlice(xSlice,\
               ySlice,\
               zSlice,\
               tSlice)

radial1D = {}

for field in mainFields:
    ccf1D.setVarName(field)
    radial1D.update(ccf1D.executeCollectAndCalc())

# Non-collects
radial1D.update({"n"    : calcN(radial1D["lnN"], not(ccf1D.convertToPhysical), ccf1D.uc)})
radial1D.update({"uIPar": calcUIPar(radial1D["momDensPar"], radial1D["n"])})
radial1D.update({"uEPar": calcUEPar(radial1D["uIPar"], radial1D["jPar"], radial1D["n"], not(ccf1D.convertToPhysical))})

plotOrder = ("lnN"  , "phi"       ,\
             "n"    , "vortD"     ,\
             "jPar" , "vort"      ,\
             "uIPar", "momDensPar",\
             "uEPar", "S"         ,\
            )

p1DRadial = PlotAnim1DRadial(".", ccf1D.convertToPhysical)
p1DRadial.setRadialData(radial1D, "mainFields", ".", plotOrder=plotOrder)
p1DRadial.plotAndSaveRadialProfile()










#%%
%load_ext autoreload
%autoreload 2
import sys
commonDir = "/home/mmag/CELMA-dev/celmaRC/common/"
sys.path.append(commonDir)
from CELMAPython.fields1D import CollectAndCalcFields1D, PlotAnim1DParallel
from CELMAPython.collectAndCalcHelpers import calcN, calcUIPar, calcUEPar

mainFields  = ("lnN"       ,\
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

xSlice = 0
ySlice = None
zSlice = 0
tSlice = None
ccf1D.setSlice(xSlice,\
               ySlice,\
               zSlice,\
               tSlice)

parallel1D = {}

for field in mainFields:
    ccf1D.setVarName(field)
    parallel1D.update(ccf1D.executeCollectAndCalc())

# Non-collects
parallel1D.update({"n"    : calcN(parallel1D["lnN"], not(ccf1D.convertToPhysical), ccf1D.uc)})
parallel1D.update({"uIPar": calcUIPar(parallel1D["momDensPar"], parallel1D["n"])})
parallel1D.update({"uEPar": calcUEPar(parallel1D["uIPar"], parallel1D["jPar"], parallel1D["n"], not(ccf1D.convertToPhysical))})

plotOrder = ("lnN"  , "phi"       ,\
             "n"    , "vortD"     ,\
             "jPar" , "vort"      ,\
             "uIPar", "momDensPar",\
             "uEPar", "S"         ,\
            )

p1DParallel = PlotAnim1DParallel(".", ccf1D.convertToPhysical)
p1DParallel.setParallelData(parallel1D, "mainFields", ".", plotOrder=plotOrder)
p1DParallel.plotAndSaveParallelProfile()

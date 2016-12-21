# FIXME: Write a proper driver based on this

# FIXME: Name changes
# FIXME: ConvertToPhysical




# # FIXME: Move maxGradRho to field2D class, and make it a ghost there
#         # Get the current scan
#         if maxGradRhoFolder:
#             if self._scanParameters:
#                 self._maxGradRhoFolder = maxGradRhoFolder
#             else:
#                 self._maxGradRhoFolder =\
#                     convertToCurrentScanParameters(dmpFolder, maxGradRhoFolder, scanParameters)
#         else:
#             self._maxGradRhoFolder = None






%load_ext autoreload
%autoreload 2
import sys
commonDir = "/home/mmag/CELMA-dev/celmaRC/common/"
sys.path.append(commonDir)
from CELMAPy.fields2D import CollectAndCalcFields2D, PlotAnim2DPerpPol
from CELMAPy.plotHelpers import getMaxMinAnimation, getLevelsAnimation

fluct = True
convertToPhysical = True
ccF2D = CollectAndCalcFields2D(".", "lnN", mode="perp", xSlice=None, zSlice=None, ySlice=16, fluct=fluct, convertToPhysical=convertToPhysical)
perp2D = ccF2D.executeCollectAndCalc()

ccF2D = CollectAndCalcFields2D(".", "lnN", mode="pol", xSlice=16, zSlice=None, ySlice=None, fluct=fluct, convertToPhysical=convertToPhysical)
pol2D = ccF2D.executeCollectAndCalc()

vmax, vmin = getMaxMinAnimation((perp2D["lnN"], pol2D["lnN"]), True, True)
levels = getLevelsAnimation(vmax, vmin, 100)

p2DPerpPol = PlotAnim2DPerpPol(".", fluct, convertToPhysical)
p2DPerpPol.setContourfArguments(vmax, vmin, levels)
p2DPerpPol.setPerpData(perp2D["X"], perp2D["Y"], perp2D["lnN"], perp2D["time"], perp2D["zPos"], "lnN", ".")
p2DPerpPol.setPolData(pol2D["X"], pol2D["Y"], pol2D["lnN"], pol2D["time"], pol2D["rhoPos"], "lnN", ".")
p2DPerpPol.plotAndSavePerpPlane()




#%%
%load_ext autoreload
%autoreload 2
import sys
commonDir = "/home/mmag/CELMA-dev/celmaRC/common/"
sys.path.append(commonDir)
from CELMAPython.fields2D import CollectAndCalcFields2D, PlotAnim2DPerpPar
from CELMAPython.plotHelpers import getMaxMinAnimation, getLevelsAnimation

fluct = True
convertToPhysical = True
ccF2D = CollectAndCalcFields2D(".", "lnN", mode="perp", xSlice=None, zSlice=None, ySlice=16, fluct=fluct, convertToPhysical=convertToPhysical)
perp2D = ccF2D.executeCollectAndCalc()

ccF2D = CollectAndCalcFields2D(".", "lnN", mode="par", xSlice=None, zSlice=0, ySlice=None, fluct=fluct, convertToPhysical=convertToPhysical)
par2D = ccF2D.executeCollectAndCalc()

vmax, vmin = getMaxMinAnimation((perp2D["lnN"], par2D["lnN"], par2D["lnNPPi"]), True, True)
levels = getLevelsAnimation(vmax, vmin, 100)

p2DPerpPar = PlotAnim2DPerpPar(".", fluct, convertToPhysical)
p2DPerpPar.setContourfArguments(vmax, vmin, levels)
p2DPerpPar.setPerpData(perp2D["X"], perp2D["Y"], perp2D["lnN"], perp2D["time"], perp2D["zPos"], "lnN", ".")
p2DPerpPar.setParData(par2D["X"], par2D["Y"], par2D["lnN"], par2D["lnNPPi"], par2D["time"], par2D["thetaPos"], "lnN", ".")
p2DPerpPar.plotAndSavePerpPlane()




#%%
%load_ext autoreload
%autoreload 2
import sys
commonDir = "/home/mmag/CELMA-dev/celmaRC/common/"
sys.path.append(commonDir)
from CELMAPython.fields2D import CollectAndCalcFields2D, PlotAnim2DPol
from CELMAPython.plotHelpers import getMaxMinAnimation, getLevelsAnimation

fluct = True
convertToPhysical = True
ccF2D = CollectAndCalcFields2D(".", "lnN", mode="pol", xSlice=16, zSlice=None, ySlice=None, fluct=fluct, convertToPhysical=convertToPhysical)
pol2D = ccF2D.executeCollectAndCalc()
vmax, vmin = getMaxMinAnimation((pol2D["lnN"],), True, True)
levels = getLevelsAnimation(vmax, vmin, 100)
p2DPol = PlotAnim2DPol(".", fluct, convertToPhysical)
p2DPol.setContourfArguments(vmax, vmin, levels)
p2DPol.setPolData(pol2D["X"], pol2D["Y"], pol2D["lnN"], pol2D["time"], pol2D["rhoPos"], "lnN", ".")
p2DPol.plotAndSavePolPlane()




#%%
%load_ext autoreload
%autoreload 2
import sys
commonDir = "/home/mmag/CELMA-dev/celmaRC/common/"
sys.path.append(commonDir)
from CELMAPython.fields2D import CollectAndCalcFields2D, PlotAnim2DPar
from CELMAPython.plotHelpers import getMaxMinAnimation, getLevelsAnimation

fluct = True
convertToPhysical = True
ccF2D = CollectAndCalcFields2D(".", "lnN", mode="par", xSlice=None, zSlice=0, ySlice=None, fluct=fluct, convertToPhysical=convertToPhysical)
par2D = ccF2D.executeCollectAndCalc()
vmax, vmin = getMaxMinAnimation((par2D["lnN"], par2D["lnNPPi"]), True, True)
levels = getLevelsAnimation(vmax, vmin, 100)
p2DPar = PlotAnim2DPar(".", fluct, convertToPhysical)
p2DPar.setContourfArguments(vmax, vmin, levels)
p2DPar.setParData(par2D["X"], par2D["Y"], par2D["lnN"], par2D["lnNPPi"], par2D["time"], par2D["thetaPos"], "lnN", ".")
p2DPar.plotAndSaveParPlane()




#%%
%load_ext autoreload
%autoreload 2
import sys
commonDir = "/home/mmag/CELMA-dev/celmaRC/common/"
sys.path.append(commonDir)
from CELMAPython.fields2D import CollectAndCalcFields2D, PlotAnim2DPerp
from CELMAPython.plotHelpers import getMaxMinAnimation, getLevelsAnimation

fluct = True
convertToPhysical = True
ccF2D = CollectAndCalcFields2D(".", "lnN", mode="perp", xSlice=None, zSlice=None, ySlice=16, fluct=fluct, convertToPhysical=convertToPhysical)
perp2D = ccF2D.executeCollectAndCalc()
vmax, vmin = getMaxMinAnimation((perp2D["lnN"],), True, True)
levels = getLevelsAnimation(vmax, vmin, 100)
p2DPerp = PlotAnim2DPerp(".", fluct, convertToPhysical)
p2DPerp.setContourfArguments(vmax, vmin, levels)
p2DPerp.setPerpData(perp2D["X"], perp2D["Y"], perp2D["lnN"], perp2D["time"], perp2D["zPos"], "lnN", ".")
p2DPerp.plotAndSavePerpPlane()

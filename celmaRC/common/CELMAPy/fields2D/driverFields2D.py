#!/usr/bin/env python

"""
Contains single driver and driver class for 2D fields
"""

from ..plotHelpers import getMaxMinAnimation, getLevelsAnimation
from ..superClasses import DriverPlotFieldsSuperClass
from ..collectAndCalcHelpers import findLargestRadialGradN
from .collectAndCalcFields2D import CollectAndCalcFields2D
from .plotFields2D import (PlotAnim2DPerp,\
                           PlotAnim2DPar,\
                           PlotAnim2DPol,\
                           PlotAnim2DPerpPar,\
                           PlotAnim2DPerpPol,\
                           )
from multiprocessing import Process
import os

#{{{driver2DFieldPerpSingle
def driver2DFieldPerpSingle(collectPaths     ,\
                            savePath         ,\
                            varName          ,\
                            convertToPhysical,\
                            xSlice           ,\
                            yInd             ,\
                            zSlice           ,\
                            tSlice           ,\
                            fluct            ,\
                            varyMaxMin       ,\
                            xguards  = False ,\
                            yguards  = False ,\
                            showPlot = False ,\
                            savePlot = True  ,\
                                            ):
    #{{{doctring
    """
    Driver for plotting a single perpendicular 2D plot

    Parameters
    ----------
    collectPaths : tuple
        Paths to collect from.
        The corresponind 't_array' of the paths must be in ascending order.
    savePath : str
        Save destination
    varName : str
        The variable to plot
    convertToPhysical : bool
        Whether or not to convert to physical units.
    xSlice : slice
        How the data will be sliced in x.
    yInd : int
        How the data will be sliced in y.
    zSlice : slice
        How the data will be sliced in z.
    tSlice : [None|slice]
        How the data will be sliced in t.
    fluct : bool
        Whether or not fluctuations are shown
    varyMaxMin : bool
        If the colorbar should be adjusted to the max and min of the
        current timestep.
        If False, the global max and min is used.
    xguards : bool
        Whether or not to collect the ghost points in x.
    yguards : bool
        Whether or not to collect the ghost points in y.
    showPlot : bool
        Whether or not the plot should be displayed.
    savePlot : bool
        Whether or no the plot should be saved.
    """
    #}}}

    # Declare mode
    mode = "perp"

    # Create collect object
    ccf2D = CollectAndCalcFields2D(collectPaths             ,\
                                   fluct             = fluct,\
                                   mode              = mode ,\
                                   convertToPhysical = convertToPhysical)

    # Set the slice
    ccf2D.setSlice(xSlice, yInd, zSlice, tSlice)
    # Set name and execute
    ccf2D.setVarName(varName)
    # Execute the collection
    perp2D = ccf2D.executeCollectAndCalc()

    # Set the plot limits
    vmax, vmin = getMaxMinAnimation((perp2D[varName],),\
                                     fluct,
                                     varyMaxMin)
    levels = getLevelsAnimation(vmax, vmin, 100)

    # Create the plotting object
    p2DPerp = PlotAnim2DPerp(collectPaths           ,\
                             savePath               ,\
                             ccf2D.convertToPhysical,\
                             fluct = fluct          ,\
                             show  = showPlot       ,\
                             save  = savePlot)
    p2DPerp.setContourfArguments(vmax, vmin, levels)
    p2DPerp.setPerpData(perp2D["X"],\
                        perp2D["Y"],\
                        perp2D[varName],\
                        perp2D["time"],\
                        perp2D["zPos"],\
                        varName)
    p2DPerp.plotAndSavePerpPlane()
#}}}

#{{{driver2DFieldParSingle
def driver2DFieldParSingle(collectPaths     ,\
                           savePath         ,\
                           varName          ,\
                           convertToPhysical,\
                           xSlice           ,\
                           ySlice           ,\
                           zInd             ,\
                           tSlice           ,\
                           fluct            ,\
                           varyMaxMin       ,\
                           xguards  = False ,\
                           yguards  = False ,\
                           showPlot = False ,\
                           savePlot = True  ,\
                                            ):
    #{{{doctring
    """
    Driver for plotting a single parallel 2D plot

    Parameters
    ----------
    collectPaths : tuple
        Paths to collect from.
        The corresponind 't_array' of the paths must be in ascending order.
    savePath : str
        Save destination
    varName : str
        The variable to plot
    convertToPhysical : bool
        Whether or not to convert to physical units.
    xSlice : slice
        How the data will be sliced in x.
    ySlice : int
        How the data will be sliced in y.
    zInd : int
        How the data will be sliced in z.
    tSlice : [None|slice]
        How the data will be sliced in t.
    fluct : bool
        Whether or not fluctuations are shown
    varyMaxMin : bool
        If the colorbar should be adjusted to the max and min of the
        current timestep.
        If False, the global max and min is used.
    xguards : bool
        Whether or not to collect the ghost points in x.
    yguards : bool
        Whether or not to collect the ghost points in y.
    showPlot : bool
        Whether or not the plot should be displayed.
    savePlot : bool
        Whether or no the plot should be saved.
    """
    #}}}

    mode = "par"

    # Create collect object
    ccf2D = CollectAndCalcFields2D(collectPaths             ,\
                                   fluct             = fluct,\
                                   mode              = mode ,\
                                   convertToPhysical = convertToPhysical)

    # Set the slice
    ccf2D.setSlice(xSlice, ySlice, zInd, tSlice)
    # Set name and execute
    ccf2D.setVarName(varName)
    # Execute the collection
    par2D = ccf2D.executeCollectAndCalc()

    # Set the plot limits
    vmax, vmin = getMaxMinAnimation((par2D[varName],),\
                                     fluct,
                                     varyMaxMin)
    levels = getLevelsAnimation(vmax, vmin, 100)

    # Create the plotting object
    p2DPar = PlotAnim2DPar(collectPaths           ,\
                           savePath               ,\
                           ccf2D.convertToPhysical,\
                           fluct = fluct          ,\
                           show  = showPlot       ,\
                           save  = savePlot)
    p2DPar.setContourfArguments(vmax, vmin, levels)
    p2DPar.setParData(par2D["X"]          ,\
                      par2D["Y"]          ,\
                      par2D[varName]      ,\
                      par2D[varName+"PPi"],\
                      par2D["time"]       ,\
                      par2D["thetaPos"]   ,\
                      varName)
    p2DPar.plotAndSaveParPlane()
#}}}

#{{{driver2DFieldPolSingle
def driver2DFieldPolSingle(collectPaths     ,\
                           savePath         ,\
                           varName          ,\
                           convertToPhysical,\
                           xInd             ,\
                           ySlice           ,\
                           zSlice           ,\
                           tSlice           ,\
                           fluct            ,\
                           varyMaxMin       ,\
                           xguards  = False ,\
                           yguards  = False ,\
                           showPlot = False ,\
                           savePlot = True  ,\
                                            ):
    #{{{doctring
    """
    Driver for plotting a single poloidal 2D plot

    Parameters
    ----------
    collectPaths : tuple
        Paths to collect from.
        The corresponind 't_array' of the paths must be in ascending order.
    savePath : str
        Save destination
    varName : str
        The variable to plot
    convertToPhysical : bool
        Whether or not to convert to physical units.
    xSlice : slice
        How the data will be sliced in x.
    ySlice : slice
        How the data will be sliced in y.
    zInd : int
        How the data will be sliced in z.
    tSlice : [None|slice]
        How the data will be sliced in t.
    fluct : bool
        Whether or not fluctuations are shown
    varyMaxMin : bool
        If the colorbar should be adjusted to the max and min of the
        current timestep.
        If False, the global max and min is used.
    xguards : bool
        Whether or not to collect the ghost points in x.
    yguards : bool
        Whether or not to collect the ghost points in y.
    showPlot : bool
        Whether or not the plot should be displayed.
    savePlot : bool
        Whether or no the plot should be saved.
    """
    #}}}

    mode = "pol"

    # Create collect object
    ccf2D = CollectAndCalcFields2D(collectPaths             ,\
                                   fluct             = fluct,\
                                   mode              = mode ,\
                                   convertToPhysical = convertToPhysical)

    # Set the slice
    ccf2D.setSlice(xInd, ySlice, zSlice, tSlice)
    # Set name and execute
    ccf2D.setVarName(varName)
    # Execute the collection
    pol2D = ccf2D.executeCollectAndCalc()

    # Set the plot limits
    vmax, vmin = getMaxMinAnimation((pol2D[varName],),\
                                     fluct,
                                     varyMaxMin)
    levels = getLevelsAnimation(vmax, vmin, 100)

    # Create the plotting object
    p2DPol = PlotAnim2DPol(collectPaths           ,\
                           savePath               ,\
                           ccf2D.convertToPhysical,\
                           fluct = fluct          ,\
                           show  = showPlot       ,\
                           save  = savePlot)
    p2DPol.setContourfArguments(vmax, vmin, levels)
    p2DPol.setPolData(pol2D["X"]     ,\
                      pol2D["Y"]     ,\
                      pol2D[varName] ,\
                      pol2D["time"]  ,\
                      pol2D["rhoPos"],\
                      varName)
    p2DPol.plotAndSavePolPlane()
#}}}

#{{{driver2DFieldPerpParSingle
def driver2DFieldPerpParSingle(collectPaths     ,\
                               savePath         ,\
                               varName          ,\
                               convertToPhysical,\
                               xSlice           ,\
                               ySlice           ,\
                               zSlice           ,\
                               yInd             ,\
                               zInd             ,\
                               tSlice           ,\
                               fluct            ,\
                               varyMaxMin       ,\
                               xguards  = False ,\
                               yguards  = False ,\
                               showPlot = False ,\
                               savePlot = True  ,\
                                                ):
    #{{{doctring
    """
    Driver for plotting a single perpendicular and parallel 2D plot

    Parameters
    ----------
    collectPaths : tuple
        Paths to collect from.
        The corresponind 't_array' of the paths must be in ascending order.
    savePath : str
        Save destination
    varName : str
        The variable to plot
    convertToPhysical : bool
        Whether or not to convert to physical units.
    xSlice : slice
        How the data will be sliced in x.
    ySlice : slice
        How the data will be sliced in y.
    zSlice : slice
        How the data will be sliced in z.
    yInd : int
        Constant y index for perpendicular plot.
    zInd : int
        Constant z index for parallel plot.
    tSlice : [None|slice]
        How the data will be sliced in t.
    fluct : bool
        Whether or not fluctuations are shown
    varyMaxMin : bool
        If the colorbar should be adjusted to the max and min of the
        current timestep.
        If False, the global max and min is used.
    xguards : bool
        Whether or not to collect the ghost points in x.
    yguards : bool
        Whether or not to collect the ghost points in y.
    showPlot : bool
        Whether or not the plot should be displayed.
    savePlot : bool
        Whether or no the plot should be saved.
    """
    #}}}

    # Pependicular collection
    ccf2D = CollectAndCalcFields2D(collectPaths              ,\
                                   fluct             = fluct ,\
                                   mode              = "perp",\
                                   convertToPhysical = convertToPhysical)
    ccf2D.setSlice(xSlice, yInd, zSlice, tSlice)
    ccf2D.setVarName(varName)
    perp2D = ccf2D.executeCollectAndCalc()

    # Parallel collection
    ccf2D = CollectAndCalcFields2D(collectPaths             ,\
                                   fluct             = fluct,\
                                   mode              = "par",\
                                   convertToPhysical = convertToPhysical)
    ccf2D.setSlice(xSlice, ySlice, zInd, tSlice)
    ccf2D.setVarName(varName)
    par2D = ccf2D.executeCollectAndCalc()

    # Set the plot limits
    vmax, vmin = getMaxMinAnimation((perp2D[varName],\
                                     par2D[varName]),\
                                     fluct,
                                     varyMaxMin)
    levels = getLevelsAnimation(vmax, vmin, 100)

    # Create the plotting object
    p2DPerpPar = PlotAnim2DPerpPar(collectPaths           ,\
                                   savePath               ,\
                                   ccf2D.convertToPhysical,\
                                   fluct = fluct          ,\
                                   show  = showPlot       ,\
                                   save  = savePlot)
    p2DPerpPar.setContourfArguments(vmax, vmin, levels)
    p2DPerpPar.setPerpData(perp2D["X"]    ,\
                           perp2D["Y"]    ,\
                           perp2D[varName],\
                           perp2D["time"] ,\
                           perp2D["zPos"] ,\
                           varName)
    p2DPerpPar.setParData(par2D["X"]          ,\
                          par2D["Y"]          ,\
                          par2D[varName]      ,\
                          par2D[varName+"PPi"],\
                          par2D["time"]       ,\
                          par2D["thetaPos"]   ,\
                          varName)
    p2DPerpPar.plotAndSavePerpParPlane()
#}}}

#{{{driver2DFieldPerpPolSingle
def driver2DFieldPerpPolSingle(collectPaths     ,\
                               savePath         ,\
                               varName          ,\
                               convertToPhysical,\
                               xSlice           ,\
                               ySlice           ,\
                               zSlice           ,\
                               xInd             ,\
                               yInd             ,\
                               tSlice           ,\
                               fluct            ,\
                               varyMaxMin       ,\
                               xguards  = False ,\
                               yguards  = False ,\
                               showPlot = False ,\
                               savePlot = True  ,\
                                                ):
    #{{{doctring
    """
    Driver for plotting a single perpendicular and parallel 2D plot

    Parameters
    ----------
    collectPaths : tuple
        Paths to collect from.
        The corresponind 't_array' of the paths must be in ascending order.
    savePath : str
        Save destination
    varName : str
        The variable to plot
    convertToPhysical : bool
        Whether or not to convert to physical units.
    xSlice : slice
        How the data will be sliced in x.
    ySlice : slice
        How the data will be sliced in y.
    zSlice : slice
        How the data will be sliced in z.
    xInd : int
        Constant x index for poloidal plot.
    yInd : int
        Constant y index for perpendicular plot.
    tSlice : [None|slice]
        How the data will be sliced in t.
    fluct : bool
        Whether or not fluctuations are shown
    varyMaxMin : bool
        If the colorbar should be adjusted to the max and min of the
        current timestep.
        If False, the global max and min is used.
    xguards : bool
        Whether or not to collect the ghost points in x.
    yguards : bool
        Whether or not to collect the ghost points in y.
    showPlot : bool
        Whether or not the plot should be displayed.
    savePlot : bool
        Whether or no the plot should be saved.
    """
    #}}}

    # Pependicular collection
    ccf2D = CollectAndCalcFields2D(collectPaths              ,\
                                   fluct             = fluct ,\
                                   mode              = "perp",\
                                   convertToPhysical = convertToPhysical)
    ccf2D.setSlice(xSlice, yInd, zSlice, tSlice)
    ccf2D.setVarName(varName)
    perp2D = ccf2D.executeCollectAndCalc()

    # Poloidal collection
    ccf2D = CollectAndCalcFields2D(collectPaths             ,\
                                   fluct             = fluct,\
                                   mode              = "pol",\
                                   convertToPhysical = convertToPhysical)
    ccf2D.setSlice(xInd, ySlice, zSlice, tSlice)
    ccf2D.setVarName(varName)
    pol2D = ccf2D.executeCollectAndCalc()

    # Set the plot limits
    vmax, vmin = getMaxMinAnimation((perp2D[varName],\
                                     pol2D[varName]),\
                                     fluct,
                                     varyMaxMin)
    levels = getLevelsAnimation(vmax, vmin, 100)

    # Create the plotting object
    p2DPerpPol = PlotAnim2DPerpPol(collectPaths           ,\
                                   savePath               ,\
                                   ccf2D.convertToPhysical,\
                                   fluct = fluct          ,\
                                   show  = showPlot       ,\
                                   save  = savePlot)
    p2DPerpPol.setContourfArguments(vmax, vmin, levels)
    p2DPerpPol.setPerpData(perp2D["X"]    ,\
                           perp2D["Y"]    ,\
                           perp2D[varName],\
                           perp2D["time"] ,\
                           perp2D["zPos"] ,\
                           varName)
    p2DPerpPol.setPolData(pol2D["X"]     ,\
                          pol2D["Y"]     ,\
                          pol2D[varName] ,\
                          pol2D["time"]  ,\
                          pol2D["rhoPos"],\
                          varName)
    p2DPerpPol.plotAndSavePerpPolPlane()
#}}}

#{{{Driver2DFields
class Driver2DFields(DriverPlotFieldsSuperClass):
    """
    Class for driving of the plotting of the 2D fields.
    """
    #{{{constructor
    def __init__(self                     ,\
                 *args                    ,\
                 varName           = "n"  ,\
                 fluct             = False,\
                 varyMaxMin        = False,\
                 timeStampFolder   = True ,\
                 **kwargs):
        #{{{docstring
        """
        This constructor:
            * Calls the parent class
            * Updates the savePath and makes the folder
            * Set the member data

        Parameters
        ----------
        *args : str
            See parent class for details.
        varName : str
            Name of variable to collect and plot
        fluct : bool
            If fluctuations should be plotted instead of the normal
            fields.
        varyMaxMin : bool
            If the colorbar should be adjusted to the max and min of the
            current timestep.
            If False, the global max and min is used.
        timeStampFolder : bool
            Whether or not to timestamp the folder
        **kwargs : keyword arguments
            See parent class for details.
        """
        #}}}

        # Call the constructor of the parent class
        super().__init__(*args, **kwargs)

        # Update the savePath
        firstPathPart = os.path.dirname(self._savePath)
        if timeStampFolder:
            timePath = os.path.basename(self._savePath)
        else:
            timePath = ""
        self._savePath = os.path.join(firstPathPart, "field2D", timePath)

        # Make dir if not exists
        if not os.path.exists(self._savePath):
            os.makedirs(self._savePath)

        # Set the member data
        self._varName    = varName
        self._fluct      = fluct
        self._varyMaxMin = varyMaxMin
    #}}}

    #{{{setXIndToMaxGradInN
    def setXIndToMaxGradInN(self, steadyStatePath):
        #{{{docstring
        """
        Sets the xInd to the maximum radial gradient in n

        Parameters
        ----------
        steadyStatePath : str
            The steady state path to find the maximum gradient from.
        """
        #}}}

        self._xInd = findLargestRadialGradN(steadyStatePath)
    #}}}

    #{{{driver2DFieldsPerp
    def driver2DFieldsPerp(self):
        #{{{docstring
        """
        Wrapper to driver2DFieldsPerp
        """
        #}}}
        args =  (\
                 self._collectPaths     ,\
                 self._savePath         ,\
                 self._varName          ,\
                 self._convertToPhysical,\
                 self._xSlice           ,\
                 self._yInd             ,\
                 self._zSlice           ,\
                 self._tSlice           ,\
                 self._fluct            ,\
                 self._varyMaxMin       ,\
                )
        if self._useSubProcess:
            processes = Process(target = driver2DFieldPerpSingle, args = args)
            processes.start()
        else:
            driver2DFieldPerpSingle(*args)
    #}}}

    #{{{driver2DFieldsPar
    def driver2DFieldsPar(self):
        #{{{docstring
        """
        Wrapper to driver2DFieldsPar
        """
        #}}}
        args =  (\
                 self._collectPaths     ,\
                 self._savePath         ,\
                 self._varName          ,\
                 self._convertToPhysical,\
                 self._xSlice           ,\
                 self._ySlice           ,\
                 self._zInd             ,\
                 self._tSlice           ,\
                 self._fluct            ,\
                 self._varyMaxMin       ,\
                )
        if self._useSubProcess:
            processes = Process(target = driver2DFieldParSingle, args = args)
            processes.start()
        else:
            driver2DFieldParSingle(*args)
    #}}}

    #{{{driver2DFieldsPol
    def driver2DFieldsPol(self):
        #{{{docstring
        """
        Wrapper to driver2DFieldsPol
        """
        #}}}
        args =  (\
                 self._collectPaths     ,\
                 self._savePath         ,\
                 self._varName          ,\
                 self._convertToPhysical,\
                 self._xInd             ,\
                 self._ySlice           ,\
                 self._zSlice           ,\
                 self._tSlice           ,\
                 self._fluct            ,\
                 self._varyMaxMin       ,\
                )
        if self._useSubProcess:
            processes = Process(target = driver2DFieldPolSingle, args = args)
            processes.start()
        else:
            driver2DFieldPolSingle(*args)
    #}}}

    #{{{driver2DFieldsPerpPar
    def driver2DFieldsPerpPar(self):
        #{{{docstring
        """
        Wrapper to driver2DFieldsPerpPar
        """
        #}}}
        args =  (\
                 self._collectPaths     ,\
                 self._savePath         ,\
                 self._varName          ,\
                 self._convertToPhysical,\
                 self._xSlice           ,\
                 self._ySlice           ,\
                 self._zSlice           ,\
                 self._yInd             ,\
                 self._zInd             ,\
                 self._tSlice           ,\
                 self._fluct            ,\
                 self._varyMaxMin       ,\
                )
        if self._useSubProcess:
            processes = Process(target = driver2DFieldPerpParSingle, args = args)
            processes.start()
        else:
            driver2DFieldPerpParSingle(*args)
    #}}}

    #{{{driver2DFieldsPerpPol
    def driver2DFieldsPerpPol(self):
        #{{{docstring
        """
        Wrapper to driver2DFieldsPerpPol
        """
        #}}}
        args =  (\
                 self._collectPaths     ,\
                 self._savePath         ,\
                 self._varName          ,\
                 self._convertToPhysical,\
                 self._xSlice           ,\
                 self._ySlice           ,\
                 self._zSlice           ,\
                 self._xInd             ,\
                 self._yInd             ,\
                 self._tSlice           ,\
                 self._fluct            ,\
                 self._varyMaxMin       ,\
                )
        if self._useSubProcess:
            processes = Process(target = driver2DFieldPerpPolSingle, args = args)
            processes.start()
        else:
            driver2DFieldPerpPolSingle(*args)
    #}}}
#}}}

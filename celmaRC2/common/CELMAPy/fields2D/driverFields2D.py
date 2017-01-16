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

#{{{getVmaxVminLevels
def getVmaxVminLevels(plotSuperKwargs, tupleOfArrays, fluct, varyMaxMin):
    #{{{docstring
    """
    Getter for vmax, vmin, levels

    Parameters
    ----------
    plotSuperKwargs : dict
        Keyword arguments for the plot super class
    tupleOfArrays : tuple
        Tuple of the arrays to find the max and min of.

    Returns
    -------
    vmax : tuple
        Tuple of the vmax values used in contourf
    vmin : tuple
        Tuple of the vmin values used in contourf
    levels : tuple
        Tuple of the levels values used in contourf
    fluct : bool
        If the variables are fluctuating or not (will symmetrize the max
        and the min)
        Only active if plotSuperKwargs["vmax"] is not set
    varyMaxMin : bool
        If the max and min allowed to vary
        Only active if plotSuperKwargs["vmax"] is not set
    """
    #}}}
    if "vmax" in plotSuperKwargs.keys():
        if plotSuperKwargs["vmax"] is not None:
            vmax = plotSuperKwargs.pop("vmax")
            vmin = plotSuperKwargs.pop("vmin")
        else:
            plotSuperKwargs.pop("vmax")
            plotSuperKwargs.pop("vmin")
            vmax, vmin = getMaxMinAnimation(tupleOfArrays,\
                                            fluct        ,\
                                            varyMaxMin)
    else:
        vmax, vmin = getMaxMinAnimation(tupleOfArrays,\
                                        fluct        ,\
                                        varyMaxMin)
    if "levels" in plotSuperKwargs.keys():
        if plotSuperKwargs["levels"] is not None:
            levels = plotSuperKwargs.pop("levels")
        else:
            plotSuperKwargs.pop("levels")
            levels = getLevelsAnimation(vmax, vmin, 100)
    else:
        levels = getLevelsAnimation(vmax, vmin, 100)

    return vmax, vmin, levels
#}}}

#{{{driver2DFieldPerpSingle
def driver2DFieldPerpSingle(collectPaths     ,\
                            varName          ,\
                            convertToPhysical,\
                            xSlice           ,\
                            yInd             ,\
                            zSlice           ,\
                            tSlice           ,\
                            fluct            ,\
                            varyMaxMin       ,\
                            plotSuperKwargs  ,\
                                            ):
    #{{{docstring
    """
    Driver for plotting a single perpendicular 2D plot

    Parameters
    ----------
    collectPaths : tuple
        Paths to collect from.
        The corresponind 't_array' of the paths must be in ascending order.
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
    plotSuperKwargs : dict
        Keyword arguments for the plot super class
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
    # Set name
    ccf2D.setVarName(varName)
    # Execute the collection
    perp2D = ccf2D.executeCollectAndCalc()

    # Set the plot limits
    tupleOfArrays = (perp2D[varName],)
    vmax, vmin, levels = getVmaxVminLevels(plotSuperKwargs, tupleOfArrays, fluct, varyMaxMin)

    # Create the plotting object
    p2DPerp = PlotAnim2DPerp(ccf2D.uc        ,\
                             fluct    = fluct,\
                             **plotSuperKwargs)
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
                           varName          ,\
                           convertToPhysical,\
                           xSlice           ,\
                           ySlice           ,\
                           zInd             ,\
                           tSlice           ,\
                           fluct            ,\
                           varyMaxMin       ,\
                           plotSuperKwargs  ,\
                                            ):
    #{{{docstring
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
    plotSuperKwargs : dict
        Keyword arguments for the plot super class
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
    # Set name and
    ccf2D.setVarName(varName)
    # Execute the collection
    par2D = ccf2D.executeCollectAndCalc()

    # Set the plot limits
    tupleOfArrays = (par2D[varName],)
    vmax, vmin, levels = getVmaxVminLevels(plotSuperKwargs, tupleOfArrays, fluct, varyMaxMin)

    # Create the plotting object
    p2DPar = PlotAnim2DPar(ccf2D.uc        ,\
                           fluct    = fluct,\
                           **plotSuperKwargs)
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
                           varName          ,\
                           convertToPhysical,\
                           xInd             ,\
                           ySlice           ,\
                           zSlice           ,\
                           tSlice           ,\
                           fluct            ,\
                           varyMaxMin       ,\
                           plotSuperKwargs  ,\
                                            ):
    #{{{docstring
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
    plotSuperKwargs : dict
        Keyword arguments for the plot super class
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
    # Set name
    ccf2D.setVarName(varName)
    # Execute the collection
    pol2D = ccf2D.executeCollectAndCalc()

    # Set the plot limits
    tupleOfArrays = (pol2D[varName],)
    vmax, vmin, levels = getVmaxVminLevels(plotSuperKwargs, tupleOfArrays, fluct, varyMaxMin)

    # Create the plotting object
    p2DPol = PlotAnim2DPol(ccf2D.uc        ,\
                           fluct    = fluct,\
                           **plotSuperKwargs)

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
                               plotSuperKwargs  ,\
                                                ):
    #{{{docstring
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
    plotSuperKwargs : dict
        Keyword arguments for the plot super class
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
    tupleOfArrays = (perp2D[varName], par2D[varName])
    vmax, vmin, levels = getVmaxVminLevels(plotSuperKwargs, tupleOfArrays, fluct, varyMaxMin)

    # Create the plotting object
    p2DPerpPar = PlotAnim2DPerpPar(ccf2D.uc        ,\
                                   fluct    = fluct,\
                                   **plotSuperKwargs)
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
                               plotSuperKwargs  ,\
                                                ):
    #{{{docstring
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
    savePlot : bool
        Whether or no the plot should be saved.
    plotSuperKwargs : dict
        Keyword arguments for the plot super class
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
    tupleOfArrays = (perp2D[varName], pol2D[varName])
    vmax, vmin, levels = getVmaxVminLevels(plotSuperKwargs, tupleOfArrays, fluct, varyMaxMin)

    # Create the plotting object
    p2DPerpPol = PlotAnim2DPerpPol(ccf2D.uc        ,\
                                   fluct    = fluct,\
                                   **plotSuperKwargs)
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
    def __init__(self                       ,\
                 dmp_folders                ,\
                 plotSuperKwargs            ,\
                 guardSlicesAndIndicesKwargs,\
                 varName           = "n"    ,\
                 fluct             = False  ,\
                 varyMaxMin        = False  ,\
                 **kwargs):
        #{{{docstring
        """
        This constructor:
            * Calls the parent class
            * Set the member data
            * Updates the plotSuperKwargs

        Parameters
        ----------
        dmp_folders : tuple
            Tuple of the dmp_folder (output from bout_runners).
        plotSuperKwargs : dict
            Keyword arguments for the plot super class.
        guardSlicesAndIndicesKwargs : dict
            Keyword arguments for the slices, guards and indices.
        varName : str
            Name of variable to collect and plot
        fluct : bool
            If fluctuations should be plotted instead of the normal
            fields.
        varyMaxMin : bool
            If the colorbar should be adjusted to the max and min of the
            current timestep.
            If False, the global max and min is used.
        **kwargs : keyword arguments
            See parent class for details.
        """
        #}}}

        # Call the constructor of the parent class
        super().__init__(dmp_folders, **kwargs)

        # Set the member data
        self._varName    = varName
        self._fluct      = fluct
        self._varyMaxMin = varyMaxMin

        # Update the plotSuperKwargs dict
        plotSuperKwargs.update({"plotType":"field2D"})
        plotSuperKwargs.update({"dmp_folders":dmp_folders})
        self._plotSuperKwargs = plotSuperKwargs

        # Set the guards, slices and indices
        self.setGuardSlicesAndIndices(**guardSlicesAndIndicesKwargs)
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
                 self._collectPaths    ,\
                 self._varName         ,\
                 self.convertToPhysical,\
                 self._xSlice          ,\
                 self._yInd            ,\
                 self._zSlice          ,\
                 self._tSlice          ,\
                 self._fluct           ,\
                 self._varyMaxMin      ,\
                 self._plotSuperKwargs ,\
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
                 self._collectPaths    ,\
                 self._varName         ,\
                 self.convertToPhysical,\
                 self._xSlice          ,\
                 self._ySlice          ,\
                 self._zInd            ,\
                 self._tSlice          ,\
                 self._fluct           ,\
                 self._varyMaxMin      ,\
                 self._plotSuperKwargs ,\
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
                 self._collectPaths    ,\
                 self._varName         ,\
                 self.convertToPhysical,\
                 self._xInd            ,\
                 self._ySlice          ,\
                 self._zSlice          ,\
                 self._tSlice          ,\
                 self._fluct           ,\
                 self._varyMaxMin      ,\
                 self._plotSuperKwargs ,\
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
                 self._collectPaths    ,\
                 self._varName         ,\
                 self.convertToPhysical,\
                 self._xSlice          ,\
                 self._ySlice          ,\
                 self._zSlice          ,\
                 self._yInd            ,\
                 self._zInd            ,\
                 self._tSlice          ,\
                 self._fluct           ,\
                 self._varyMaxMin      ,\
                 self._plotSuperKwargs ,\
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
                 self._collectPaths    ,\
                 self._varName         ,\
                 self.convertToPhysical,\
                 self._xSlice          ,\
                 self._ySlice          ,\
                 self._zSlice          ,\
                 self._xInd            ,\
                 self._yInd            ,\
                 self._tSlice          ,\
                 self._fluct           ,\
                 self._varyMaxMin      ,\
                 self._plotSuperKwargs ,\
                )
        if self._useSubProcess:
            processes = Process(target = driver2DFieldPerpPolSingle, args = args)
            processes.start()
        else:
            driver2DFieldPerpPolSingle(*args)
    #}}}
#}}}

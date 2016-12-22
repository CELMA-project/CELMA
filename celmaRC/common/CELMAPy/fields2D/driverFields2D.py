#!/usr/bin/env python

"""
Contains single driver and driver class for 2D fields
"""

from ..plotHelpers import getMaxMinAnimation, getLevelsAnimation
from .collectAndCalcFields2D import CollectAndCalcFields2D
from .plotFields2D import (PlotAnim2DPerp,\
                           PlotAnim2DPar,\
                           PlotAnim2DPol,\
                           PlotAnim2DPerpPar,\
                           PlotAnim2DPerpPol,\
                           )

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
def driver2DFieldPerpParSingle(collectPaths ,\
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

    mode = "perpPar"

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


#    mode : ["perp"|"par"|"pol"|"perpPar"|"perpPol"]
#        * "perp"    - Only perpendicular plane will be plotted
#        * "par"     - Only parallel plane will be plotted
#        * "pol"     - Only poloidal plane will be plotted
#        * "perpPar" - The prependicular and parallel plan will be plotted
#        * "perpPol" - The prependicular and poloidal plan will be plotted

#    xSlice : [int|slice]
#        How the data will be sliced in x.
#        If "mode" is "pol" this must be an int.
#    ySlice : [int|slice]
#        How the data will be sliced in y.
#        If "mode" is "perp" this must be an int.
#    zSlice : [int|slice]
#        How the data will be sliced in z.
#        If "mode" is "pol" this must be an int.


#           #CollectAndCalcFields2D
#           fluct  = False ,\
#           mode   = "perp",\
#           #CollectAndCalcFieldsSuperClass
#           collectPaths              ,\
#           convertToPhysical = True  ,\
#           xguards           = False ,\
#           yguards           = False ,\
#           uc                = None  ,\
#           dh                = None  ,\
#
#           #PlotAnim2DPerp
#           pltSize
#           #PlotAnim2DSuperClass
#           fluct = None
#           #PlotAnimSuperClass
#           collectPaths     ,\
#           savePath         ,\
#           convertToPhysical,\
#           show = False     ,\
#           save = True      ,\










#
##{{{getCollectFieldsAndPlotOrder
#def getCollectFieldsAndPlotOrder(fieldPlotType, hyperIncluded):
#    #{{{doctring
#    """
#    Returns a tuple of the fields to collect and how to organize them
#
#    Parameters
#    ----------
#    fieldPlotType : str
#        What predefined fieldPlotType to plot.
#    hyperIncluded : bool
#        If hyper viscosities are used.
#
#    Returns
#    -------
#    collectFields : tuple
#        The fields to collect
#    plotOrder : tuple
#        The plot order
#    """
#    #}}}
#
#    plotOrder = None
#    if fieldPlotType == "mainFields":
#        collectFields  = ("lnN"       ,\
#                          "jPar"      ,\
#                          "phi"       ,\
#                          "vort"      ,\
#                          "vortD"     ,\
#                          "momDensPar",\
#                          "S"         ,\
#                         )
#        plotOrder = ("lnN"  , "phi"       ,\
#                     "n"    , "vortD"     ,\
#                     "jPar" , "vort"      ,\
#                     "uIPar", "momDensPar",\
#                     "uEPar", "S"         ,\
#                    )
#    elif fieldPlotType == "lnN":
#        collectFields  = ("ddt(lnN)"      ,\
#                          "lnNAdv"        ,\
#                          "gradUEPar"     ,\
#                          "lnNUeAdv"      ,\
#                          "srcN"          ,\
#                          "lnNRes"        ,\
#                          "lnNPerpArtVisc",\
#                          "lnNParArtVisc" ,\
#                         )
#    elif fieldPlotType == "jPar":
#        collectFields  = ("ddt(jPar)"      ,\
#                          "jParAdv"        ,\
#                          "uIParAdvSum"    ,\
#                          "uEParDoubleAdv" ,\
#                          "jParRes"        ,\
#                          "gradPhiLnN"     ,\
#                          "neutralERes"    ,\
#                          "neutralIRes"    ,\
#                          "jParPerpArtVisc",\
#                          "jParParArtVisc" ,\
#                         )
#    elif fieldPlotType == "momDensPar":
#        collectFields  = ("ddt(momDensPar)"   ,\
#                          "momDensAdv"        ,\
#                          "uIParAdvSum"       ,\
#                          "elPressure"        ,\
#                          "neutralIRes"       ,\
#                          "densDiffusion"     ,\
#                          "momDensPerpArtVisc",\
#                          "momDensParArtVisc" ,\
#                         )
#    elif fieldPlotType == "vortD":
#        collectFields  = ["ddt(vortD)"                ,\
#                          "vortNeutral"               ,\
#                          "potNeutral"                ,\
#                          "parDerDivUIParNGradPerpPhi",\
#                          "vortDAdv"                  ,\
#                          "kinEnAdvN"                 ,\
#                          "divParCur"                 ,\
#                          "vortDParArtVisc"           ,\
#                          "vortDPerpArtVisc"          ,\
#                         ]
#        if hyperIncluded:
#            collectFields.append("vortDHyperVisc")
#        collectFields = tuple(collectFields)
#    elif fieldPlotType == "vort":
#        collectFields  = ["ddt(vort)"               ,\
#                          "vortNeutral"             ,\
#                          "potNeutral"              ,\
#                          "DDYGradPerpPhiGradPerpUI",\
#                          "vortAdv"                 ,\
#                          "vortParAdv"              ,\
#                          "divParCur"               ,\
#                          "divSourcePhi"            ,\
#                          "vortParArtVisc"          ,\
#                          "vortPerpArtVisc"         ,\
#                        ]
#        if hyperIncluded:
#            collectFields.append("vortHyperVisc")
#        collectFields = tuple(collectFields)
#    else:
#        raise NotImplementedError("{} is not implemented".format(fieldPlotType))
#
#    return collectFields, plotOrder
##}}}
#
##{{{Driver1DFields
#class Driver1DFields(DriverPlotFieldsSuperClass):
#    """
#    Class for plotting of the 1D fields.
#    """
#
#    #{{{static members
#    _fieldPlotTypes = (\
#                       "mainFields",\
#                       "momDensPar",\
#                       "lnN"       ,\
#                       "jPar"      ,\
#                       "momDensPar",\
#                      )
#    #}}}
#
#    #{{{constructor
#    def __init__(self                   ,\
#                 *args                  ,\
#                 timeStampFolder = True ,\
#                 boussinesq      = False,\
#                 hyperIncluded   = False,\
#                 **kwargs):
#        #{{{docstring
#        """
#        This constructor:
#            * Calls the parent class
#            * Updates the fieldPlotTypes
#            * Updates the savePath and makes the folder
#            * Sets the member data
#
#        Parameters
#        ----------
#        *args : str
#            See parent class for details.
#        timeStampFolder : bool
#            Whether or not to timestamp the folder
#        boussinesq : bool
#            Whether or not the boussinesq approximation is used
#        hyperIncluded : bool
#            If hyper viscosities are used
#        **kwargs : keyword arguments
#            See parent class for details.
#        """
#        #}}}
#
#        # Call the constructor of the parent class
#        super().__init__(*args, **kwargs)
#
#        # Update fieldPlotTypes
#        if not(boussinesq):
#            Driver1DFields._fieldPlotTypes =\
#                    list(Driver1DFields._fieldPlotTypes)
#            Driver1DFields._fieldPlotTypes.append("vortD")
#        else:
#            Driver1DFields._fieldPlotTypes =\
#                    list(Driver1DFields._fieldPlotTypes)
#            Driver1DFields._fieldPlotTypes.append("vort")
#
#        # Recast to tuple
#        Driver1DFields._fieldPlotTypes =\
#                tuple(set(Driver1DFields._fieldPlotTypes))
#
#        # Update the savePath
#        firstPathPart = os.path.dirname(self._savePath)
#        if timeStampFolder:
#            timePath = os.path.basename(self._savePath)
#        else:
#            timePath = ""
#        self._savePath = os.path.join(firstPathPart, "field1D", timePath)
#
#        # Make dir if not exists
#        if not os.path.exists(self._savePath):
#            os.makedirs(self._savePath)
#
#        # Set member data
#        self._hyperIncluded = hyperIncluded
#    #}}}
#
#    #{{{driver1DFieldsAll
#    def driver1DFieldsAll(self):
#        #{{{docstring
#        """
#        Wrapper to driver1DFieldSingle.
#
#        Drives all implemeted combinations of driver1DFieldSingle using the
#        member data.
#        """
#        #}}}
#        if self._useSubProcess:
#            parallelProcess = Process(\
#                                 target = self.driver1DFieldsParallel,\
#                                 args   = ()                         ,\
#                                 kwargs = {}                         ,\
#                                )
#
#            radialProcess = Process(\
#                                 target = self.driver1DFieldsRadial,\
#                                 args   = ()                       ,\
#                                 kwargs = {}                       ,\
#                                )
#            parallelProcess.start()
#            radialProcess  .start()
#            parallelProcess.join()
#            radialProcess  .join()
#        else:
#            self.driver1DFieldsParallel()
#            self.driver1DFieldsRadial()
#    #}}}
#
#    #{{{driver1DFieldsParallel
#    def driver1DFieldsParallel(self):
#        #{{{docstring
#        """
#        Wrapper to driver1DFieldSingle.
#
#        Drives all implemeted combinations of driver1DFieldSingle in
#        "parallel" mode using the member data.
#        """
#        #}}}
#        argTemplate =  [\
#                        self._collectPaths     ,\
#                        self._savePath         ,\
#                        self._convertToPhysical,\
#                        self._xInd             ,\
#                        self._ySlice           ,\
#                        self._zInd             ,\
#                        self._tSlice           ,\
#                        "parallel"             ,\
#                        self._hyperIncluded    ,\
#                       ]
#        if self._useSubProcess:
#            processes = {}
#            for fieldPlotType in Driver1DFields._fieldPlotTypes:
#                args = argTemplate.copy()
#                args.insert(2, fieldPlotType)
#                processes[fieldPlotType] =\
#                    Process(target = driver1DFieldSingle,\
#                            args = args                 ,\
#                            )
#                processes[fieldPlotType].start()
#            for fieldPlotType in Driver1DFields._fieldPlotTypes:
#                processes[fieldPlotType].join()
#        else:
#            for fieldPlotType in Driver1DFields._fieldPlotTypes:
#                args = argTemplate.copy()
#                args.insert(2, fieldPlotType)
#                driver1DFieldSingle(*args)
#    #}}}
#
#    #{{{driver1DFieldsRadial
#    def driver1DFieldsRadial(self):
#        #{{{docstring
#        """
#        Wrapper to driver1DFieldSingle.
#
#        Drives all implemeted combinations of driver1DFieldSingle in
#        "radial" mode using the member data.
#        """
#        #}}}
#        argTemplate =  [
#                        self._collectPaths     ,\
#                        self._savePath         ,\
#                        self._convertToPhysical,\
#                        self._xSlice           ,\
#                        self._yInd             ,\
#                        self._zInd             ,\
#                        self._tSlice           ,\
#                        "radial"               ,\
#                        self._hyperIncluded    ,\
#                       ]
#        if self._useSubProcess:
#            processes = {}
#            for fieldPlotType in Driver1DFields._fieldPlotTypes:
#                args = argTemplate.copy()
#                args.insert(2, fieldPlotType)
#                processes[fieldPlotType] =\
#                    Process(target = driver1DFieldSingle,\
#                            args = args                 ,\
#                            )
#                processes[fieldPlotType].start()
#            for fieldPlotType in Driver1DFields._fieldPlotTypes:
#                processes[fieldPlotType].join()
#        else:
#            for fieldPlotType in Driver1DFields._fieldPlotTypes:
#                args = argTemplate.copy()
#                args.insert(2, fieldPlotType)
#                driver1DFieldSingle(*args)
#    #}}}
##}}}

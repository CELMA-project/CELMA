#!/usr/bin/env python

"""
Contains single driver and driver class for 1D fields
"""

from .collectAndCalcFields1D import CollectAndCalcFields1D
from .plotFields1D import PlotAnim1DRadial, PlotAnim1DParallel
from ..collectAndCalcHelpers import calcN, calcUIPar, calcUEPar
from ..superClasses import DriverPlotFieldsSuperClass
from multiprocessing import Process

#{{{driver1DFieldSingle
def driver1DFieldSingle(collectPaths     ,\
                        fieldPlotType    ,\
                        convertToPhysical,\
                        xSlice           ,\
                        ySlice           ,\
                        zSlice           ,\
                        tSlice           ,\
                        mode             ,\
                        hyperIncluded    ,\
                        plotSuperKwargs  ,\
                        ):
    #{{{docstring
    """
    Driver for plotting a single predefined fieldPlotType plot

    Parameters
    ----------
    collectPaths : tuple
        Paths to collect from.
        The corresponind 't_array' of the paths must be in ascending order.
    fieldPlotType : str
        What predefined fieldPlotType to plot.
        See getCollectFieldsAndPlotOrder for details.
    convertToPhysical : bool
        Whether or not to convert to physical units.
    xSlice : [int|slice]
        How the data will be sliced in x.
        If "mode" is "parallel" this must be an int.
    ySlice : [int|slice]
        How the data will be sliced in y.
        If "mode" is "radial" this must be an int.
    zSlice : int
        How the data will be sliced in z.
    tSlice : [None|slice]
        How the data will be sliced in t.
    mode : ["radial"|"parallel"]
        * "radial"    - Radial profiles will be used.
        * "parallel"  - Parallel profiles will be used.
    hyperIncluded : bool
        If hyper viscosities are used.
    plotSuperKwargs : dict
        Keyword arguments for the plot super class
    """
    #}}}

    # Magic number
    processing = None

    if mode == "radial":
        PlotClass = PlotAnim1DRadial
    elif mode == "parallel":
        PlotClass = PlotAnim1DParallel

    collectFields, plotOrder =\
            getCollectFieldsAndPlotOrder(fieldPlotType, hyperIncluded)

    if fieldPlotType != "mainFields" and convertToPhysical:
        # NOTE: Normalization for each term not implemented
        convertToPhysical = False
        print("fieldPlotType is not 'mainFields', "\
              "setting 'convertToPhysical' to False")

    ccf1D = CollectAndCalcFields1D(collectPaths,\
                                   mode = mode,\
                                   processing = processing,\
                                   convertToPhysical = convertToPhysical)

    ccf1D.setSlice(xSlice,\
                   ySlice,\
                   zSlice,\
                   tSlice)

    dict1D = {}

    for field in collectFields:
        ccf1D.setVarName(field)
        dict1D.update(ccf1D.executeCollectAndCalc())

    if fieldPlotType == "mainFields":
        # Non-collects
        dict1D.update({"n"    : calcN(dict1D["lnN"],\
                                      not(ccf1D.convertToPhysical),\
                                      ccf1D.uc)})
        dict1D.update({"uIPar": calcUIPar(dict1D["momDensPar"], dict1D["n"])})
        dict1D.update({"uEPar": calcUEPar(dict1D["uIPar"],\
                                          dict1D["jPar"],\
                                          dict1D["n"],\
                                          not(ccf1D.convertToPhysical))})

    p1D = PlotClass(ccf1D.uc        ,\
                    **plotSuperKwargs)

    p1D.setData(dict1D, fieldPlotType, plotOrder=plotOrder)
    p1D.plotAndSaveProfile()
#}}}

#{{{getCollectFieldsAndPlotOrder
def getCollectFieldsAndPlotOrder(fieldPlotType, hyperIncluded):
    #{{{docstring
    """
    Returns a tuple of the fields to collect and how to organize them

    Parameters
    ----------
    fieldPlotType : str
        What predefined fieldPlotType to plot.
    hyperIncluded : bool
        If hyper viscosities are used.

    Returns
    -------
    collectFields : tuple
        The fields to collect
    plotOrder : tuple
        The plot order
    """
    #}}}

    plotOrder = None
    if fieldPlotType == "mainFields":
        collectFields  = ("lnN"       ,\
                          "jPar"      ,\
                          "phi"       ,\
                          "vort"      ,\
                          "vortD"     ,\
                          "momDensPar",\
                          "S"         ,\
                         )
        plotOrder = ("lnN"  , "phi"       ,\
                     "n"    , "vortD"     ,\
                     "jPar" , "vort"      ,\
                     "uIPar", "momDensPar",\
                     "uEPar", "S"         ,\
                    )
    if fieldPlotType == "mainFieldsBoussinesq":
        collectFields  = ("lnN"       ,\
                          "jPar"      ,\
                          "phi"       ,\
                          "vort"      ,\
                          "momDensPar",\
                          "S"         ,\
                         )
        plotOrder = ("lnN"  , "phi"       ,\
                     "n"    , "vort"      ,\
                     "jPar" , "momDensPar",\
                     "uIPar", "S"         ,\
                     "uEPar",\
                    )
    elif fieldPlotType == "lnN":
        collectFields  = ("ddt(lnN)"      ,\
                          "lnNAdv"        ,\
                          "gradUEPar"     ,\
                          "lnNUeAdv"      ,\
                          "srcN"          ,\
                          "lnNRes"        ,\
                          "lnNPerpArtVisc",\
                          "lnNParArtVisc" ,\
                         )
    elif fieldPlotType == "jPar":
        collectFields  = ("ddt(jPar)"      ,\
                          "jParAdv"        ,\
                          "uEParAdv"       ,\
                          "uIParAdv"       ,\
                          "jParParAdv"     ,\
                          "gradPhiLnN"     ,\
                          "jParRes"        ,\
                          "neutralERes"    ,\
                          "neutralIRes"    ,\
                          "jParPerpArtVisc",\
                          "jParParArtVisc" ,\
                         )
    elif fieldPlotType == "momDensPar":
        collectFields  = ("ddt(momDensPar)"   ,\
                          "momDensAdv"        ,\
                          "uIFluxAdv"         ,\
                          "uIParAdv"          ,\
                          "elPressure"        ,\
                          "densDiffusion"     ,\
                          "neutralIRes"       ,\
                          "neutralEResMu"     ,\
                          "momDensPerpArtVisc",\
                          "momDensParArtVisc" ,\
                         )
    elif fieldPlotType == "vortD":
        collectFields  = ["ddt(vortD)"                ,\
                          "vortNeutral"               ,\
                          "potNeutral"                ,\
                          "parDerDivUIParNGradPerpPhi",\
                          "vortDAdv"                  ,\
                          "kinEnAdvN"                 ,\
                          "divParCur"                 ,\
                          "vortDParArtVisc"           ,\
                          "vortDPerpArtVisc"          ,\
                         ]
        if hyperIncluded:
            collectFields.append("vortDHyperVisc")
        collectFields = tuple(collectFields)
    elif fieldPlotType == "vort":
        collectFields  = ["ddt(vort)"               ,\
                          "vortNeutral"             ,\
                          "potNeutral"              ,\
                          "DDYGradPerpPhiGradPerpUI",\
                          "vortAdv"                 ,\
                          "vortParAdv"              ,\
                          "divParCur"               ,\
                          "divSourcePhi"            ,\
                          "vortParArtVisc"          ,\
                          "vortPerpArtVisc"         ,\
                        ]
        if hyperIncluded:
            collectFields.append("vortHyperVisc")
        collectFields = tuple(collectFields)
    else:
        message = "{} is not implemented".format(fieldPlotType)
        raise NotImplementedError(message)

    return collectFields, plotOrder
#}}}

#{{{Driver1DFields
class Driver1DFields(DriverPlotFieldsSuperClass):
    """
    Class for the driving of the plotting of the 1D fields.
    """

    #{{{static members
    _fieldPlotTypes = (\
                       "momDensPar",\
                       "lnN"       ,\
                       "jPar"      ,\
                       "momDensPar",\
                      )
    #}}}

    #{{{constructor
    def __init__(self                       ,\
                 dmp_folders                ,\
                 plotSuperKwargs            ,\
                 guardSlicesAndIndicesKwargs,\
                 timeStampFolder = True     ,\
                 boussinesq      = False    ,\
                 hyperIncluded   = False    ,\
                 **kwargs):
        #{{{docstring
        """
        This constructor:
            * Calls the parent class
            * Updates the fieldPlotTypes
            * Sets the member data
            * Updates the plotSuperKwargs

        Parameters
        ----------
        dmp_folders : tuple
            Tuple of the dmp_folder (output from bout_runners).
        plotSuperKwargs : dict
            Keyword arguments for the plot super class.
        timeStampFolder : bool
            Whether or not to timestamp the folder
        boussinesq : bool
            Whether or not the boussinesq approximation is used
        hyperIncluded : bool
            If hyper viscosities are used
        **kwargs : keyword arguments
            See parent class for details.
        """
        #}}}

        # Call the constructor of the parent class
        super().__init__(dmp_folders, **kwargs)

        # Update fieldPlotTypes
        if not(boussinesq):
            Driver1DFields._fieldPlotTypes =\
                    list(Driver1DFields._fieldPlotTypes)
            Driver1DFields._fieldPlotTypes.append("vortD")
            Driver1DFields._fieldPlotTypes.append("mainFields")
        else:
            Driver1DFields._fieldPlotTypes =\
                    list(Driver1DFields._fieldPlotTypes)
            Driver1DFields._fieldPlotTypes.append("vort")
            Driver1DFields._fieldPlotTypes.append("mainFieldsBoussinesq")

        # Recast to tuple
        Driver1DFields._fieldPlotTypes =\
                tuple(set(Driver1DFields._fieldPlotTypes))

        # Set member data
        self._hyperIncluded = hyperIncluded

        # Update the plotSuperKwargs dict
        plotSuperKwargs.update({"plotType":"field1D"})
        plotSuperKwargs.update({"dmp_folders":dmp_folders})
        self._plotSuperKwargs = plotSuperKwargs

        # Set the guards, slices and indices
        self.setGuardSlicesAndIndices(**guardSlicesAndIndicesKwargs)
    #}}}

    #{{{driver1DFieldsAll
    def driver1DFieldsAll(self):
        #{{{docstring
        """
        Wrapper to driver1DFieldSingle.

        Drives all implemeted combinations of driver1DFieldSingle using the
        member data.
        """
        #}}}
        if self._useSubProcess:
            parallelProcess = Process(\
                                 target = self.driver1DFieldsParallel,\
                                 args   = ()                         ,\
                                 kwargs = {}                         ,\
                                )

            radialProcess = Process(\
                                 target = self.driver1DFieldsRadial,\
                                 args   = ()                       ,\
                                 kwargs = {}                       ,\
                                )
            parallelProcess.start()
            radialProcess  .start()
            parallelProcess.join()
            radialProcess  .join()
        else:
            self.driver1DFieldsParallel()
            self.driver1DFieldsRadial()
    #}}}

    #{{{driver1DFieldsParallel
    def driver1DFieldsParallel(self):
        #{{{docstring
        """
        Wrapper to driver1DFieldSingle.

        Drives all implemeted combinations of driver1DFieldSingle in
        "parallel" mode using the member data.
        """
        #}}}
        argTemplate =  [\
                        self._collectPaths    ,\
                        self.convertToPhysical,\
                        self._xInd            ,\
                        self._ySlice          ,\
                        self._zInd            ,\
                        self._tSlice          ,\
                        "parallel"            ,\
                        self._hyperIncluded   ,\
                        self._plotSuperKwargs ,\
                       ]
        if self._useSubProcess:
            processes = {}
            for fieldPlotType in Driver1DFields._fieldPlotTypes:
                args = argTemplate.copy()
                args.insert(1, fieldPlotType)
                processes[fieldPlotType] =\
                    Process(target = driver1DFieldSingle,\
                            args = args                 ,\
                            )
                processes[fieldPlotType].start()
            for fieldPlotType in Driver1DFields._fieldPlotTypes:
                processes[fieldPlotType].join()
        else:
            for fieldPlotType in Driver1DFields._fieldPlotTypes:
                args = argTemplate.copy()
                args.insert(1, fieldPlotType)
                driver1DFieldSingle(*args)
    #}}}

    #{{{driver1DFieldsRadial
    def driver1DFieldsRadial(self):
        #{{{docstring
        """
        Wrapper to driver1DFieldSingle.

        Drives all implemeted combinations of driver1DFieldSingle in
        "radial" mode using the member data.
        """
        #}}}
        argTemplate =  [
                        self._collectPaths    ,\
                        self.convertToPhysical,\
                        self._xSlice          ,\
                        self._yInd            ,\
                        self._zInd            ,\
                        self._tSlice          ,\
                        "radial"              ,\
                        self._hyperIncluded   ,\
                        self._plotSuperKwargs ,\
                       ]
        if self._useSubProcess:
            processes = {}
            for fieldPlotType in Driver1DFields._fieldPlotTypes:
                args = argTemplate.copy()
                args.insert(1, fieldPlotType)
                processes[fieldPlotType] =\
                    Process(target = driver1DFieldSingle,\
                            args = args                 ,\
                            )
                processes[fieldPlotType].start()
            for fieldPlotType in Driver1DFields._fieldPlotTypes:
                processes[fieldPlotType].join()
        else:
            for fieldPlotType in Driver1DFields._fieldPlotTypes:
                args = argTemplate.copy()
                args.insert(1, fieldPlotType)
                driver1DFieldSingle(*args)
    #}}}
#}}}

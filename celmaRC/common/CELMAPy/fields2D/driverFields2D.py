#!/usr/bin/env python

"""
Contains single driver and driver class for 1D fields
"""

from .collectAndCalcFields1D import CollectAndCalcFields1D
from .plotFields1D import PlotAnim1DRadial, PlotAnim1DParallel
from ..collectAndCalcHelpers import calcN, calcUIPar, calcUEPar
from ..superClasses import DriverPlotFieldsSuperClass
import os


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

#{{{driver1DFieldSingle
def driver1DFieldSingle(collectPaths,\
                        fieldPlotType,\
                        savePath,\
                        convertToPhysical,\
                        xSlice,\
                        ySlice,\
                        zSlice,\
                        tSlice,\
                        mode):
    #{{{doctring
    """
    Driver for plotting a single predefined fieldPlotType plot

    Parameters
    ----------
    collectPaths : tuple
        Paths to collect from.
        The corresponind 't_array' of the paths must be in ascending order.
    fieldPlotType : str
        What predefined fieldPlotType to plot.
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
        * "radial"    - Radial profiles will be used
        * "parallel"  - Parallel profiles will be used
    """
    #}}}

    # Magic number
    processing = None

    if mode == "radial":
        PlotClass = PlotAnim1DRadial
    elif mode == "parallel":
        PlotClass = PlotAnim1DParallel

    collectFields, plotOrder = getCollectFieldsAndPlotOrder(fieldPlotType)

    if fieldPlotType != "mainField" and convertToPhysical:
        # NOTE: Normalization for each term not implemented
        convertToPhysical = False
        print("fieldPlotType is not 'mainField', "\
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

    if fieldPlotType == "mainField":
        # Non-collects
        dict1D.update({"n"    : calcN(dict1D["lnN"],\
                                      not(ccf1D.convertToPhysical),\
                                      ccf1D.uc)})
        dict1D.update({"uIPar": calcUIPar(dict1D["momDensPar"], dict1D["n"])})
        dict1D.update({"uEPar": calcUEPar(dict1D["uIPar"],\
                                          dict1D["jPar"],\
                                          dict1D["n"],\
                                          not(ccf1D.convertToPhysical))})

    p1D = PlotClass(savePath, ccf1D.convertToPhysical)
    p1D.setRadialData(dict1D, fieldPlotType, savePath, plotOrder=plotOrder)
    p1D.plotAndSaveRadialProfile()
#}}}

#{{{getCollectFieldsAndPlotOrder
def getCollectFieldsAndPlotOrder(fieldPlotType):
    #{{{doctring
    """
    Returns a tuple of the fields to collect and how to organize them

    Parameters
    ----------
    fieldPlotType : str
        What predefined fieldPlotType to plot

    Returns
    -------
    collectFields : tuple
        The fields to collect
    plotOrder : tuple
        The plot order
    """
    #}}}

    if fieldPlotType == "mainField":
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

    return collectFields, plotOrder
#}}}

#{{{Driver1DFields
class Driver1DFields(DriverPlotFieldsSuperClass):
    """
    Class for plotting of the 1D fields
    """

    #{{{static members
    _fieldPlotTypes = (\
                   "mainField",\
                  )
    #}}}

    #{{{constructor
    def __init__(self,
                 *args      ,\
                 xInd = None,\
                 yInd = None,\
                 **kwargs):
        #{{{docstring
        """
        This constructor:
            * Calls the parent class
            * Updates the savePath
            * Sets the common member data

        Parameters
        ----------
        *args : str
            See parent class for details.
        xInd : int
            Fixed rho index for parallel plotting.
        yInd : int
            Fixed z index for radial plotting.
        **kwargs : keyword arguments
            See parent class for details.
        """
        #}}}

        # Call the constructor of the parent class
        super().__init__(*args, **kwargs)

        # Update the savePath
        self._savePath = os.path.join(self._savePath, "field1D")

        # Set the member data
        self._xInd = xInd
        self._yInd = yInd
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
        self.driver1DFieldParallel()
        self.driver1DFieldRadial()
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
        for fieldPlotType in Driver1DFields._fieldPlotTypes:
            driver1DFieldSingle(\
                        self._collectPaths     ,\
                        fieldPlotType          ,\
                        self._savePath         ,\
                        self._convertToPhysical,\
                        self._xInd             ,\
                        self._ySlice           ,\
                        self._zInd             ,\
                        self._tSlice           ,\
                        "parallel")
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
        for fieldPlotType in Driver1DFields._fieldPlotTypes:
            driver1DFieldSingle(\
                        self._collectPaths     ,\
                        fieldPlotType          ,\
                        self._savePath         ,\
                        self._convertToPhysical,\
                        self._xSlice           ,\
                        self._yInd             ,\
                        self._zInd             ,\
                        self._tSlice           ,\
                        "parallel")
    #}}}
#}}}

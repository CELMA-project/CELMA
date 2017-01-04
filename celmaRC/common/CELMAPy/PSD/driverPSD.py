#!/usr/bin/env python

"""
Contains single driver and driver class for the power spectral density
"""

from ..superClasses import DriverPointsSuperClass
from ..collectAndCalcHelpers import getSizes
from .collectAndCalcPSD import CollectAndCalcPSD
from .plotPSD import PlotPSD
import numpy as np
from multiprocessing import Process

#{{{driverPSD
def driverPSD(collectPaths     ,\
              varName          ,\
              convertToPhysical,\
              mode             ,\
              indicesArgs      ,\
              indicesKwargs    ,\
              plotSuperKwargs  ,\
             ):
    #{{{docstring
    """
    Driver for plotting power spectral density.

    Parameters
    ----------
    collectPaths : tuple
        Paths to collect from.
        The corresponind 't_array' of the paths must be in ascending order.
    varName : str
        The variable name which will be used.
    convertToPhysical : bool
        Whether or not to convert to physical units.
    mode : ["normal"|"fluct"]
        If mode is "normal" the raw data is given as an output.
        If mode is "fluct" the fluctuations are given as an output.
    indicesArgs : tuple
        Contains xInd, yInd and zInd.
        See CollectAndCalcPointsSuperClass.setIndices for details.
    indicesKwargs : dict
        Contains tslice, nPoints, equallySpace and steadyStatePath.
        See CollectAndCalcPointsSuperClass.setIndices for details.
    plotSuperKwargs : dict
        Keyword arguments for the plot super class.
    """
    #}}}

    # Create collect object
    ccPSD = CollectAndCalcPSD(collectPaths                         ,\
                              mode              = mode             ,\
                              convertToPhysical = convertToPhysical,\
                             )

    # Set the slice
    ccPSD.setIndices(*indicesArgs, **indicesKwargs)

    # Set name
    ccPSD.setVarName(varName)

    # Execute the collection
    tt = ccPSD.executeCollectAndCalc()
    tt = ccPSD.convertTo1D(tt)

    # Calculate the PSD
    PSD = ccPSD.calcPSD(tt)

    # Plot
    pPSD = PlotPSD(ccPSD.uc         ,\
                   **plotSuperKwargs)
    pPSD.setData(PSD, mode)
    pPSD.plotSaveShowPSD()
#}}}

#{{{driverPSD2D
def driverPSD2D(collectPaths     ,\
                varName          ,\
                convertToPhysical,\
                mode             ,\
                indicesArgs      ,\
                indicesKwargs    ,\
                plotSuperKwargs  ,\
               ):
    #{{{docstring
    """
    Driver for plotting power spectral density.

    NOTE: In order to be efficient, one should have collected with the
          field1D class rather than with the points collecter.
          In any case, this gets the job done

    Parameters
    ----------
    collectPaths : tuple
        Paths to collect from.
        The corresponind 't_array' of the paths must be in ascending order.
    varName : str
        The variable name which will be used.
    convertToPhysical : bool
        Whether or not to convert to physical units.
    mode : ["normal"|"fluct"]
        If mode is "normal" the raw data is given as an output.
        If mode is "fluct" the fluctuations are given as an output.
    indicesArgs : tuple
        Contains xInd, yInd and zInd.
        See CollectAndCalcPointsSuperClass.setIndices for details.
    indicesKwargs : dict
        Contains tslice, nPoints, equallySpace and steadyStatePath.
        See CollectAndCalcPointsSuperClass.setIndices for details.
    plotSuperKwargs : dict
        Keyword arguments for the plot super class.
    """
    #}}}

    # Create collect object
    ccPSD = CollectAndCalcPSD(collectPaths                         ,\
                              mode              = mode             ,\
                              convertToPhysical = convertToPhysical,\
                             )

    # Get all the xInds
    nx = getSizes(collectPaths[0], "x")
    xInds = tuple(range(nx))
    indicesArgs = list(indicesArgs)
    indicesArgs[0] = xInds
    indicesArgs = tuple(indicesArgs)

    # Set the slice
    ccPSD.setIndices(*indicesArgs, **indicesKwargs)

    # Set name
    ccPSD.setVarName(varName)

    # Execute the collection
    tt = ccPSD.executeCollectAndCalc()
    tt = ccPSD.convertTo1D(tt)

    # Calculate the PSD
    PSD = ccPSD.calcPSD(tt)

    # Recast to 2d np.array with dimensions (t, nx)
    keys = tuple(sorted(list(PSD.keys())))
    # Make the rho axis
    rho = np.zeros(len(keys))
    for nr in range(len(rho)):
        rho[nr] = float(keys[nr].split(",")[0])

    # Make the frequency axis
    freq = PSD[keys[0]]["{}PSDX".format(varName)]

    # Make the 2D PSD matrix
    PSDMat = np.zeros((len(PSD[keys[0]]["{}PSDY".format(varName)]) ,len(keys)))
    for xInd in range(len(rho)):
        PSDMat[:, xInd] = PSD[keys[xInd]]["{}PSDY".format(varName)]

    import pdb; pdb.set_trace()
    # Plot
    pPSD = PlotPSD(ccPSD.uc         ,\
                   **plotSuperKwargs)
    pPSD.setData(PSD, mode)
    pPSD.plotSaveShowPSD()
#}}}

#{{{DriverPSD
class DriverPSD(DriverPointsSuperClass):
    """
    Class for driving of the plotting of the power spectral density.
    """

    #{{{Constructor
    def __init__(self                       ,\
                 dmp_folders                ,\
                 indicesArgs                ,\
                 indicesKwargs              ,\
                 plotSuperKwargs            ,\
                 varName           = "n"    ,\
                 mode              = "fluct",\
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
        indicesArgs : tuple
            Contains xInd, yInd and zInd.
            See CollectAndCalcPointsSuperClass.setIndices for details.
        indicesKwargs : dict
            Contains tslice, nPoints, equallySpace and steadyStatePath.
            See CollectAndCalcPointsSuperClass.setIndices for details.
        plotSuperKwargs : dict
            Keyword arguments for the plot super class.
        varName : str
            Name of variable to collect and plot
        mode : ["normal"|"fluct"]
            If mode is "normal" the raw data is given as an output.
            If mode is "fluct" the fluctuations are given as an output.
        **kwargs : keyword arguments
            See parent class for details.
        """
        #}}}

        # Call the constructor of the parent class
        super().__init__(dmp_folders, **kwargs)

        # Set the member data
        self._varName       = varName
        self._mode          = mode
        self._indicesArgs   = indicesArgs
        self._indicesKwargs = indicesKwargs

        # Update the plotSuperKwargs dict
        plotSuperKwargs.update({"dmp_folders":dmp_folders})
        plotSuperKwargs.update({"plotType"   :"PSDs"})
        self._plotSuperKwargs = plotSuperKwargs
    #}}}

    #{{{driverPSD
    def driverPSD(self):
        #{{{docstring
        """
        Wrapper to driverPSD
        """
        #}}}
        args =  (\
                 self._collectPaths    ,\
                 self._varName         ,\
                 self.convertToPhysical,\
                 self._mode            ,\
                 self._indicesArgs     ,\
                 self._indicesKwargs   ,\
                 self._plotSuperKwargs ,\
                )
        if self._useSubProcess:
            processes = Process(target = driverPSD, args = args)
            processes.start()
        else:
            driverPSD(*args)
    #}}}
#}}}

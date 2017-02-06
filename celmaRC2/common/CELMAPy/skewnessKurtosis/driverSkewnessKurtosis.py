#!/usr/bin/env python

"""
Contains single driver and driver class for the power spectral density
"""

from ..superClasses import DriverPointsSuperClass
from ..collectAndCalcHelpers import getGridSizes
from .collectAndCalcSkewnessKurtosis import CollectAndCalcSkewnessKurtosis
from .plotSkewnessKurtosis import PlotSkewnessKurtosis
import numpy as np
from multiprocessing import Process

#{{{driverSkewnessKurtosis
def driverSkewnessKurtosis(collectPaths     ,\
                           varName          ,\
                           convertToPhysical,\
                           mode             ,\
                           indicesArgs      ,\
                           indicesKwargs    ,\
                           plotSuperKwargs  ,\
                          ):
    #{{{docstring
    """
    Driver for plotting the skewness and kurtosis.

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
    ccSK = CollectAndCalcSkewnessKurtosis(\
                collectPaths                         ,\
                mode              = mode             ,\
                convertToPhysical = convertToPhysical,\
               )

    # Get all the xInds
    nx = getGridSizes(collectPaths[0], "x")
    xInds = tuple(range(nx))
    indicesArgs = list(indicesArgs)
    indicesArgs[0] = xInds
    indicesArgs = tuple(indicesArgs)

    # Set the slice
    ccSK.setIndices(*indicesArgs, **indicesKwargs)

    # Set name
    ccSK.setVarName(varName)

    # Execute the collection
    tt = ccSK.executeCollectAndCalc()
    tt = ccSK.convertTo1D(tt)

    # Calculate the skewness and kurtosis
    skewKurt = ccSK.calcSkewnessKurtosis(tt)

    # Recast to position
    keys    = tuple(sorted(list(skewKurt.keys())))
    skewKey = "{}Skew".format(varName)
    kurtKey = "{}Kurt".format(varName)
    # Make the rho, skew and kurt axes
    rho   = np.zeros(len(keys))
    skews = np.zeros(len(keys))
    kurts = np.zeros(len(keys))
    for nr in range(len(rho)):
        rho  [nr] = float(keys[nr].split(",")[0])
        skews[nr] = skewKurt[keys[nr]][skewKey]
        kurts[nr] = skewKurt[keys[nr]][kurtKey]

    # Make a new skewKurt
    skewKurt = {\
                "skew"     : skews                ,\
                "kurt"     : kurts                ,\
                "rho"      : rho                  ,\
                "thetaPos" : keys[0].split(",")[1],\
                "zPos"     : keys[0].split(",")[2],\
                "varName"  : varName              ,\
                }

    # Plot
    pSK = PlotSkewnessKurtosis(ccSK.uc         ,\
                               **plotSuperKwargs)
    pSK.setData(skewKurt, mode)
    pSK.plotSaveShowSkewnessKurtosis()
#}}}

#{{{DriverSkewnessKurtosis
class DriverSkewnessKurtosis(DriverPointsSuperClass):
    """
    Class for driving of the plotting of the skewness and kurtosis.
    """

    #{{{Constructor
    def __init__(self                      ,\
                 dmp_folders               ,\
                 indicesArgs               ,\
                 indicesKwargs             ,\
                 plotSuperKwargs           ,\
                 varName          = "n"    ,\
                 mode             = "fluct",\
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
        plotSuperKwargs.update({"plotType"   :"skewKurt"})
        self._plotSuperKwargs = plotSuperKwargs
    #}}}

    #{{{driverSkewnessKurtosis
    def driverSkewnessKurtosis(self):
        #{{{docstring
        """
        Wrapper to driverSkewnessKurtosis2D
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
        if self._useMultiProcess:
            processes = Process(target = driverSkewnessKurtosis, args = args)
            processes.start()
        else:
            driverSkewnessKurtosis(*args)
    #}}}
#}}}

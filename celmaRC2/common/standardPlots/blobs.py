#!/usr/bin/env python

"""Post-processor for blobs"""

import os, sys
# If we add to sys.path, then it must be an absolute path
commonDir = os.path.abspath("./../common")
# Sys path is a list of system paths
sys.path.append(commonDir)

from CELMAPy.blobs import DriverBlobs

# Global data without encapsulation
xInd   = 26
yInd   = 16
zInd   = 0

pctPadding = 400
normed     = False
convertToPhysical = True
useSubProcess = False
plotAll = False

#{{{blobRadialFlux
def blobRadialFlux(*args):
    #{{{docstring
    """
    Plots the radialFlux.

    Parameters
    ----------
    See getBlobDriver for details.
    """
    #}}}

    dB = getBlobDriver(*args)

    dB.driverRadialFlux()
#}}}

#{{{blobWaitingTimePulsePlot
def blobWaitingTimePulsePlot(*args):
    #{{{docstring
    """
    Plots the waiting time and pulse width.

    Parameters
    ----------
    See getBlobDriver for details.
    """
    #}}}

    dB = getBlobDriver(*args)

    dB.driverWaitingTimePulse()
#}}}

#{{{blobTimeTracesPlot
def blobTimeTracesPlot(*args):
    #{{{docstring
    """
    Plots a timetrace of the blob and its individual bins.

    Parameters
    ----------
    See getBlobDriver for details.
    """
    #}}}

    dB = getBlobDriver(*args)

    dB.driverBlobTimeTraces()
#}}}

#{{{blob2DPlot
def blob2DPlot(*args, mode=None, fluct=None):
    #{{{docstring
    """
    Plots a timetrace of the blob and its individual bins.

    Parameters
    ----------
    *arg : positional arguments
        See getBlobDriver for details.
    mode : ["perp"|"par"|"pol"]
        The 2D mode to use.
    fluct : bool
        Whether or not to use fluctuation.
    """
    #}}}

    dB = getBlobDriver(*args)
    dB.setMode(mode)
    dB.setFluct(fluct)
    dB.driverPlot2DData()
#}}}

#{{{getBlobDriver
def getBlobDriver(dmp_folders    ,\
                  collectPaths   ,\
                  plotSuperKwargs,\
                  tSlice         ,\
                  condition      ,\
                  plotAll        ,\
                 ):
    #{{{docstring
    """
    Returns the blob driver.

    Parameters
    ----------
    dmp_folders : tuple
        Tuple of the dmp_folders.
    collectPaths : tuple
        Tuple of the paths to collect from.
    plotSuperKwargs : dict
        Keyword arguments for the plot super class.
    tSlice : slice
        How to slice the time.
    condition : float
        The condition in the conditional average will be set to
        flux.std()*condition
    plotAll : bool
           If True: The individual blobs will be plotted.

    Returns
    -------
    dB : DriverBlobs
        The blob driver.
    """
    #}}}

    slices = (xInd, yInd, zInd, tSlice)

    dB = DriverBlobs(\
                     # DriverBlobs
                     dmp_folders          ,\
                     slices               ,\
                     pctPadding           ,\
                     convertToPhysical    ,\
                     plotSuperKwargs      ,\
                     condition = condition,\
                     normed    = normed   ,\
                     plotAll   = plotAll  ,\
                     # DriverSuperClass
                     collectPaths  = collectPaths ,\
                     useSubProcess = useSubProcess,\
                    )

    return dB
#}}}

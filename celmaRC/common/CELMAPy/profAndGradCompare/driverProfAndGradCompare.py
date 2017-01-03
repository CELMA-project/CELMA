#!/usr/bin/env python

"""
Contains single driver and driver class for 1D fields
"""

from ..fields1D import CollectAndCalcFields1D
from ..collectAndCalcHelpers import (calcN, calcUIPar, calcUEPar,\
                                     polAvg, timeAvg)

#{{{driverProfAndGradCompare
def driverProfAndGradCompare(varName          ,\
                             collectPaths     ,\
                             steadyStatePath  ,\
                             convertToPhysical,\
                             xSlice           ,\
                             ySlice           ,\
                             zSlice           ,\
                             tSlice           ,\
                             plotSuperKwargs  ,\
                            ):
    #{{{docstring
    """
    Driver for plotting a single predefined fieldPlotType plot

    Parameters
    ----------
    varName : str
        Variable to collect.
    collectPaths : tuple
        Paths to collect from.
        The corresponind 't_array' of the paths must be in ascending order.
    steadyStatePath : str
        The steady state path.
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
    plotSuperKwargs : dict
        Keyword arguments for the plot super class
    """
    #}}}

    # Collect the steady state
    slices = (None, ySlice, zSlice, None)
    dict1D = collectWrapper((steadyStatePath,), varName,\
                            None, slices, convertToPhysical)

    # Extract the steady state variable at the last time
    steadyVar = dict1D[varName][-1,:]

    # Collect the variable
    # We set zSlice to 0 in order to avoid errors, note though that the
    # number is irrelevant as long as it is within the range
    slices = (None, ySlice, 0, None)
    dict1D = collectWrapper(collectPaths, varName,\
                            "polAndTimeAvg", slices,convertToPhysical)
    # Extract the average variable at the last time
    avgVar = dict1D[varName][-1,:]

    # Collect the fluctuations
    dict1D = collectWrapper(collectPaths, varName,\
                            "polAndTimeAvgFluct", slices, convertToPhysical)

    # Recast dict
    # Remove not needed
    dict1D.pop("time")
    dict1D.pop("thetaPos")
    # Find the standard deviation of the fluctuations
    varFluct = dict1D.pop(varName)
    stdDev   = np.sqrt(polAvg(timeAvg(varFluct**2.0)))

    dict1D["varName"]     = varName
    dict1D["steadyState"] = steadyVar
    dict1D["average"]     = avgVar
    dict1D["stdDev"]      = stdDev

# FIXME: You are here: Start plotting
    import pdb; pdb.set_trace()
    a =1
#    # Collect the fluctuating part
#    ccf1D = CollectAndCalcFields1D(collectPaths,\
#                                   mode = "radial",\
#                                   processing = None,\
#                                   convertToPhysical = convertToPhysical)
#
#    ccf1D.setSlice(xSlice, ySlice, zSlice, tSlice)
#
#    ccf1D.setVarName("lnN")
#    outDict = ccf1D.executeCollectAndCalc()
#    outDict.update({"n" : calcN(dict1D["lnN"],\
#                                not(ccf1D.convertToPhysical),\
#                                ccf1D.uc)})
#}}}

#{{{specialCollect
def specialCollect(varName, ccf1D):
    #{{{docstring
    """
    Calculate non-collectable variables

    Parameters
    ----------
    varName : str
        Name of the variable to calculate
    ccf1D : CollectAndCalcFields1D
        Object to perform helper collections

    Returns
    -------
    dict1D : dict
        Dictionary with the added variable.
        For details, see the documentation of the CollectAndCalcFields1D
        class.
    """
    #}}}

    ccf1D.setVarName("lnN")
    dict1D = ccf1D.executeCollectAndCalc()
    dict1D.update({"n" : calcN(dict1D["lnN"],\
                               not(ccf1D.convertToPhysical),\
                               ccf1D.uc)})
    if varName != "n":
        ccf1D.setVarName("momDensPar")
        dict1D.update(ccf1D.executeCollectAndCalc())
        dict1D.update({"uIPar":calcUIPar(dict1D["momDensPar"],dict1D["n"])})
    if varName == "uEPar":
        ccf1D.setVarName("jPar")
        dict1D.update(ccf1D.executeCollectAndCalc())
        dict1D.update({"uEPar": calcUEPar(dict1D["uIPar"],\
                                          dict1D["jPar"],\
                                          dict1D["n"],\
                                          not(ccf1D.convertToPhysical))})
    return dict1D
#}}}

#{{{collectWrapper
def collectWrapper(paths, varName, processing, slices, convertToPhysical):
    #{{{docstring
    """
    Wrapper function around the collect routine

    Parameters
    ----------
    paths : tuple
        Tuple of the collect paths
    varName : str
        Name of the variable to collect
    slices : tuple
        Tuple of the slices to use
    convertToPhysical : bool
        Whether or not to convert to physical units.

    Returns
    -------
    dict1D : dict
        Dictionary with the added variable.
        For details, see the documentation of the CollectAndCalcFields1D
        class.
    """
    #}}}
    ccf1D = CollectAndCalcFields1D(paths,\
                                   mode = "radial",\
                                   processing = processing,\
                                   convertToPhysical = convertToPhysical)

    ccf1D.setSlice(*slices)
    specialCollects = ("n", "uIPar", "uEPar")
    if not(varName in specialCollects):
        ccf1D.setVarName(varName)
        dict1D = ccf1D.executeCollectAndCalc()
    else:
        dict1D = specialCollect(varName, ccf1D)

    return dict1D
#}}}

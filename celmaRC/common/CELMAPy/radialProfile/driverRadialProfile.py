#!/usr/bin/env python

""" Contains single driver for radial profiles """

#{{{driverProfAndGradCompare
def driverProfAndGradCompare(varName          ,\
                             collectPaths     ,\
                             steadyStatePath  ,\
                             convertToPhysical,\
                             yInd             ,\
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
    yInd : int
        Fixed position in the parallel direction.
    plotSuperKwargs : dict
        Keyword arguments for the plot super class
    """
    #}}}

    ccRP = CollectAndCalcRadialProfile(yInd, convertToPhysical)

    # Collect the steady state variable
    dict1D = ccRP.collectWrapper((steadyStatePath,), varName)
    # Extract the steady state variable at the last time (but keep the
    # 4d)
    steadyVar = dict1D[varName][-2:-1,:,:,:]

    # Collect the variable
    dict1D = collectWrapper(collectPaths, varName)
    # Extract the variable
    var = dict1D[varName]

    varAvg, varFluct, varStd = ccPR.calcAvgFluctStd(var)

    # Calculate the derivatives
    dx = getUniformSpacing(steadyStatePath, "x")
    DDXSteadyVar = DDX(steadyVar, dx)
    DDXVar = DDX(var, dx)

    DDXVarAvg, DDXVarFluct, DDXVarStd = ccPR.calcAvgFluctStd(var)

    # Recast to dict
    # Remove not needed
    dict1D.pop("time")
    dict1D.pop("thetaPos")
    dict1D["varName"]   = varName
# FIXME: Check me
    import pdb; pdb.set_trace()
    dict1D["steadyVar"]    = steadyVar   [0,:,0,0]
    dict1D["DDXSteadyVar"] = DDXSteadyVar[0,:,0,0]
    dict1D["varAvg"]       = varAvg      [0,:,0,0]
    dict1D["DDXVarAvg"]    = DDXVarAvg   [0,:,0,0]
    dict1D["varAvgStd"]    = varAvgStd   [0,:,0,0]
    dict1D["DDXVarAvgStd"] = DDXVarAvgStd[0,:,0,0]

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

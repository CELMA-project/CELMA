#!/usr/bin/env python

""" Contains single driver for radial profiles """

from ..superClasses import DriverSuperClass
from ..collectAndCalcHelpers import DDX
from .collectAndCalcRadialProfile import CollectAndCalcRadialProfile
from .plotRadialProfile import PlotProfAndGradCompare

#{{{driverProfAndGradCompare
def driverProfAndGradCompare(varName          ,\
                             collectPaths     ,\
                             steadyStatePath  ,\
                             convertToPhysical,\
                             yInd             ,\
                             tSlice           ,\
                             plotSuperKwargs  ,\
                            ):
    #{{{docstring
    """
    Driver for plotting the comparison between profiles and gradients

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
    tSlice : slice
        How the data will be sliced in time
    plotSuperKwargs : dict
        Keyword arguments for the plot super class
    """
    #}}}

    ccrp = CollectAndCalcRadialProfile(yInd, tSlice, convertToPhysical)

    # Collect the steady state variable
    rp = ccrp.collectWrapper((steadyStatePath,), varName)
    # Extract the steady state variable at the last time (but keep the 4d)
    steadyVar = rp[varName][-2:-1,:,:,:]

    # Collect the variable
    rp = ccrp.collectWrapper(collectPaths, varName)
    # Extract the variable
    var = rp[varName]

    varAvg, varAvgFluct, varAvgStd = ccrp.calcAvgFluctStd(var)

    # Calculate the derivatives
    dx = ccrp.dh.dx
    DDXSteadyVar = DDX(steadyVar, dx)
    DDXVar = DDX(var, dx)

    DDXVarAvg, DDXVarAvgFluct, DDXVarAvgStd = ccrp.calcAvgFluctStd(DDXVar)

    # Recast to dict
    # Remove not needed
    rp.pop("time")
    rp.pop("thetaPos")
    rp["varName"]      = varName
    rp["steadyVar"]    = steadyVar   [0,:,0,0]
    rp["DDXSteadyVar"] = DDXSteadyVar[0,:,0,0]
    rp["varAvg"]       = varAvg      [0,:,0,0]
    rp["varAvgStd"]    = varAvgStd   [0,:,0,0]
    rp["DDXVarAvg"]    = DDXVarAvg   [0,:,0,0]
    rp["DDXVarAvgStd"] = DDXVarAvgStd[0,:,0,0]

    # Plot
    prp = PlotProfAndGradCompare(ccrp.uc, **plotSuperKwargs)
    prp.setData(rp)
    prp.plotSaveShowRadialProfiles()
#}}}

#{{{driverPosOfFluct
def driverPosOfFluct(var1Name         ,\
                     var2Name         ,\
                     collectPaths     ,\
                     steadyStatePath  ,\
                     convertToPhysical,\
                     yInd             ,\
                     tSlice           ,\
                     plotSuperKwargs  ,\
                    ):
    #{{{docstring
    """
    Driver for plotting the position of the amplitude of the
    fluctuations.

    Parameters
    ----------
    var1Name : str
        Variable to plot the profile, derivative and position of
        fluctuation of.
    var2Name : str
        Variable to plot position of fluctuation of.
    collectPaths : tuple
        Paths to collect from.
        The corresponind 't_array' of the paths must be in ascending order.
    steadyStatePath : str
        The steady state path.
    convertToPhysical : bool
        Whether or not to convert to physical units.
    yInd : int
        Fixed position in the parallel direction.
    tSlice : slice
        How the data will be sliced in time
    plotSuperKwargs : dict
        Keyword arguments for the plot super class
    """
    #}}}

    ccrp = CollectAndCalcRadialProfile(yInd, tSlice, convertToPhysical)
    # Collect the steady state variable
    rp = ccrp.collectWrapper((steadyStatePath,), var1Name)
    # Extract the steady state variable at the last time (but keep the 4d)
    steadyVar = rp[var1Name][-2:-1,:,:,:]

    # Collect first variable
    rp = ccrp.collectWrapper(collectPaths, var1Name)
    # Extract the variable
    var = rp[var1Name]
    varAvg, varAvgFluct, varAvgStd = ccrp.calcAvgFluctStd(var)

    # Calculate the derivatives of first variable
    dx = getUniformSpacing(steadyStatePath, "x")[0,0]
    DDXSteadyVar = DDX(steadyVar, dx)
    DDXVar = DDX(var, dx)

    DDXVarAvg, DDXVarAvgFluct, DDXVarAvgStd = ccrp.calcAvgFluctStd(DDXVar)

    # Collect phi
    rp = ccrp.collectWrapper(collectPaths, var2Name)
    # Extract the variable
    var = rp[var2Name]
    _, _, var2AvgStd = ccrp.calcAvgFluctStd(var)

    # Recast to dict
    # Remove not needed
    rp.pop("time")
    rp.pop("thetaPos")
    rp["varName"]      = var1Name
    rp["var2Name"]     = var2Name
    rp["steadyVar"]    = steadyVar   [0,:,0,0]
    rp["DDXSteadyVar"] = DDXSteadyVar[0,:,0,0]
    rp["varAvg"]       = varAvg      [0,:,0,0]
    rp["varAvgStd"]    = varAvgStd   [0,:,0,0]
    rp["var2AvgStd"]   = var2AvgStd  [0,:,0,0]
    rp["DDXVarAvg"]    = DDXVarAvg   [0,:,0,0]
    rp["DDXVarAvgStd"] = DDXVarAvgStd[0,:,0,0]

    # Plot
    prp = PlotProfAndGradCompare(ccrp.uc, **plotSuperKwargs)
    prp.setData(rp)
    prp.setDataForSecondVar(rp)
    prp.plotSaveShowPosOfFluct()
#}}}

#{{{DriverRadialProfile
class DriverRadialProfile(DriverSuperClass):
    """
    Class for driving of the plotting of the radial profiles.
    """

    #{{{Constructor
    def __init__(self                     ,\
                 dmp_folders              ,\
                 steadyStatePath          ,\
                 yInd                     ,\
                 tSlice                   ,\
                 plotSuperKwargs          ,\
                 varName           = "n"  ,\
                 var2Name          = "phi",\
                 convertToPhysical = True ,\
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
        steadyStatePath : str
            The steady state path.
        yInd : int
            Fixed position in the parallel direction.
        tSlice : slice
            How the data will be sliced in time.
        plotSuperKwargs : dict
            Keyword arguments for the plot super class.
        varName : str
            Variable to plot the profile, derivative and position of
            fluctuation of.
        var2Name : str
            Variable to plot position of fluctuation of.
        convertToPhysical : bool
            Whether or not to convert to physical units.
        **kwargs : keyword arguments
            See parent class for details.
        """
        #}}}

        # Call the constructor of the parent class
        super().__init__(dmp_folders, **kwargs)

        # Set the member data
        self._varName          = varName
        self._var2Name         = var2Name
        self._steadyStatePath  = steadyStatePath
        self._yInd             = yInd
        self._tSlice           = tSlice
        self.convertToPhysical = convertToPhysical

        # Update the plotSuperKwargs dict
        plotSuperKwargs.update({"dmp_folders":dmp_folders})
        plotSuperKwargs.update({"plotType"   :"radialProfiles"})
        self._plotSuperKwargs = plotSuperKwargs
    #}}}

    #{{{driverProfAndGradCompare
    def driverProfAndGradCompare(self):
        #{{{docstring
        """
        Wrapper to driverProfAndGradCompare
        """
        #}}}
        args =  (\
                 self._varName         ,\
                 self._collectPaths    ,\
                 self._steadyStatePath ,\
                 self.convertToPhysical,\
                 self._yInd            ,\
                 self._tSlice          ,\
                 self._plotSuperKwargs ,\
                )
        if self._useSubProcess:
            processes = Process(target = driverProfAndGradCompare, args = args)
            processes.start()
        else:
            driverProfAndGradCompare(*args)
    #}}}

    #{{{driverPosOfFluct
    def driverPosOfFluct(self):
        #{{{docstring
        """
        Wrapper to driverPosOfFluct
        """
        #}}}
        args =  (\
                 self._varName         ,\
                 self._var2Name        ,\
                 self._collectPaths    ,\
                 self._steadyStatePath ,\
                 self.convertToPhysical,\
                 self._yInd            ,\
                 self._tSlice          ,\
                 self._plotSuperKwargs ,\
                )
        if self._useSubProcess:
            processes = Process(target = driverPosOfFluct, args = args)
            processes.start()
        else:
            driverPosOfFluct(*args)
    #}}}
#}}}

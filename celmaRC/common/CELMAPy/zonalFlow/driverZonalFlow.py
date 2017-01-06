#!/usr/bin/env python

""" Contains single driver for zonal flows """

from ..superClasses import DriverSuperClass
from ..collectAndCalcHelpers import DDX
from ..radialProfile import CollectAndCalcRadialProfile
from .collectAndCalcZonalFlow import CollectAndCalcZonalFlow
from .plotZonalFlow import PlotZonalFlow
from multiprocessing import Process
import numpy as np

#{{{driverZonalFlow
def driverZonalFlow(\
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

    cczf = CollectAndCalcZonalFlow(yInd, tSlice, convertToPhysical)

    # Calculate the steady state variable
    polExBSS = cczf.calcPoloidalExB((steadyStatePath,))

    # Set rho and dx (obtained after call to calcPoloidalExB)
    rho = cczf.dh.rho
    dx = cczf.dh.dx

    # Extract the steady state variable at the last time (but keep the 4d)
    sSVar = polExBSS["uExBPoloidal"][-2:-1,:,:,:]
    # Expand rho in order to have a defined arithmetic operation with sSVar
    tmp = np.expand_dims(\
            np.expand_dims(\
                np.expand_dims(rho, axis=0),\
            axis=2),\
          axis=3)
    rho    = np.empty(sSVar.shape)
    rho[:] = tmp
    # Calculate the angular velocity
    angVelSSVar =sSVar/rho
    # Calculate the shear
    angVelSSVarShear = DDX(angVelSSVar, dx)

    # Calculate the poloidal ExB in the saturarted turbulence state
    polExB = cczf.calcPoloidalExB(collectPaths)
    var = polExB["uExBPoloidal"]
    # Calculate average, fluct and std
    varAvg, _, varStd = CollectAndCalcRadialProfile.calcAvgFluctStd(var)
    # Expand rho in order to have a defined arithmetic operation with sSVar
    rho    = np.empty(varAvg.shape)
    rho[:] = tmp
    # Calculate the angular velocity (the error propagation is linear)
    angVelVar    = var   /rho
    angVelVarAvg = varAvg/rho
    angVelVarStd = varStd/rho

    # Calculate the shear
    angVelShearVar = DDX(angVelVar, dx)
    # Calculate average, fluct and std of the shear
    angVelShearVarAvg, _, angVelShearVarStd =\
            CollectAndCalcRadialProfile.calcAvgFluctStd(angVelShearVar)

    # Recast to dict
    # Remove not needed
    polExB.pop("time")
    polExB["polExBSS"]          = sSVar            [0,:,0,0]
    polExB["angPolExBSS"]       = angVelSSVar      [0,:,0,0]
    polExB["angPolExBShearSS"]  = angVelSSVarShear [0,:,0,0]
    polExB["polExBAvg"]         = varAvg           [0,:,0,0]
    polExB["polExBStd"]         = varStd           [0,:,0,0]
    polExB["angPolExBAvg"]      = angVelVarAvg     [0,:,0,0]
    polExB["angPolExBStd"]      = angVelVarStd     [0,:,0,0]
    polExB["angPolExBShearAvg"] = angVelShearVarAvg[0,:,0,0]
    polExB["angPolExBShearStd"] = angVelShearVarStd[0,:,0,0]
    polExB["rho"]               = rho              [0,:,0,0]

    # Plot
    prp = PlotZonalFlow(cczf.uc, **plotSuperKwargs)
    prp.setData(polExB)
    prp.plotSaveShowZonalFlow()
#}}}

#{{{DriverZonalFlow
class DriverZonalFlow(DriverSuperClass):
    """
    Class for driving of the plotting of the zonal flows.
    """

    #{{{Constructor
    def __init__(self                     ,\
                 dmp_folders              ,\
                 steadyStatePath          ,\
                 yInd                     ,\
                 tSlice                   ,\
                 plotSuperKwargs          ,\
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
        convertToPhysical : bool
            Whether or not to convert to physical units.
        **kwargs : keyword arguments
            See parent class for details.
        """
        #}}}

        # Call the constructor of the parent class
        super().__init__(dmp_folders, **kwargs)

        # Set the member data
        self._steadyStatePath  = steadyStatePath
        self._yInd             = yInd
        self._tSlice           = tSlice
        self.convertToPhysical = convertToPhysical

        # Update the plotSuperKwargs dict
        plotSuperKwargs.update({"dmp_folders":dmp_folders})
        plotSuperKwargs.update({"plotType"   :"zonalFlows"})
        self._plotSuperKwargs = plotSuperKwargs
    #}}}

    #{{{driverZonalFlow
    def driverZonalFlow(self):
        #{{{docstring
        """
        Wrapper to driverPosOfFluct
        """
        #}}}
        args =  (\
                 self._collectPaths    ,\
                 self._steadyStatePath ,\
                 self.convertToPhysical,\
                 self._yInd            ,\
                 self._tSlice          ,\
                 self._plotSuperKwargs ,\
                )
        if self._useSubProcess:
            processes = Process(target = driverZonalFlow, args = args)
            processes.start()
        else:
            driverZonalFlow(*args)
    #}}}
#}}}

#!/usr/bin/env python

""" Contains single driver for poloidal flows """

from ..superClasses import DriverSuperClass
from ..collectAndCalcHelpers import DDX
from ..radialProfile import CollectAndCalcRadialProfile
from .collectAndCalcPoloidalFlow import CollectAndCalcPoloidalFlow
from .plotPoloidalFlow import PlotPoloidalFlow
from multiprocessing import Process
import numpy as np

#{{{driverPoloidalFlow
def driverPoloidalFlow(\
                    collectPaths     ,\
                    steadyStatePath  ,\
                    convertToPhysical,\
                    yInd             ,\
                    tSlice           ,\
                    mode             ,\
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
    mode : ["plain"|"wShear"]
        If "plain", only the poloidal ExB should will be plotted.
        If "wShear", only the angular frequency is plotted together with
        its shear.
    plotSuperKwargs : dict
        Keyword arguments for the plot super class
    """
    #}}}

    # Initialize the collector
    cczf = CollectAndCalcPoloidalFlow(yInd, tSlice, convertToPhysical)

    # Calculate the steady state variable
    polExBSS = cczf.calcPoloidalExB((steadyStatePath,))

    # Obtain rho (must be after calcPoloidalExB)
    rho = cczf.dh.rho

    # Extract the steady state variable at the last time (but keep the 4d)
    sSVar = polExBSS["uExBPoloidal"][-1:,:,:,:]

    # Calculate the poloidal ExB in the saturarted turbulence state
    polExB = cczf.calcPoloidalExB(collectPaths)
    var = polExB["uExBPoloidal"]
    # Calculate average, fluct and std
    varAvg, _, varStd = CollectAndCalcRadialProfile.calcAvgFluctStd(var)

    # Remove not needed
    polExB.pop("time")
    # Recast to dict
    polExB["polExBSS"]  = sSVar [0,:,0,0]
    polExB["polExBAvg"] = varAvg[0,:,0,0]
    polExB["polExBStd"] = varStd[0,:,0,0]
    polExB["rho"]       = rho

    if mode == "plain":
        # Plot
        prp = PlotPoloidalFlow(cczf.uc, **plotSuperKwargs)
        prp.setData(polExB)
        prp.plotSaveShowPoloidalFlow()

    elif mode == "wShear":
        # Obtain dx
        dx = cczf.dh.dx

        # Expand rho in order to have a defined arithmetic operation with sSVar
        tmp = np.expand_dims(\
                np.expand_dims(\
                    np.expand_dims(rho, axis=0),\
                axis=2),\
              axis=3)
        rho    = np.empty(sSVar.shape)
        rho[:] = tmp
        # Calculate the angular frequency
        angFreqSSVar =sSVar/rho
        # Calculate the shear
        angFreqSSVarShear = DDX(angFreqSSVar, dx)

        # Expand rho in order to have a defined arithmetic operation with sSVar
        rho    = np.empty(varAvg.shape)
        rho[:] = tmp
        # Calculate the angular frequency (the error propagation is linear)
        angFreqVar    = var   /rho
        angFreqVarAvg = varAvg/rho
        angFreqVarStd = varStd/rho

        # Calculate the shear
        angFreqShearVar = DDX(angFreqVar, dx)
        # Calculate average, fluct and std of the shear
        angFreqShearVarAvg, _, angFreqShearVarStd =\
                CollectAndCalcRadialProfile.calcAvgFluctStd(angFreqShearVar)

        # Recast to dict
        polExB["angPolExBSS"]       = angFreqSSVar      [0,:,0,0]
        polExB["angPolExBShearSS"]  = angFreqSSVarShear [0,:,0,0]
        polExB["angPolExBAvg"]      = angFreqVarAvg     [0,:,0,0]
        polExB["angPolExBStd"]      = angFreqVarStd     [0,:,0,0]
        polExB["angPolExBShearAvg"] = angFreqShearVarAvg[0,:,0,0]
        polExB["angPolExBShearStd"] = angFreqShearVarStd[0,:,0,0]

        # Plot
        prp = PlotPoloidalFlow(cczf.uc, **plotSuperKwargs)
        prp.setDataWShear(polExB)
        prp.plotSaveShowPoloidalFlowWShear()
#}}}

#{{{DriverPoloidalFlow
class DriverPoloidalFlow(DriverSuperClass):
    """
    Class for driving of the plotting of the poloidal flows.
    """

    #{{{Constructor
    def __init__(self                       ,\
                 dmp_folders                ,\
                 steadyStatePath            ,\
                 yInd                       ,\
                 tSlice                     ,\
                 plotSuperKwargs            ,\
                 mode              = "plain",\
                 convertToPhysical = True   ,\
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
        mode : ["plain"|"wShear"]
            If "plain", only the poloidal ExB should will be plotted.
            If "wShear", only the angular frequency is plotted together with
            its shear.
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
        self._mode             = mode
        self.convertToPhysical = convertToPhysical

        # Update the plotSuperKwargs dict
        plotSuperKwargs.update({"dmp_folders":dmp_folders})
        plotSuperKwargs.update({"plotType"   :"poloidalFlows"})
        self._plotSuperKwargs = plotSuperKwargs
    #}}}

    #{{{driverPoloidalFlow
    def driverPoloidalFlow(self):
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
                 self._mode            ,\
                 self._plotSuperKwargs ,\
                )
        if self._useMultiProcess:
            processes = Process(target = driverPoloidalFlow, args = args)
            processes.start()
        else:
            driverPoloidalFlow(*args)
    #}}}
#}}}

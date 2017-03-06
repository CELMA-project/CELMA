#!/usr/bin/env python

"""
Contains class to calculate the phase shift between n and phi.
"""

from ..collectAndCalcHelpers import getScanValue
from ..fourierModes import CollectAndCalcFourierModes
from ..timeTrace import getTimeTrace
from .collectAndCalcAnalyticGrowthRates import\
     CollectAndCalcAnalyticGrowthRates
from scipy import signal
import pandas as pd
import numpy as np

#{{{CollectAndCalcPhaseShift
class CollectAndCalcPhaseShift(object):
    """
    Class for collecting and calculating the phase shifts
    """

    #{{{constructor
    def __init__(self):
        """
        This constructor will set notSet parameters
        """

        self._notCalled = ["setCollectArgs", "setGetDataArgs"]
    #}}}

    #{{{setCollectArgs
    def setCollectArgs(self,\
                       scanCollectPaths,\
                       steadyStatePaths,\
                       tSlices,\
                       scanParameter):
        #{{{docstring
        """
        Sets the arguments used to collect.

        NOTE: These arguments can be made with
               DriverGrowthRates.makeCollectArgs

        Parameters
        ----------
        scanCollectPaths : tuple of tuple of strings
            One tuple of strings for each value of the scan which will
            be used in collective collect.
        steadyStatePaths : tuple
            Path to the steady state simulation.
            The tuple must be ordered according to the scan order in
            scanCollectPaths.
        tSlices : tuple of slices
            The time slices to use for each folder in the scan order.
            This must manually be selected to the linear phase found in
            from the fourier moedes.
            The tuple must be ordered according to the scan order in
            scanCollectPaths.
        scanParameter : str
            String segment representing the scan
        """
        #}}}

        if "setCollectArgs":
            self._notCalled.remove("setCollectArgs")

        # Set the member data
        self._scanCollectPaths = scanCollectPaths
        self._steadyStatePaths = steadyStatePaths
        self._tSlices          = tSlices
        self._scanParameter    = scanParameter
    #}}}

    #{{{setGetDataArgs
    def setGetDataArgs(self,\
                       convertToPhysical,\
                       indicesArgs      ,\
                       indicesKwargs    ,\
                       ):
        #{{{docstring
        """
        Sets additional parameters used in self.getData.

        NOTE: These arguments can be made with
               DriverGrowthRates.makeGetDataArgs[1:-1]

        Parameters
        ----------
        convertToPhysical : bool
            Whether or not to convert to physical units.
        indicesArgs : tuple
            Contains (xInd, yInd, zInd)
            NOTE: Only one spatial point should be used.
        indicesKwargs : dict
            Keyword arguments to use when setting the indices for
            collecting.
            Contains the keys
                * "tSlice"          - Will be updated
                * "nPoints"         - Should only be 1
                * "equallySpace"    - Should be "x"
                * "steadyStatePath" - The steady state path
        """
        #}}}

        if "setGetDataArgs":
            self._notCalled.remove("setGetDataArgs")

        # Set the data
        self._convertToPhysical = convertToPhysical
        self._indicesArgs       = indicesArgs
        self._yInd              = indicesArgs[1]
        self._indicesKwargs     = indicesKwargs
    #}}}

    #{{{getData
    def getData(self):
        #{{{docstring
        """
        Makes a DataFrame of the phase shifts obtained from the simulations.

        Returns
        -------
        phaseShiftDataFrame : DataFrame
            DataFrame consisting of the variables (measured properties):
                * "phaseShift"
            over the observation "Scan" over the observation "modeNr"
        phaseShiftDataFrame : DataFrame
            DataFrame consisting of the variables (measured properties):
                * "phaseShift"
            over the observation "Scan"
        positionTuple : tuple
            The tuple containing (rho, theta, z).
            Needed in the plotting routine.
        uc : Units Converter
            The units converter used when obtaining the fourier modes.
            Needed in the plotting routine.
        """
        #}}}

        # Guard
        if len(self._notCalled) > 0:
            message = "The following functions were not called:\n{}".\
                        format("\n".join(self._notCalled))
            raise RuntimeError(message)

        # Collect the analytic phase shift
        # Create collect object
        ccagr = CollectAndCalcAnalyticGrowthRates(self._steadyStatePaths,\
                                                  self._scanParameter,\
                                                  self._yInd)
        # Obtain the data
        analyticalGRDataFrame, _, _, uc =\
            ccagr.getData()

        # Recast the data frame
        analyticalGRDataFrame.drop("growthRate", axis=1, inplace=True)
        analyticalGRDataFrame.drop("angularFrequency", axis=1, inplace=True)
        analyticalPhaseShiftDataFrame =\
            analyticalGRDataFrame.\
                rename(columns={"phaseShiftNPhi":"phaseShift"})

        # Get the levels correct
        analyticalPhaseShiftDataFrame = analyticalPhaseShiftDataFrame.swaplevel()

        # Recast the data frame
        phaseShiftDataFrame = 0

        dataFrameDict = {"phaseShift":[]}
        scanValues    = []

        # Collect the phase shift from the simulations
        loopOver = zip(self._scanCollectPaths,\
                       self._steadyStatePaths,\
                       self._tSlices         ,\
                       )

        # Loop over the folders
        for scanPaths, steadyStatePath, tSlice in loopOver:

            # Obtain the scan value
            scanValue = getScanValue(scanPaths, self._scanParameter)

            # Update with the correct tSlice
            self._indicesKwargs.update({"tSlice" : tSlice})

            # Obtain teh time traces
            n, phi, positionTuple = self._getTimeTraces(scanPaths)

            # Obtain the cross spectral density
            # NOTE: If this is below the number of samples, a smoothing
            #       will occur
            nperseg = len(n)
            # NOTE: The triangular window corresponds to the periodogram
            #       estimate of the spectral density
            # NOTE: The first output (frequency) is not used
            _, csd = signal.csd(n, phi, window="triang", nperseg=nperseg)

            maxInd = self._getMaxIndOfMagnitude(csd)

            # FIXME: Not sure why, but there seem to be a sign error
            avgPhaseShiftNPhi = -np.angle(csd[maxInd])

            scanValues.append(scanValue)
            dataFrameDict["phaseShift"].append(avgPhaseShiftNPhi)

        # Make the data frame
        phaseShiftDataFrame = pd.DataFrame(dataFrameDict, index=scanValues)
        phaseShiftDataFrame.index.name = self._scanParameter

        return analyticalPhaseShiftDataFrame,\
               phaseShiftDataFrame,\
               positionTuple,\
               uc
    #}}}

    #{{{_getTimeTraces
    def _getTimeTraces(self, collectPaths):
        #{{{docstring
        """
        Returns the n and phi time traces

        Parameters
        ----------
        collectPaths : tuple
            Tuple of strings containing the paths to collect from.

        Returns
        -------
        n : array 1-d
            The time trace of n.
        phi : array 1-d
            The time trace of phi.
        positionTuple : tuple
            Tuple with the rho, theta and z value.
        """
        #}}}
        mode = "fluct"

        tt = tuple(\
                   getTimeTrace(collectPaths           ,\
                                varName                ,\
                                self._convertToPhysical,\
                                mode                   ,\
                                self._indicesArgs      ,\
                                self._indicesKwargs    ,\
                               )\
                   for varName in ("n", "phi")\
                   )

        # Extract n, phi and the position tuple
        nDict         = tt[0][0]
        positionTuple = list(nDict.keys())[0]
        n             = nDict[positionTuple]["n"]
        phi           = tt[1][0][positionTuple]["phi"]
        # Cast the string to a tuple
        positionTuple = tuple(float(val) for val in positionTuple.split(","))

        return n, phi, positionTuple
    #}}}

    #{{{_getMaxIndOfMagnitude
    def _getMaxIndOfMagnitude(self, csd):
        #{{{
        """
        Returns the max index of the magnitude of the cross spectral density.

        Parameters
        ----------
        csd : array 1-d
            The cross spectral density.

        Returns
        -------
        maxInd : int
            The index at the maximum amplitude
        """
        #}}}

        # Find the magnitude by first casting csd into a format which
        # CollectAndCalcGrowthRates.calcMagnitude(csdDict) accepts

        # Expand the time dimension (only one point)
        csd = np.expand_dims(csd,axis=0)
        # Put into dict so in can be used in calcMagnitude
        # NOTE: The position is not required here, so we use "0" as a
        #       place holder
        csdDict = {"0":{"csd":csd}}

        # Obtain the magnitudes
        csdDictWMagnitues = CollectAndCalcFourierModes.calcMagnitude(csdDict)

        # Extract Magnitudes
        magnitudes = csdDictWMagnitues["0"]["csdMagnitude"]

        # NOTE: Only first occurence is returned from numpy argmax
        maxInd = np.argmax(magnitudes)

        return maxInd
    #}}}
#}}}

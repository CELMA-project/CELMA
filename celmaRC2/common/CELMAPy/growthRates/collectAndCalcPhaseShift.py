#!/usr/bin/env python

"""
Contains class to calculate the phase shift between n and phi.
"""

from ..fourierModes import CollectAndCalcFourierModes
from ..timeTrace import getTimeTrace
from scipy import signal
import pandas as pd

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
    def setCollectArgs(scanCollectPaths,\
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
    def setGetDataArgs(convertToPhysical,\
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
            Contains (xInd, yInd)
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
            over the observation "modeNr" over the observation "Scan"
        phaseShiftDataFrame : DataFrame
            DataFrame consisting of the variables (measured properties):
                * "phaseShift"
            over the observation "Scan"
        positionTuple : tuple
            The tuple containing (rho, z).
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
        analyticalGRDataFrame, _, positionTuple, uc =\
            ccagr.getData()

        # Recast the data frame
        import pdb; pdb.set_trace()
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

            import pdb; pdb.set_trace()
            # FIXME: What are these scanPaths???
            n, phi = self._getTimeTraces(scanPaths)

            # Obtain the cross spectral density
            # NOTE: The triangular window corresponds to the periodogram
            #       estimate of the spectral density
            csd = signal.csd(n, phi, window="triang")

            maxInd = self._getMaxIndOfMagnitude(csd)

            avgPhaseShiftNPhi = np.angle(csd[maxInd])

            scanValues.append(scanValue)
            dataFrameDict["phaseShift"].append(avgPhaseShiftNPhi)

        # Make the data frame
        phaseShiftDataFrame = pd.DataFrame(dataFrameDict, index=scanValues)
        phaseShiftDataFrame.index.name = self._scanParameter

        return analyticalPhaseShiftDataFrame, phaseShiftDataFrame, positionTuple, uc
    #}}}

    #{{{_getTimeTraces
    def _getTimeTraces(self):
        #{{{
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

        # Extract n and phi
        n   = tt[0][0]
        phi = tt[1][0]

        return n, phi
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
        magnitudes = csdDict["0"]["csdMagnitude"]

        # NOTE: Only first occurence is returned from numpy argmax
        maxInd = np.argmax(magnitudes)

        return maxInd
    #}}}
#}}}

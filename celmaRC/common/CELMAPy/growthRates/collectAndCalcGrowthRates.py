#!/usr/bin/env python

"""
Contains class to calculate the and the growth rate and angular velocity
of a time trace of a spatial FFT.
"""

from ..collectAndCalcHelpers import linRegOfExp, getScanValue
from ..fourierModes import CollectAndCalcFourierModes
import pandas as pd

#{{{CollectAndCalcGrowthRates
class CollectAndCalcGrowthRates(object):
    """
    Class for collecting and calculating the growth rates
    """

    #{{{constructor
    def __init__(self            ,\
                 scanCollectPaths,\
                 steadyStatePaths,\
                 startInds       ,\
                 endInds         ,\
                 scanParameter):
        #{{{docstring
        """
        This constructor will:
            * Set the member data

        Parameters
        ----------
        scanCollectPaths : tuple of tuple of strings
            One tuple of strings for each value of the scan which will
            be used in collective collect.
        steadyStatePaths : tuple
            Path to the steady state simulation.
            The tuple must be ordered according to the scan order in
            scanCollectPaths.
        startInds : tuple
            Index indicating what time index to start the calculations
            from.
            The zeroth index is the first index of the first path being
            used in the tuple which will be used in the collective
            collect.
            The tuple must be ordered according to the scan order in
            scanCollectPaths.
        endInds : tuple
            Index indicating what time index to end the calculations
            at.
            The zeroth index is the first index of the first path being
            used in the tuple which will be used in the collective
            collect.
            The tuple must be ordered according to the scan order in
            scanCollectPaths.
        scanParameter : str
            String segment representing the scan
        """
        #}}}

        # Set the member data
        self._scanCollectPaths = scanCollectPaths
        self._steadyStatePaths = steadyStatePaths
        self._startInds        = startInds
        self._endInds          = endInds
        self._scanParameter    = scanParameter
    #}}}

    #{{{calcSlopeAndSpread
    @staticmethod
    def calcSlopeAndSpread(magnitudes, time, startInd=None, endInd=None):
        #{{{docstring
        """
        Calculates the slope and the spread for the given input.

        Parameters
        ----------
        magnitudes : array 2d
            Array of the time trace of the magnitudes of a fourier
            transformed poloidal profile on the form (t, mode). Where
            the first mode in the array is the 1st mode (i.e. the 0th
            mode has been excluded)
        time : array
            The time corresponding to the timeTrace of the magnitudes
        startInd : [None|int]
            Start index to calcuate from.
        endInd : [None|int]
            End index to calcuate to.

        Returns
        -------
        slopes : tuple
            Tuple of all the calculated slopes.
            The tuple is ordered after ascending modenumber (starting at 1).
            The slope is calculated using linear regression of an exponential.
        spreads : tuple
            Tuple of all the calculated spread.
            The tuple is ordered after ascending modenumber (starting at 1).
            The spread is calculated using linear regression of an exponential.
        """
        #}}}

        # Slice the magnitudes and the time
        magnitudes  = magnitudes[startInd:endInd, :]
        time        = time      [startInd:endInd]
        modeIndices = range(magnitudes.shape[1])

        slopes  = []
        spreads = []

        for modeNr in modeIndices:
            # Calculate the slope and the spread
            curSlope, curSpread = linRegOfExp(time, magnitudes[:,modeNr])
            slopes .append(curSlope)
            spreads.append(curSlope)

        return tuple(slopes), tuple(spreads)
    #}}}

    #{{{calcAvgAngularVelocityAndSpread
    @staticmethod
    def calcAvgAngularVelocityAndSpread(angularVelocity,\
                                        startInd=None,\
                                        endInd=None):
        #{{{docstring
        """
        Calculates the average angular velocity and the spread for the
        given input.

        Parameters
        ----------
        angularVelocity : array 2d
            Array of the time trace of the angular velocity of a fourier
            transformed poloidal profile on the form (t-1, mode). Where
            the first mode in the array is the 1st mode (i.e. the 0th
            mode has been excluded)
        startInd : [None|int]
            Start index to calcuate from.
        endInd : [None|int]
            End index to calcuate to.

        Returns
        -------
        avgAngVel : tuple
            Tuple of all the calculated averaged angular velocities.
            The tuple is ordered after ascending modenumber (starting at 1).
        spread : tuple
            Tuple of all the calculated spread.
            The tuple is ordered after ascending modenumber (starting at 1).
        """
        #}}}

        # Slice the angularVelocity
        angularVelocity  = angularVelocity[startInd:endInd, :]
        modeIndices      = range(angularVelocity.shape[1])

        avgAngVels = []
        spreads    = []

        for modeNr in modeIndices:
            # Calculate the average and the spread
            avgAngVels.append(angularVelocity.mean())
            spreads   .append(angularVelocity.std())

        return tuple(avgAngVels), tuple(spreads)
    #}}}

    #{{{makeDataFrame
    def makeDataFrame(self             ,\
                      varName          ,\
                      convertToPhysical,\
                      indicesArgs      ,\
                      indicesKwargs    ,\
                      nModes = 7       ,\
                      ):
        #{{{docstring
        """
        Makes a DataFrame of the growth rates and angular velocities.

        Parameters
        ----------
        varName : str
            Name of variable to find the growth rates of.
        convertToPhysical : bool
            Whether or not to convert to physical units.
        indicesArgs : tuple
            Tuple of indices to use when collecting.
            NOTE: Only one spatial point should be used.
        indicesKwargs : dict
            Keyword arguments to use when setting the indices for
            collecting.
            NOTE: Only one spatial point should be used.
        nModes : int
            Number of modes.

        Returns
        -------
        growthRateDataFrame : DataFrame
            DataFrame consisting of the variables (measured properties):
                * "growthRate"
                * "growthRateStd"
                * "averageAngularVelocity"
                * "averageAngularVelocityStd"
            over the observation "Mode number" over the observation
            "Scan"
        """
        #}}}

        # SEE:
        # https://github.com/pandas-dev/pandas/blob/master/doc/cheatsheet/Pandas_Cheat_Sheet.pdf
        # http://stackoverflow.com/questions/10715965/add-one-row-in-a-pandas-dataframe

        # Something like this
        # How to know number of modes apriori?

        multiTuples = []
        fullDict = {"growthRate":[], "growthRateStd":[],\
                    "averageAngularVelocity":[], "averageAngularVelocityStd":[]}

        loopOver = zip(self._scanCollectPaths,\
                       self._steadyStatePaths,\
                       self._startInds       ,\
                       self._endInds         ,\
                       )

        # Loop over the folders
        for scanPaths, steadyStatePath, startInd, endInd in loopOver:
            # Obtain the scan value
            scanValue = getScanValue(scanPaths, self._scanParameter)

            # Obtain the magnitudes and angular velocities
            ccfm = CollectAndCalcFourierModes(\
                                scanPaths                         ,\
                                convertToPhysical = convertToPhysical,\
                                             )
            # Set the slice
            indicesKwargs.update({"steadyStatePath" : steadyStatePath})
            ccfm.setIndices(*indicesArgs, **indicesKwargs)

            # Set name
            ccfm.setVarName(varName)

            # Execute the collection
            fm = ccfm.executeCollectAndCalc()
            fm = ccfm.convertTo2D(fm)
            fm = ccfm.calcMagnitude(fm)
            fm = ccfm.calcAngularVelocity(fm)

            # NOTE: We skip the offset mode.
            #       Thus, we add 1 in the range in order to look at
            #       nModes modes
            modeStart = 1
            modeEnd   = nModes+1

            # Get the keys
            firstKey = tuple(fm.keys())[0]

            # Obtain the time, magitude and the angular velocity
            time            = fm[firstKey]["time"]
            magnitudes      = fm[firstKey][varName+"Magnitude"]
            angularVelocity = fm[firstKey][varName+"AngularVelocity"]

            slope, slopeStd =\
                self.calcSlopeAndSpread(magnitudes[modeStart:modeEnd],\
                                        time                         ,\
                                        startInd = startInd          ,\
                                        endInd   = endInd            ,\
                                       )

            avgAngVel, avgAngVelStd =\
                self.calcAvgAngularVelocityAndSpread(\
                                        angularVelocity[modeStart:modeEnd],\
                                        startInd=startInd                 ,\
                                        endInd=endInd                     ,\
                                        )

            for modeInd in range(len(slope)):
                # Fill the multiIndexTuple and the dict
                multiTuples.append((scanValue, modeInd + 1))

                fullDict["growthRate"               ] = slope       [modeInd]
                fullDict["growthRateStd"            ] = slopeStd    [modeInd]
                fullDict["averageAngularVelocity"   ] = avgAngVel   [modeInd]
                fullDict["averageAngularVelocityStd"] = avgAngVelStd[modeInd]


        # Make the data frame
        growthRateDataFrame =\
            pd.DataFrame(fullDict,\
                         index=pd.MultiIndex.from_tuples(\
                            multiTuples, names=(self._scanParameter,"modes")))

        return growthRateDataFrame
    #}}}
#}}}

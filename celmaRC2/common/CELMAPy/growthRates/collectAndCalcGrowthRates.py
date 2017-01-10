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
                 tSlices         ,\
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

        # Set the member data
        self._scanCollectPaths = scanCollectPaths
        self._steadyStatePaths = steadyStatePaths
        self._tSlices          = tSlices
        self._scanParameter    = scanParameter
    #}}}

    @staticmethod
    #{{{calcSlopeAndSpread
    def calcSlopeAndSpread(magnitudes, time):
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
        modeIndices = range(magnitudes.shape[1])

        slopes  = []
        spreads = []

        for modeNr in modeIndices:
            # Calculate the slope and the spread
            curSlope, curSpread = linRegOfExp(time, magnitudes[:,modeNr])
            slopes .append(curSlope)
            spreads.append(curSpread)

        return tuple(slopes), tuple(spreads)
    #}}}

    @staticmethod
    #{{{calcAvgAngularVelocityAndSpread
    def calcAvgAngularVelocityAndSpread(angularVelocity):
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
        modeIndices      = range(angularVelocity.shape[1])

        avgAngVels = []
        spreads    = []

        for modeNr in modeIndices:
            # Calculate the average and the spread
            avgAngVels.append(angularVelocity[:,modeNr].mean())
            spreads   .append(angularVelocity[:,modeNr].std() )

        return tuple(avgAngVels), tuple(spreads)
    #}}}

    #{{{getData
    def getData(self             ,\
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
            over the observation "modeNr" over the observation "Scan"
        positionTuple : tuple
            The tuple containing (rho, z).
            Needed in the plotting routine.
        uc : Units Converter
            The units converter used when obtaining the fourier modes.
            Needed in the plotting routine.
        """
        #}}}

        # For ideas on how to append a DataFrame, see:
        # https://github.com/pandas-dev/pandas/blob/master/doc/cheatsheet/Pandas_Cheat_Sheet.pdf
        # http://stackoverflow.com/questions/10715965/add-one-row-in-a-pandas-dataframe

        multiTuples = []
        fullDict =\
            {"growthRate":[], "growthRateStd":[],\
             "averageAngularVelocity":[], "averageAngularVelocityStd":[]}

        loopOver = zip(self._scanCollectPaths,\
                       self._steadyStatePaths,\
                       self._tSlices         ,\
                       )

        # Loop over the folders
        for scanPaths, steadyStatePath, tSlice in loopOver:
            # Obtain the scan value
            scanValue = getScanValue(scanPaths, self._scanParameter)

            # Update with the correct tSlice
            indicesKwargs.update({"tSlice" : tSlice})

            fm, positionTuple, uc = \
                self._collectAndCalcFourierModes(scanPaths        ,\
                                                 varName          ,\
                                                 convertToPhysical,\
                                                 steadyStatePath  ,\
                                                 indicesArgs      ,\
                                                 indicesKwargs    ,\
                                                )

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

            slopes, slopesStd =\
                self.calcSlopeAndSpread(magnitudes[:, modeStart:modeEnd],\
                                        time                            ,\
                                       )

            avgAngVels, avgAngVelsStd =\
                self.calcAvgAngularVelocityAndSpread(\
                                        angularVelocity[:, modeStart:modeEnd])

            for modeInd in range(len(slopes)):
                # Fill the multiIndexTuple and the dict
                multiTuples.append((scanValue, modeInd + 1))

                fullDict["growthRate"               ].\
                        append(slopes[modeInd])
                fullDict["growthRateStd"            ].\
                        append(slopesStd[modeInd])
                fullDict["averageAngularVelocity"   ].\
                        append(avgAngVels[modeInd])
                fullDict["averageAngularVelocityStd"].\
                        append(avgAngVelsStd[modeInd])

        # Make the data frame
        growthRateDataFrame =\
            pd.DataFrame(fullDict,\
                         index=pd.MultiIndex.from_tuples(\
                            multiTuples,
                            names=(self._scanParameter,"modeNr")))

        return growthRateDataFrame, positionTuple, uc
    #}}}

    #{{{_collectAndCalcFourierModes
    def _collectAndCalcFourierModes(self             ,\
                                    scanPaths        ,\
                                    varName          ,\
                                    convertToPhysical,\
                                    steadyStatePath  ,\
                                    indicesArgs      ,\
                                    indicesKwargs    ,\
                                   ):
        #{{{docstring
        """
        Collects and calculates the fourier modes.

        Parameters
        ----------
        scanPaths : tuple
            Tuple of the scan paths.
        varName : str
            Name of variable to find the growth rates of.
        convertToPhysical : bool
            Whether or not to convert to physical units.
        steadyStatePath : str
            String containing the steady state path
        indicesArgs : tuple
            Tuple containing the indices.
        indicesKwargs : dict
            Keyword arguments to use when setting the indices for
            collecting.

        Returns
        -------
        fm : dict
            The dictionary containing the fourier modes.
        positionTuple : tuple
            The tuple containing (rho, z).
            Needed in the plotting routine.
        uc : Units Converter
            The units converter used when obtaining the fourier modes.
            Needed in the plotting routine.
        """
        #}}}

        # Obtain the magnitudes and angular velocities
        ccfm = CollectAndCalcFourierModes(\
                            scanPaths                            ,\
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

        # NOTE: Theta remains unspecified as we have done a fourier transform
        firstKey = tuple(fm.keys())[0]
        rho, z   = firstKey.split(",")

        positionTuple = (rho, z)

        uc = ccfm.uc

        return fm, positionTuple, uc
    #}}}
#}}}

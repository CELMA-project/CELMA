#!/usr/bin/env python

"""
Contains drivers for calculation of the growth rates
"""

from .statsAndSignalsDriver import StatsAndSignalsDrivers
from ..statsAndSignals import PerpPlaneProbes, calcGrowthRate, PlotGrowthRates
from multiprocessing import Process, Pool
from collections import ChainMap
import pandas as pd
import numpy as np

#{{{DriversGrowthRates
class DriversGrowthRates(StatsAndSignalsDrivers):
    """
    Class which handles calculation and plotting of the growth rates
    """

    #{{{Constructor
    def __init__(self                   ,\
                 *args                  ,\
                 var              = None,\
                 scanParam        = None,\
                 yInd             = None,\
                 steadyStatePaths = None,\
                 maxMode          = 7   ,\
                 **kwargs):
        #{{{docstring
        """
        This constructor:

        1. Calls the parent constructor
        2. Sets additional member data.

        NOTE: The order of steadyStatePaths must agree with the order of
              self._paths

        Parameters
        ----------
        *args : positional arguments
            See the constructor of StatsAndSignalsDrivers for details.
        var : str
            Variable to collect
        scanParam : str
            Parameter which is scanned
        yInd : int
            yInd of the simulation data to use
        steadyStatePaths : iterable
            Path to the steady states.
            NOTE! The order of steadyStatePaths must agree with the
                  order of self._paths.
        maxMode : int
            Maximum modes to make a dispersion relation diagram from
        **kwargs : keyword arguments
            See the constructor of StatsAndSignalsDrivers and
FIXME:
            PlotGrowthRates for details.
        """
        #}}}

        # Call the constructor of the parent class
        super().__init__(*args, **kwargs)

        # Set member data
        self._var              = var
        self._scanParam        = scanParam
        self._yInd             = yInd
        self._maxMode          = maxMode
        self._steadyStatePaths = steadyStatePaths

        # Placeholder for the growthRates
        self._growthRatesDict = None
    #}}}

    #{{{calcProbes
    def calcProbes(self):
        """ Calculates the statistics of the probes """

        # Make the growthRatesDict a dictionary
        self._growthRatesDict = {}

        # Create the probes (one for each scan point)
        # Sort the folders (ensures that they come in the right order)
        self._dmp_folders       = list(self._dmp_folders)
        self._steadyStatePaths = list(self._steadyStatePaths)
        self._dmp_folders.sort()
        self._steadyStatePaths.sort()
        self._dmp_folders       = tuple(self._dmp_folders)
        self._steadyStatePaths = tuple(self._steadyStatePaths)

        # Obtain the growth rates dictionary
        pathsAndSteadyStatePathsTuple =\
                tuple(zip(self._paths, self._steadyStatePaths))
        # Appendable list
        listOfDicts = []
        if self._useSubProcess:
            with Pool(len(pathsAndSteadyStatePathsTuple)) as p:
                listOfDicts =\
                    p.map(self._getGrowthRatesDict,\
                          pathsAndSteadyStatePathsTuple)
        else:
            for curTuple in pathsAndSteadyStatePathsTuple:
                listOfDicts.append(\
                        self._getGrowthRatesDict(curTuple))

        # Combine the dictionaries to one
        self._growthRatesDict = dict(ChainMap(*listOfDicts))
    #}}}

    #{{{_getGrowthRatesDict
    def _getGrowthRatesDict(self, pathsAndSteadyStatePathsTuple):
        """
        Obtain the growth rates

        Parameters
        ----------
        pathsAndSteadyStatePathsTuple : tuple
            Tuple of paths and steady state paths

        Returns
        -------
        growthRatesDict : dict
            Dictionary of the growth rates, where the scan parameter is
            the outermost dict, and the modenumber is at the next level
        """

        # Unpack
        path, steadyStatePath = pathsAndSteadyStatePathsTuple
        # We only find the growth rates at position of highest gradient,
        # as we from theory expect the highest growth rates to be there
        curProbe = PerpPlaneProbes(\
                      self._var                                    ,\
                      paths               = path                   ,\
                      yInd                = self._yInd             ,\
                      convertToPhysical   = self.convertToPhysical,\
                      steadyStatePath     = steadyStatePath        ,\
                      nProbes             = 1                      ,\
                      radialProbesIndices = None                   ,\
                     )

        # Initialize the probes in order to calculate
        curProbe.initializeInputOutput(curProbe.radialProbesIndices,\
                                       [curProbe.yInd],\
                                       [0])

        # Calculate the FFT
        curProbe.calcFFTs()

        # Get the positionKey of the probe
        positionKey = tuple(curProbe.results.keys())[0]

        # Obtain the current scan value
        if hasattr(path, "__iter__") and type(path) != str:
            # The path is a sequence
            scanValue = path[0].split("_")
        else:
            # The path is a string
            scanValue = path.split("_")
        # +1 as the value is immediately after the scan parameter
        scanValue = scanValue[scanValue.index(self._scanParam)+1]
        # The time and the time is clipped with three to clip away
        # where initial modes decay
        initClip = 3
        # We will also clip so that approximately only the linear
        # mode is present (less data to process)
        linClip = curProbe.results[positionKey]["zFFTLinearIndex"]
        if linClip <= initClip:
            message = ("{0}{1}WARNING: "\
                       "Could not find a proper startpoint for the "
                       "linear stage{1}{0}")
            linClip = None
            print(message.format("\n"*2, "!"*5))
        # Clipping modes and time
        modes = curProbe.results[positionKey]["zFFT"][initClip:linClip, :]
        time  = curProbe.time[initClip:linClip]

        # Find the growth rates
        growthRatesDict =\
            {"{}={}".format(self._scanParam, scanValue) :\
                calcGrowthRate(modes   = modes,\
                               time    = time ,\
                               maxMode = self._maxMode)\
            }

        return growthRatesDict
    #}}}

    #{{{plotGrowthRates
    def plotGrowthRates(self):
        """Plot the growth rates"""

        # Calculate the probes if not already done
        if self._growthRatesDict == None:
            self.calcProbes()

        #{{{ Cast to dataframe in order to plot more easily
        # We start by cleaning up the dict by removing "allGrowthRates"
        # and replace "None" with NaN
        # Loop over the scan
        for scan in self._growthRatesDict.keys():
            # Loop over the mode numbers
            for mNr in self._growthRatesDict[scan].keys():
                # Replace None with NaN
                if self._growthRatesDict[scan][mNr] is None:
                    message = ("{0}{1}WARNING: No hits in mode nr {2} in {3}. "
                               "Setting to NaN{1}{0}")
                    print(message.format("\n", "!"*4, mNr, scan))
                    self._growthRatesDict[scan][mNr] =\
                        {"growthRate"   :np.nan,\
                         "growthRateStd":np.nan,\
                         "angFreq"      :np.nan,\
                         "angFreqStd"   :np.nan}
                else:
                    # Remove "allGrowthRates" by pop (pop returns the
                    # key)
                    # The second argument specifies what should be
                    # returned if the key is not found. If unspecified
                    # as keyError will be thrown
                    self._growthRatesDict[scan][mNr].pop("allGrowthRates",None)

        # We would now like to write the mode numbers as modeNr=mNr,
        # rather than just mNr
        # In order to avoid memory problems, it will be done in the
        # following way
        firstDict = tuple(self._growthRatesDict.keys())
        secondDict = tuple(self._growthRatesDict[firstDict[0]].keys())
        for f in firstDict:
            for s in secondDict:
                self._growthRatesDict[f]["modeNr={}".format(s)] =\
                    self._growthRatesDict[f].pop(s)

        # Reform a dict of dict of dict so that we get a dict with tuple keys
        # http://stackoverflow.com/questions/24988131/nested-dictionary-to-multiindex-dataframe-where-dictionary-keys-are-column-label
        reform = {(outerKey, innerKey): values\
                  for outerKey, innerDict in self._growthRatesDict.items()\
                  for innerKey, values in innerDict.items()}

        # This can now be made to a data frame
        growthRatesDF = pd.DataFrame.from_dict(reform)

        # NOTE: The format is now
        #       cols = ("Scan", "Mode") and the indices are the
        #       different angFreq, growthRates etc (hereby refered to as
        #       the "Data")
        # NOTE: We would rather like
        #       Columns to be either "Scan" or "mNr" (depending on what
        #       we would like to plot) and make the rest indices When
        #       this is done we can easily slice the data.
        # We start by transposing the indices and columns
        growthRatesDF=growthRatesDF.T
        # Then we add names to the frame (easier unstack and stack)
        growthRatesDF.index.names = ("Scan", "Mode")
        growthRatesDF.columns.name="Data"
        #}}}

        # Finally we make either "Mode" or "Scan" the column and the "Data"
        # the index
        # http://nikgrozev.com/2015/07/01/reshaping-in-pandas-pivot-pivot-table-stack-and-unstack-explained-with-pictures/
        # This can be done by stack: Rotating column index to row index
        # And unstack              : Rotating row index to column index
        if self._useSubProcess:
            #{{{ Function call through subprocess
            # We would like to have "Scan" on the x axis
            df = growthRatesDF.unstack("Scan").stack("Data")
            gRPlotterScan = PlotGrowthRates(\
                                self._paths                                ,\
                                growthRates       = df                     ,\
                                convertToPhysical = self.convertToPhysical,\
                                showPlot          = self._showPlot         ,\
                                savePlot          = self._savePlot         ,\
                                extension         = self._extension        ,\
                                savePath          = self._savePath         ,\
                                pltSize           = self._pltSize          ,\
                                )
            # NOTE: If we had growthRatesDF.unstack("Scan").stack("Data")
            #       and would've liked to swap "Mode" and "Scan" with the
            #       current dataframe, we would have to
            #
            #       growthRatesDF =
            #       growthRatesDF.unstack("Mode").stack("Scan").swaplevel()
            #
            #       The swap needed to make the "Data" the innermost level
            # We would like to thave "Mode" on the x axis
            df = growthRatesDF.unstack("Mode").stack("Data")
            gRPlotterMNr = PlotGrowthRates(\
                                self._paths                                ,\
                                growthRates       = df                     ,\
                                convertToPhysical = self.convertToPhysical,\
                                showPlot          = self._showPlot         ,\
                                savePlot          = self._savePlot         ,\
                                extension         = self._extension        ,\
                                savePath          = self._savePath         ,\
                                pltSize           = self._pltSize          ,\
                                )

            proc = {}
            # Create process
            proc["scanOnXAxis"] =\
                Process(\
                    target = gRPlotterScan.plotGrowthRates(),\
                    args   = ()                             ,\
                    kwargs = {})
            proc["modeOnXAxis"] =\
                Process(\
                    target = gRPlotterMNr.plotGrowthRates(),\
                    args   = ()                            ,\
                    kwargs = {})

            # Start processes
            for key in proc.keys():
                proc[key].start()
            # Wait for process to finish
            for key in proc.keys():
                proc[key].join()
            #}}}
        else:
            #{{{Normal function call
            # We would like to have "Scan" on the x axis
            df = growthRatesDF.unstack("Scan").stack("Data")
            gRPlotter = PlotGrowthRates(\
                            self._paths                                ,\
                            growthRates       = df                     ,\
                            convertToPhysical = self.convertToPhysical,\
                            showPlot          = self._showPlot         ,\
                            savePlot          = self._savePlot         ,\
                            extension         = self._extension        ,\
                            savePath          = self._savePath         ,\
                            pltSize           = self._pltSize          ,\
                            )

            gRPlotter.plotGrowthRates()
            # NOTE: If we had growthRatesDF.unstack("Scan").stack("Data")
            #       and would've liked to swap "Mode" and "Scan" with the
            #       current dataframe, we would have to
            #
            #       growthRatesDF =
            #       growthRatesDF.unstack("Mode").stack("Scan").swaplevel()
            #
            #       The swap needed to make the "Data" the innermost level
            # We would like to thave "Mode" on the x axis
            df = growthRatesDF.unstack("Mode").stack("Data")
            gRPlotter = PlotGrowthRates(\
                            self._paths                                ,\
                            growthRates       = df                     ,\
                            convertToPhysical = self.convertToPhysical,\
                            showPlot          = self._showPlot         ,\
                            savePlot          = self._savePlot         ,\
                            extension         = self._extension        ,\
                            savePath          = self._savePath         ,\
                            pltSize           = self._pltSize          ,\
                            )
            gRPlotter.plotGrowthRates()
            #}}}
    #}}}
#}}}

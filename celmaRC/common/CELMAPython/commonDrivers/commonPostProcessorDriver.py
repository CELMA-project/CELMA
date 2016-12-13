#!/usr/bin/env python

"""
Contains the common post processing driver class
"""

from .saveFolderFuncs import scanWTagSaveFunc
from matplotlib.pylab import plt
import re
import datetime
import os

#{{{CommonPostProcessingDriver
class CommonPostProcessingDriver(object):
    """
    The parent driver of all BOUT++ post processing functions

    * Sets the memberdata.
    * Sets the save path.
    """

    #{{{Constructor
    def __init__(self                     ,\
                 dmp_folder               ,\
                 convertToPhysical = False,\
                 fluctuation       = False,\
                 showPlot          = False,\
                 savePlot          = True ,\
                 saveFolder        = None ,\
                 saveFolderFunc    = None ,\
                 useSubProcess     = True ,\
                 extension         = "png",\
                 scanParameters    = None ,\
                 **kwargs):
        #{{{docstring
        """
        This constructor sets the memberdata, and sets the savePath.

        Parameters
        ----------
        dmp_folder : [str|sequence of str]
            The path to collect from. This is the output dmp_folder from
            bout_runners.
        convertToPhysical : bool
            If the physical or normalized units should be plotted.
        fluctuation : bool
            If the poloidal average should be subtracted from the data
        showPlot : bool
            If the plot should be displayed.
        savePlot : bool
            If plot should be saved.
        saveFolder : str
            Name of the folder to save the plots in.
        saveFolderFunc : str
            Name of an implemented function which returns the name of
            the folder to save plots.
        useSubProcess : bool
            Whether each plot will be made by a new sub process, or the
            plots should be made in series.
        extension : str
            The extension to use when saving non-animated plots
        scanParameters : [None|sequence (not string)]
            Sequence of parameters changed in the scan. If this is not None,
            calls to convertToCurrentScanParameters will be triggered in
            the child classes.
        **kwargs : keyword arguments
            Additional keyword arguments given as input to saveFolderFunc.
        """
        #}}}

        # Set the member data
        self._dmp_folders       = dmp_folder
        self._convertToPhysical = convertToPhysical
        self._fluctuation       = fluctuation
        self._showPlot          = showPlot
        self._savePlot          = savePlot
        self._saveFolder        = saveFolder
        self._useSubProcess     = useSubProcess
        self._extension         = extension
        self._scanParameters    = scanParameters

        #{{{Set the saveFolder
        if saveFolderFunc is not None:
            # FIXME: Check if it is possible to change the API here. Would
            # be nice if could send in a function instead
            if saveFolderFunc == 'scanWTagSaveFunc':
                saveFolder = scanWTagSaveFunc(dmp_folder                       ,\
                                    convertToPhysical = convertToPhysical,\
                                    **kwargs)
            else:
                message  = "{0}Warning: saveFolderFunc '{1}' not found, "
                message += "falling back to standard implementation{0}"
                print(message.format("\n"*3, saveFolderFunc))
                saveFolder = "-".join(dmp_folder.split('/')[::-1])
        else:
            if saveFolder is None:
                saveFolder = "-".join(dmp_folder.split('/')[::-1])

        # Get the timefolder
        self._timeFolder = self._getTime()


        # Differentiate whether or not there are several folders
        if hasattr(dmp_folder, "__iter__") and type(dmp_folder) != str:
            dmp_folder = dmp_folder[0]

        # Create the savepath (based on the first dmp_folder string)
        visualizationType = "Physical" if convertToPhysical else "Normalized"
        saveDirs = [os.path.normpath(dmp_folder).split(os.sep)[0],\
                    "visualization{}".format(visualizationType),\
                    saveFolder,\
                    self._timeFolder]
        if self._fluctuation:
            saveDirs.append("fluctuation")
        self._savePath = ""
        for saveDir in saveDirs:
            self._savePath = os.path.join(self._savePath, saveDir)

        # Make dir if not exists
        if not os.path.exists(self._savePath):
            os.makedirs(self._savePath)
        #}}}

        if self._useSubProcess:
            #{{{ The multiprocess currently only works with the Agg backend
            # Qt4Agg currently throws
            # [xcb] Unknown sequence number while processing queue
            # [xcb] Most likely this is a multi-threaded client and XInitThreads has not been called
            # [xcb] Aborting, sorry about that.
            # python: ../../src/xcb_io.c:274: poll_for_event: Assertion
            # `!xcb_xlib_threads_sequence_lost' failed.
            #}}}
            plt.switch_backend('Agg')
    #}}}

    #{{{_getTime
    def _getTime(self, depth = 'second'):
        """
        Gets the current time, and returns it as a string

        Parameters
        ----------
        depth : ['hour' | 'minute' | 'second']
            String giving the temporal accuracy of the output string.

        Returns
        -------
        nowStr : str
            The string containing the current time
        """
        now = datetime.datetime.now()
        nowStr = "{}-{:02d}-{:02d}".format(now.year, now.month, now.day)

        if depth == 'hour' or depth == 'minute' or depth == 'second':
            nowStr += "-{:02d}".format(now.hour)
        if depth == 'minute' or depth == 'second':
            nowStr += "-{:02d}".format(now.minute)
        if depth == 'second':
            nowStr += "-{:02d}".format(now.second)

        return nowStr
    #}}}

    #{{{_toggleSubPolAvg
    def _toggleSubPolAvg(self, value):
        """
        Sets fluctuation and modifies savePath

        Parameters
        ----------
        value : bool
            Value to set to fluctuation
        """

        if self._fluctuation != value:
            if self._fluctuation == True and value == False:
                # Remove fluctuation from the path
                paths = tuple(os.path.split(self._savePath))
                paths.remove("fluctuation")
                self._savePath = os.path.join(*paths)
                self._fluctuation = value
            elif self._fluctuation == False and value == True:
                self._savePath = os.path.join(self._savePath, "fluctuation")
                self._fluctuation = value

            # Make dir if not exists
            if not os.path.exists(self._savePath):
                os.makedirs(self._savePath)
    #}}}

    #{{{_convertToCurrentScanParameters
    def _convertToCurrentScanParameters(self, aScanPath):
        #{{{docstring
        """
        Function which converts a path belonging to one paths in a scan
        to the path belonging to the current scan.

        The function obtains the current scan parameters from
        self._dmp_folders (should be tuple with one element), and
        inserts the current scan parameters into aScanPath (the function
        input which is one of the paths belonging to the scan).

        Parameters
        ----------
        aScanPath : str
            One of the paths from the simulations.

        Returns
        -------
        scanPath : str
            aScanPath converted to the scan parameters of the current run.
        """
        #}}}

        # Make a template string of aScanPath
        scanPathTemplate = aScanPath
        for scanParameter in self._scanParameters:
            hits = [m.start() for m in \
                    re.finditer(scanParameter, scanPathTemplate)]
            while(len(hits) > 0):
                # Replace the values with {}
                # The value is separated from the value by 1 character
                value_start = hits[0] + len(scanParameter) + 1
                # Here we assume that the value is not separated by an
                # underscore
                value_len = len(scanPathTemplate[value_start:].split("_")[0])
                value_end = value_start + value_len
                # Replace the values with {}
                scanPathTemplate =\
                    "{}{{0[{}]}}{}".format(\
                        scanPathTemplate[:value_start],\
                        scanParameter,\
                        scanPathTemplate[value_end:])
                # Update hits
                hits.remove(hits[0])

        # Get the values from the current self._dmp_folders
        values = {}
        for scanParameter in self._scanParameters:
            hits = [m.start() for m in \
                    re.finditer(scanParameter, self._dmp_folders[0])]
            # Choose the first hit to get the value from (again we assume
            # that the value does not contain a _)
            value_start = hits[0] + len(scanParameter) + 1
            # Here we assume that the value is not separated by an
            # underscore
            values[scanParameter]=self._dmp_folders[0][value_start:].split("_")[0]

        # Insert the values
        scanPath = scanPathTemplate.format(values)

        return scanPath
    #}}}
#}}}

#!/usr/bin/env python

"""
Contains the common post processing driver class
"""

from ..driverHelpers import scanWTagSaveFunc
import matplotlib.pyplot as plt
import datetime
import os

#{{{DriverPostProcessingSuperClass
class DriverPostProcessingSuperClass(object):
    """
    The parent driver of all BOUT++ post processing functions

    * Sets common member data.
    * Sets the save path.
    """

    #{{{Constructor
    def __init__(self                     ,\
                 dmp_folders              ,\
                 collectPaths      = None ,\
                 convertToPhysical = False,\
                 showPlot          = False,\
                 savePlot          = True ,\
                 saveFolder        = None ,\
                 saveFolderFunc    = None ,\
                 useSubProcess     = True ,\
                 extension         = "png",\
                 **kwargs):
        #{{{docstring
        """
        This constructor sets the memberdata, and sets the savePath.

        Parameters
        ----------
        dmp_folders: tuple
            This is the output dmp_folder from bout_runners.
            Typically, these are the folders in a given scan
        collectPaths : [None|tuple]
            Paths to collect from.
            If None dmp_folders will be set to collectPaths
        convertToPhysical : bool
            If the physical or normalized units should be plotted.
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
        if collectPaths is None:
            self._collectPaths  = tuple(dmp_folders)
        else:
            self._collectPaths = collectPaths

        self._convertToPhysical = convertToPhysical
        self._showPlot          = showPlot
        self._savePlot          = savePlot
        self._saveFolder        = saveFolder
        self._useSubProcess     = useSubProcess
        self._extension         = extension

        # Get the firs dmp_folder (will be used to create the savepath)
        dmp_folder = dmp_folders[0]

        # Set the saveFolder
        if saveFolderFunc is not None:
            # FIXME: Check if it is possible to change the API here.
            #        Would be nice if could send in a function instead.
            if saveFolderFunc == "scanWTagSaveFunc":
                saveFolder = scanWTagSaveFunc(dmp_folder, **kwargs)
            else:
                message  = "{0}Warning: saveFolderFunc '{1}' not found, "
                message += "falling back to standard implementation{0}"
                print(message.format("\n"*3, saveFolderFunc))
                saveFolder = "-".join(dmp_folder.split("/")[::-1])
        else:
            if saveFolder is None:
                saveFolder = "-".join(dmp_folder.split("/")[::-1])

        # Get the timefolder
        self._timeFolder = self._getTime()

        # Create the savepath (based on the first dmp_folder string)
        visualizationType = "Physical" if convertToPhysical else "Normalized"
        saveDirs = [os.path.normpath(dmp_folder).split(os.sep)[0],\
                    "visualization{}".format(visualizationType),\
                    saveFolder,\
                    self._timeFolder]

        self._savePath = ""
        for saveDir in saveDirs:
            self._savePath = os.path.join(self._savePath, saveDir)

        if self._useSubProcess:
            #{{{ The multiprocess currently only works with the Agg backend
            # Qt4Agg currently throws
            # [xcb] Unknown sequence number while processing queue
            # [xcb] Most likely this is a multi-threaded client and XInitThreads has not been called
            # [xcb] Aborting, sorry about that.
            # python: ../../src/xcb_io.c:274: poll_for_event: Assertion
            # `!xcb_xlib_threads_sequence_lost' failed.
            #}}}
            plt.switch_backend("Agg")
    #}}}

    #{{{_getTime
    def _getTime(self, depth = "second"):
        """
        Gets the current time, and returns it as a string

        Parameters
        ----------
        depth : ["hour" | "minute" | "second"]
            String giving the temporal accuracy of the output string.

        Returns
        -------
        nowStr : str
            The string containing the current time
        """
        now = datetime.datetime.now()
        nowStr = "{}-{:02d}-{:02d}".format(now.year, now.month, now.day)

        if depth == "hour" or depth == "minute" or depth == "second":
            nowStr += "-{:02d}".format(now.hour)
        if depth == "minute" or depth == "second":
            nowStr += "-{:02d}".format(now.minute)
        if depth == "second":
            nowStr += "-{:02d}".format(now.second)

        return nowStr
    #}}}
#}}}

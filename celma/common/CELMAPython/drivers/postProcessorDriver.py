#!/usr/bin/env python

"""
Contains the post processing driver class
"""

from matplotlib.pylab import plt
import datetime
import os

#{{{PostProcessorDriver
class PostProcessorDriver(object):
    """
    The parent driver of all BOUT++ post processing functions

    Sets the memberdata
    """

    #{{{Constructor
    def __init__(self                             ,\
                 path                             ,\
                 xguards           = False        ,\
                 yguards           = False        ,\
                 xSlice            = slice(0,None),\
                 ySlice            = slice(0,None),\
                 zSlice            = slice(0,None),\
                 tSlice            = None         ,\
                 convertToPhysical = False        ,\
                 subPolAvg         = False        ,\
                 showPlot          = False        ,\
                 savePlot          = True         ,\
                 saveFolder        = None         ,\
                 saveFolderFunc    = None         ,\
                 useSubProcess     = True         ,\
                 **kwargs
                ):
        #{{{docstring
        """
        This constructor sets the memberdata, and sets the savePath.

        Parameters
        ----------
        path : str
            The path to collect from.
        xguards : bool
            If xguards should be included when collecting.
        yguards : bool
            If yguards should be included when collecting.
        xSlice : slice
            How the data will be sliced in x.
        ySlice : slice
            How the data will be sliced in y.
        zSlice : slice
            How the data will be sliced in z.
        tSlice : slice
            How the data will be sliced in t.
        convertToPhysical : bool
            If the physical or normalized units should be plotted.
        subPolAvg : bool
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
        **kwargs : keyword arguments
            Additional keyword arguments given as input to saveFolderFunc.
        """
        #}}}

        # Set the member data
        self._path              = path
        self._xguards           = xguards
        self._yguards           = yguards
        self._xSlice            = xSlice
        self._ySlice            = ySlice
        self._zSlice            = zSlice
        self._tSlice            = tSlice
        self._convertToPhysical = convertToPhysical
        self._subPolAvg         = subPolAvg
        self._showPlot          = showPlot
        self._savePlot          = savePlot
        self._saveFolder        = saveFolder
        self._useSubProcess     = useSubProcess

        #{{{Set the saveFolder
        if saveFolderFunc is not None:
            # FIXME: Check if it is possible to change the API here. Would
            # be nice if could send in a function instead
            if saveFolder == 'scanWTagSaveFunc':
                saveFolderFunc = scanWTagSaveFunc(path                   ,\
                                                  xguards    = xguards   ,\
                                                  yguards    = yguards   ,\
                                                  xSlice     = xSlice    ,\
                                                  ySlice     = ySlice    ,\
                                                  zSlice     = zSlice    ,\
                                                  tSlice     = tSlice    ,\
                                                  saveFolder = saveFolder,\
                                                  **kwargs)
            else:
                message  = "{0}Warning: saveFolderFunc '{1}' not found, "
                message += "falling back to standard implementation{0}"
                print(message.format("\n"*3, saveFolderFunc))
                saveFolder = "-".join(path.split('/')[::-1])
        else:
            if saveFolder is None:
                saveFolder = "-".join(path.split('/')[::-1])

        # Get the timefolder
        self._timeFolder = self._getTime()

        # Create the savepath
        saveDirs = [os.path.dirname(self._path),\
                    'visualization',\
                    saveFolder,\
                    self._timeFolder]
        if self._subPolAvg:
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

    #{{{_setSubPolAvg
    def _setSubPolAvg(self, value):
        """
        Sets subPolAvg and modifies savePath

        Parameters
        ----------
        value : bool
            Value to set to subPolAvg
        """

        if self._subPolAvg != value:
            if self._subPolAvg == True and value == False:
                # Remove fluctuation from the path
                paths = list(os.path.split(self._savePath))
                paths.remove("fluctuation")
                self._savePath = os.path.join(*paths)
                self._subPolAvg = value
            elif self._subPolAvg == False and value == True:
                self._savePath = os.path.join(self._savePath, "fluctuation")
                self._subPolAvg = value

            # Make dir if not exists
            if not os.path.exists(self._savePath):
                os.makedirs(self._savePath)
    #}}}
#}}}

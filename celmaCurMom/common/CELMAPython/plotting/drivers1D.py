#!/usr/bin/env python

"""
Contains drivers for plotting 1D plots
"""

from .plotters import Plot1D
from .getStrings import getSaveString
from .lineGetters import (getMainFields,\
                          getLnNFields,\
                          getJParFields,\
                          momDensParFields,\
                          getVortDFields,\
                         )
import matplotlib.pyplot as plt
import numpy as np
from multiprocessing import Process

#{{{Drivers
#{{{parPerpDriver
def parPerpDriver(path                ,\
                  xSlice        = 0   ,\
                  ySlice        = 0   ,\
                  zSlice        = 0   ,\
                  saveFolder    = None,\
                  useSubProcess = True,\
                  timeFolder    = None,\
                  **kwargs):
    #{{{docstring
    """
    Wrapper function for the parallel and perpendicular plots.

    Specific parPerpDriver input
    useSubProcess - Each plot will be made by a new sub process
    timeFolder    - The name of the timefolder
                    Enables several plots to be put into same timeFolder

    Specific parPerpDriver output:
    timeFolder - The name of the timefolder
                 Enables several plots to be put into same timeFolder

    See docstring of single1DDriver for more info.
    """
    #}}}

    if saveFolder == None:
        saveFolder = 'parPerp'

    if useSubProcess:
        #{{{ The multiprocess currently only works with the Agg backend
        # Qt4Agg currently throws
        # [xcb] Unknown sequence number while processing queue
        # [xcb] Most likely this is a multi-threaded client and XInitThreads has not been called
        # [xcb] Aborting, sorry about that.
        # python: ../../src/xcb_io.c:274: poll_for_event: Assertion
        # `!xcb_xlib_threads_sequence_lost' failed.
        #}}}
        plt.switch_backend('Agg')
        #{{{ Function call through subprocess
        if timeFolder is None:
            # Create the savepath if not already set
            prePaths = ['visualization', saveFolder]
            saveString, timeFolder = getSaveString(''                     ,\
                                                   path                   ,\
                                                   timeFolder = timeFolder,\
                                                   prePaths   = prePaths  ,\
                                                   )

        p1 = Process(\
                    target = perpDriver                    ,\
                    args   = (path,)                       ,\
                    kwargs = {'ySlice'       :ySlice       ,\
                              'zSlice'       :zSlice       ,\
                              'saveFolder'   :saveFolder   ,\
                              'timeFolder'   :timeFolder   ,\
                              'useSubProcess':useSubProcess,\
                              **kwargs
                             }
                    )

        p2 = Process(\
                    target = parDriver                     ,\
                    args   = (path,)                       ,\
                    kwargs = {'xSlice'       :xSlice       ,\
                              'zSlice'       :zSlice       ,\
                              'saveFolder'   :saveFolder   ,\
                              'timeFolder'   :timeFolder   ,\
                              'useSubProcess':useSubProcess,\
                              **kwargs
                             }
                    )

        p1.start()
        p2.start()
        #}}}
    else:
        #{{{ Normal function call
        timeFolder = perpDriver(path                         ,\
                                ySlice        = ySlice       ,\
                                zSlice        = zSlice       ,\
                                saveFolder    = saveFolder   ,\
                                timeFolder    = timeFolder   ,\
                                useSubProcess = useSubProcess,\
                                **kwargs)

        timeFolder = parDriver(path                         ,\
                               xSlice        = xSlice       ,\
                               zSlice        = zSlice       ,\
                               saveFolder    = saveFolder   ,\
                               timeFolder    = timeFolder   ,\
                               useSubProcess = useSubProcess,\
                               **kwargs)
        #}}}

    return timeFolder
#}}}

#{{{perpDriver
def perpDriver(path                ,\
               ySlice        = 0   ,\
               zSlice        = 0   ,\
               saveFolder    = None,\
               skipPlots     = None,\
               useSubProcess = True,\
               timeFolder    = None,\
               **kwargs):
    #{{{docstring
    """
    Wrapper function for the perpendicular plots.

    Specific perpDriver input
    skipPlots     - List of plots to skip
    useSubProcess - Each plot will be made by a new sub process
    timeFolder    - The name of the timefolder
                    Enables several plots to be put into same timeFolder

    Specific perpDriver output:
    timeFolder - The name of the timefolder
                 Enables several plots to be put into same timeFolder

    See docstring of single1DDriver for more info.
    """
    #}}}

    if saveFolder == None:
        saveFolder = 'perpendicular'

    # Initiate the pltNames and timeFolder
    pltNames   = ["mainFields" ,\
                  "lnNFields"  ,\
                  "jParFields",\
                  "momDensParFields",\
                  "vortDFields",\
                  ]

    # Take out the skip plots
    if skipPlots is not None:
        for skip in skipPlots:
            pltNames.remove(skip)

    if useSubProcess:
        #{{{ The multiprocess currently only works with the Agg backend
        # Qt4Agg currently throws
        # [xcb] Unknown sequence number while processing queue
        # [xcb] Most likely this is a multi-threaded client and XInitThreads has not been called
        # [xcb] Aborting, sorry about that.
        # python: ../../src/xcb_io.c:274: poll_for_event: Assertion
        # `!xcb_xlib_threads_sequence_lost' failed.
        #}}}
        plt.switch_backend('Agg')
        #{{{ Function call through subprocess
        if timeFolder is None:
            # Create the savepath if not already set
            prePaths = ['visualization', saveFolder]
            saveString, timeFolder = getSaveString(''                     ,\
                                                   path                   ,\
                                                   timeFolder = timeFolder,\
                                                   prePaths   = prePaths  ,\
                                                   )

        # Make an array of processes
        processes = []
        for pltName in pltNames:
            processes.append(\
                Process(\
                        target = single1DDriver          ,\
                        args   = (path,)                 ,\
                        kwargs = {'ySlice'    :ySlice    ,\
                                  'zSlice'    :zSlice    ,\
                                  'saveFolder':saveFolder,\
                                  'timeFolder':timeFolder,\
                                  'pltName'   :pltName   ,\
                                  **kwargs
                                 }
                       )
            )

        for process in processes:
            process.start()
        #}}}
    else:
        #{{{Normal function call
        for pltName in pltNames:
            timeFolder = single1DDriver(path                   ,\
                                        ySlice     = ySlice    ,\
                                        zSlice     = zSlice    ,\
                                        saveFolder = saveFolder,\
                                        timeFolder = timeFolder,\
                                        pltName    = pltName   ,\
                                        **kwargs)
        #}}}
    return timeFolder
#}}}

#{{{parDriver
def parDriver(path                ,\
              xSlice        = 0   ,\
              zSlice        = 0   ,\
              saveFolder    = None,\
              skipPlots     = None,\
              useSubProcess = True,\
              timeFolder    = None,\
              **kwargs):
    #{{{docstring
    """
    Wrapper function for the parallel plots.

    Specific parDriver input
    skipPlots     - List of plots to skip
    useSubProcess - Each plot will be made by a new sub process
    timeFolder    - The name of the timefolder
                    Enables several plots to be put into same timeFolder

    Specific parDriver output:
    timeFolder - The name of the timefolder
                 Enables several plots to be put into same timeFolder

    See docstring of single1DDriver for more info.
    """
    #}}}

    if saveFolder == None:
        saveFolder = 'parallel'

    # Initiate the pltNames and timeFolder
    pltNames   = ["mainFields" ,\
                  "lnNFields"  ,\
                  "jParFields",\
                  "momDensParFields",\
                  "vortDFields",\
                  ]

    # Take out the skip plots
    if skipPlots is not None:
        for skip in skipPlots:
            pltNames.remove(skip)

    if useSubProcess:
        #{{{ The multiprocess currently only works with the Agg backend
        # Qt4Agg currently throws
        # [xcb] Unknown sequence number while processing queue
        # [xcb] Most likely this is a multi-threaded client and XInitThreads has not been called
        # [xcb] Aborting, sorry about that.
        # python: ../../src/xcb_io.c:274: poll_for_event: Assertion
        # `!xcb_xlib_threads_sequence_lost' failed.
        #}}}
        plt.switch_backend('Agg')
        #{{{ Function call through subprocess
        if timeFolder is None:
            # Create the savepath if not already set
            prePaths = ['visualization', saveFolder]
            saveString, timeFolder = getSaveString(''                     ,\
                                                   path                   ,\
                                                   timeFolder = timeFolder,\
                                                   prePaths   = prePaths  ,\
                                                   )

        # Make an array of processes
        processes = []
        for pltName in pltNames:
            processes.append(\
                Process(\
                        target = single1DDriver          ,\
                        args   = (path,)                 ,\
                        kwargs = {'xSlice'    :xSlice    ,\
                                  'zSlice'    :zSlice    ,\
                                  'saveFolder':saveFolder,\
                                  'timeFolder':timeFolder,\
                                  'pltName'   :pltName   ,\
                                  **kwargs
                                 }
                       )
            )

        for process in processes:
            process.start()
        #}}}
    else:
        #{{{ Normal function call
        for pltName in pltNames:
            timeFolder = single1DDriver(path                   ,\
                                        xSlice     = xSlice    ,\
                                        zSlice     = zSlice    ,\
                                        saveFolder = saveFolder,\
                                        timeFolder = timeFolder,\
                                        pltName    = pltName   ,\
                                        **kwargs)
        #}}}

    return timeFolder
#}}}

#{{{single1DDriver
def single1DDriver(path                      ,\
                   xguards    = False        ,\
                   yguards    = False        ,\
                   marker     = 'o'          ,\
                   xSlice     = slice(0,None),\
                   ySlice     = slice(0,None),\
                   zSlice     = slice(0,None),\
                   tSlice     = None         ,\
                   subPolAvg  = False        ,\
                   physicalU  = False        ,\
                   showPlot   = False        ,\
                   savePlot   = True         ,\
                   saveFolder = None         ,\
                   pltName    = None         ,\
                   timeFolder = None         ,\
                   ):
    #{{{docstring
    """
    Driver for a single plot.

    Input:
    path       - The simulation folder
    xguards    - If the xguards are to be collected
    yguards    - If the yguards are to be collected
    marker     - Marker in the plots
    xSlice     - How to slice in x
    ySlice     - How to slice in y
    zSlice     - How to slice in z
    tSlice     - How to slice in t
    subPolAvg  - Whether or not the poloidal average should be
                 subtracted from the data
    physicalU  - If the physical units should be plotted
    showPlot   - If the plot is to be displayed
    savePlot   - If the plot is to be saved
    saveFolder - Name of save folder
    pltName    - Name of plot to make
    timeFolder - The name of the timefolder
                 Enables several plots to be put into same timeFolder

    Output:
    timeFolder - The name of the timefolder
                 Enables several plots to be put into same timeFolder
    """
    #}}}

    # Make the plotter object
    plotter = Plot1D(path                   ,\
                     xguards    = xguards   ,\
                     yguards    = yguards   ,\
                     marker     = marker    ,\
                     xSlice     = xSlice    ,\
                     ySlice     = ySlice    ,\
                     zSlice     = zSlice    ,\
                     tSlice     = tSlice    ,\
                     physicalU  = physicalU ,\
                     subPolAvg  = subPolAvg ,\
                     showPlot   = showPlot  ,\
                     savePlot   = savePlot  ,\
                     saveFolder = saveFolder,\
                    )

    # Select plt type
    if pltName == None:
        pltName = saveFolder

    if pltName == "mainFields":
        orgObj = getMainFields(path)
    elif pltName == "lnNFields":
        orgObj = getLnNFields(path)
    elif pltName == "jParFields":
        orgObj = getJParFields(path)
    elif pltName == "momDensParFields":
        orgObj = momDensParFields(path)
    elif pltName == "vortDFields":
        orgObj = getVortDFields(path)
    else:
        raise NotImplementedError("{} is not implemented".format(pltName))

    # Prepare the lines for plotting
    fig = orgObj.pltPrepare()

    # Collect with the plotter object
    for line in orgObj.lines:
        plotter.collectLine(line)

    if orgObj.useCombinedPlot:
        orgObj.makeCombinedLine()

    # Treatment of extra lines
    if pltName == "mainFields":
        # Find lnN, uEPar and uIPar
        for line in orgObj.lines:
            if line.name == 'lnN':
                lnN = line.field
            if line.name == 'uEPar':
                uEPar = line.field
            if line.name == 'uIPar':
                uIPar = line.field
        # Create the n line (second extra line)
        orgObj.extraLines['n'].field = np.exp(lnN)
        if orgObj.extraLines['n'].plotPos:
            orgObj.lines.insert(orgObj.extraLines['n'].plotPos,\
                                orgObj.extraLines['n'])
        else:
            orgObj.lines.append(orgObj.extraLines['n'])

    # Get the correct units and numbers
    for line in orgObj.lines:
        line.field =\
                plotter._getUnitsAndSetPhysical(line.name, line.field)

        if plotter.physicalU:
            if line.name == "momDensPar":
                line.label = "$m_i$" + line.label
            if pltName == "mainFields":
                line.label += r" $[{}]$".format(plotter._units)

    # Do the plot
    timeFolder = plotter.plotDriver(fig, orgObj, timeFolder=timeFolder)

    return timeFolder
#}}}
#}}}

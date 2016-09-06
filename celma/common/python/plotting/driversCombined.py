#!/usr/bin/env python

"""
Contains drivers for plotting 1D and 2D plots at the same time
"""

from .drivers1D import parPerpDriver
from .drivers2D import allMainFields2DDriver
from .getStrings import getSaveString
from .saveFolderFuncs import scanWTagSaveFunc
import matplotlib.pyplot as plt
from multiprocessing import Process

#{{{combinedDriver
def combinedDriver(path                          ,\
                   xguards        = False        ,\
                   yguards        = False        ,\
                   marker         = 'o'          ,\
                   xSlice         = slice(0,None),\
                   ySlice         = slice(0,None),\
                   zSlice         = slice(0,None),\
                   tSlice         = None         ,\
                   showPlot       = False        ,\
                   savePlot       = True         ,\
                   saveFolder     = None         ,\
                   saveFolderFunc = None         ,\
                   varMax         = None         ,\
                   varMin         = None         ,\
                   varyMaxMin     = None         ,\
                   useSubProcess  = True         ,\
                   timeFolder     = None         ,\
                   **kwargs):
    #{{{docstring
    """
    A wrapper function for combined combined1D2D

    Input:
    path           - The simulation folder
    xguards        - If the xguards are to be collected
    yguards        - If the yguards are to be collected
    marker         - Marker in the plots
    xSlice         - How to slice in x
    ySlice         - How to slice in y
    zSlice         - How to slice in z
    tSlice         - How to slice in t
    showPlot       - If the plot is to be displayed
    savePlot       - If the plot is to be saved
    saveFolder     - Name of save folder
    saveFolderFunc - Function which takes this input as an input and
                     returns the saveFolder. To simplify the calling
                     sequence when running on a server, this is given
                     as a string.
    varMax         - Setting a hard upper limit z-axis in the plot
    varMin         - Setting a hard lower limit z-axis in the plot
    varyMaxMin     - Whether or not the limits of the z-axis should be set
                     to the max/min of the current timestep or not
    useSubProcess  - Each plot will be made by a new sub process
    timeFolder     - The name of the timefolder
                     Enables several plots to be put into same timeFolder
    **kwargs       - Extra keyword arguments used in input functions

    Output:
    timeFolder - The name of the timefolder
                 Enables several plots to be put into same timeFolder
    """
    #}}}

    if saveFolderFunc is not None:
        if saveFolderFunc == 'scanWTagSaveFunc':
            saveFolder = scanWTagSaveFunc(path                   ,\
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
    else:
        if saveFolder is None:
            saveFolder = "-".join(path.split('/')[::-1])

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

        # Plots without poloidal average
        p1 = Process(\
                    target = combined1D2D                   ,\
                    args   = (path,)                          ,\
                    kwargs = {'xguards'        :xguards       ,\
                              'yguards'        :yguards       ,\
                              'marker'         :marker        ,\
                              'xSlice'         :xSlice        ,\
                              'ySlice'         :ySlice        ,\
                              'zSlice'         :zSlice        ,\
                              'tSlice'         :tSlice        ,\
                              'showPlot'       :showPlot      ,\
                              'savePlot'       :savePlot      ,\
                              'saveFolder'     :saveFolder    ,\
                              'saveFolderFunc' :saveFolderFunc,\
                              'varMax'         :varMax        ,\
                              'varMin'         :varMin        ,\
                              'varyMaxMin'     :varyMaxMin    ,\
                              'useSubProcess'  :useSubProcess ,\
                              'timeFolder'     :timeFolder    ,\
                              'polAvg'         :False         ,\
                              **kwargs
                             }
                    )

        # Plots with poloidal average
        p2 = Process(\
                    target = combined1D2D                   ,\
                    args   = (path,)                          ,\
                    kwargs = {'xguards'        :xguards       ,\
                              'yguards'        :yguards       ,\
                              'marker'         :marker        ,\
                              'xSlice'         :xSlice        ,\
                              'ySlice'         :ySlice        ,\
                              'zSlice'         :zSlice        ,\
                              'tSlice'         :tSlice        ,\
                              'showPlot'       :showPlot      ,\
                              'savePlot'       :savePlot      ,\
                              'saveFolder'     :saveFolder    ,\
                              'saveFolderFunc' :saveFolderFunc,\
                              'varMax'         :varMax        ,\
                              'varMin'         :varMin        ,\
                              'varyMaxMin'     :varyMaxMin    ,\
                              'useSubProcess'  :useSubProcess ,\
                              'timeFolder'     :timeFolder    ,\
                              'polAvg'         :True          ,\
                              **kwargs
                             }
                    )

        p1.start()
        p2.start()
        #}}}
    else:
        #{{{ Normal function call
        # Plots without poloidal average
        timeFolder = combined1D2D(path                            ,\
                                  xguards         = xguards       ,\
                                  yguards         = yguards       ,\
                                  marker          = marker        ,\
                                  xSlice          = xSlice        ,\
                                  ySlice          = ySlice        ,\
                                  zSlice          = zSlice        ,\
                                  tSlice          = tSlice        ,\
                                  showPlot        = showPlot      ,\
                                  savePlot        = savePlot      ,\
                                  saveFolder      = saveFolder    ,\
                                  saveFolderFunc  = saveFolderFunc,\
                                  varMax          = varMax        ,\
                                  varMin          = varMin        ,\
                                  varyMaxMin      = varyMaxMin    ,\
                                  useSubProcess   = useSubProcess ,\
                                  timeFolder      = timeFolder    ,\
                                  polAvg          = False         ,\
                                  **kwargs)

        # Plots with poloidal average
        timeFolder = combined1D2D(path                            ,\
                                  xguards         = xguards       ,\
                                  yguards         = yguards       ,\
                                  marker          = marker        ,\
                                  xSlice          = xSlice        ,\
                                  ySlice          = ySlice        ,\
                                  zSlice          = zSlice        ,\
                                  tSlice          = tSlice        ,\
                                  showPlot        = showPlot      ,\
                                  savePlot        = savePlot      ,\
                                  saveFolder      = saveFolder    ,\
                                  saveFolderFunc  = saveFolderFunc,\
                                  varMax          = varMax        ,\
                                  varMin          = varMin        ,\
                                  varyMaxMin      = varyMaxMin    ,\
                                  useSubProcess   = useSubProcess ,\
                                  timeFolder      = timeFolder    ,\
                                  polAvg          = True          ,\
                                  **kwargs)
        #}}}

    return timeFolder
#}}}

#{{{combined1D2D
def combined1D2D(path                          ,\
                 xguards        = False        ,\
                 yguards        = False        ,\
                 marker         = 'o'          ,\
                 xSlice         = slice(0,None),\
                 ySlice         = slice(0,None),\
                 zSlice         = slice(0,None),\
                 tSlice         = None         ,\
                 polAvg         = False         ,\
                 showPlot       = False        ,\
                 savePlot       = True         ,\
                 saveFolder     = None         ,\
                 saveFolderFunc = None         ,\
                 varMax         = None         ,\
                 varMin         = None         ,\
                 varyMaxMin     = None         ,\
                 useSubProcess  = True         ,\
                 timeFolder     = None         ,\
                 **kwargs):
    #{{{docstring
    """
    Driver for combined 1D and 2D plots.

    Input:
    path           - The simulation folder
    xguards        - If the xguards are to be collected
    yguards        - If the yguards are to be collected
    marker         - Marker in the plots
    xSlice         - How to slice in x
    ySlice         - How to slice in y
    zSlice         - How to slice in z
    tSlice         - How to slice in t
    polAvg         - Whether or not to perform a poloidal average of
                     the data
    showPlot       - If the plot is to be displayed
    savePlot       - If the plot is to be saved
    saveFolder     - Name of save folder
    saveFolderFunc - Function which takes this input as an input and
                     returns the saveFolder. To simplify the calling
                     sequence when running on a server, this is given
                     as a string.
    varMax         - Setting a hard upper limit z-axis in the plot
    varMin         - Setting a hard lower limit z-axis in the plot
    varyMaxMin     - Whether or not the limits of the z-axis should be set
                     to the max/min of the current timestep or not
    useSubProcess  - Each plot will be made by a new sub process
    timeFolder     - The name of the timefolder
                     Enables several plots to be put into same timeFolder
    **kwargs       - Extra keyword arguments used in input functions

    Output:
    timeFolder - The name of the timefolder
                 Enables several plots to be put into same timeFolder
    """
    #}}}

    if saveFolderFunc is not None:
        if saveFolderFunc == 'scanWTagSaveFunc':
            saveFolder = scanWTagSaveFunc(path                   ,\
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
    else:
        if saveFolder is None:
            saveFolder = "-".join(path.split('/')[::-1])

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
                                                   prePaths   = prePaths,\
                                                   )

        # Do the 1D plots
        p1 = Process(\
                    target = parPerpDriver                 ,\
                    args   = (path,)                       ,\
                    kwargs = {'xguards'      :xguards      ,\
                              'yguards'      :yguards      ,\
                              'marker'       :marker       ,\
                              'xSlice'       :xSlice       ,\
                              'ySlice'       :ySlice       ,\
                              'zSlice'       :zSlice       ,\
                              'tSlice'       :tSlice       ,\
                              'polAvg'       :polAvg       ,\
                              'showPlot'     :showPlot     ,\
                              'savePlot'     :savePlot     ,\
                              'saveFolder'   :saveFolder   ,\
                              'timeFolder'   :timeFolder   ,\
                              'useSubProcess':useSubProcess,\
                             }
                    )

        # Do the 2D plots
        p2 = Process(\
                    target = allMainFields2DDriver         ,\
                    args   = (path,)                       ,\
                    kwargs = {'xguards'      :xguards      ,\
                              'yguards'      :yguards      ,\
                              'xSlice'       :xSlice       ,\
                              'ySlice'       :ySlice       ,\
                              'zSlice'       :zSlice       ,\
                              'tSlice'       :tSlice       ,\
                              'polAvg'       :polAvg       ,\
                              'showPlot'     :showPlot     ,\
                              'savePlot'     :savePlot     ,\
                              'saveFolder'   :saveFolder   ,\
                              'varMax'       :varMax       ,\
                              'varMin'       :varMin       ,\
                              'varyMaxMin'   :varyMaxMin   ,\
                              'timeFolder'   :timeFolder   ,\
                              'useSubProcess':useSubProcess,\
                             }
                    )

        p1.start()
        p2.start()
        #}}}
    else:
        #{{{ Normal function call
        # Do the 1D plots
        timeFolder = parPerpDriver(path                         ,\
                                   xguards       = xguards      ,\
                                   yguards       = yguards      ,\
                                   marker        = marker       ,\
                                   xSlice        = xSlice       ,\
                                   ySlice        = ySlice       ,\
                                   zSlice        = zSlice       ,\
                                   tSlice        = tSlice       ,\
                                   polAvg        = polAvg       ,\
                                   showPlot      = showPlot     ,\
                                   savePlot      = savePlot     ,\
                                   saveFolder    = saveFolder   ,\
                                   timeFolder    = timeFolder   ,\
                                   useSubProcess = useSubProcess,\
                                   )
        # Do the 2D plots
        timeFolder = allMainFields2DDriver(path                         ,\
                                           xguards       = xguards      ,\
                                           yguards       = yguards      ,\
                                           xSlice        = xSlice       ,\
                                           ySlice        = ySlice       ,\
                                           zSlice        = zSlice       ,\
                                           tSlice        = tSlice       ,\
                                           polAvg        = polAvg       ,\
                                           showPlot      = showPlot     ,\
                                           savePlot      = savePlot     ,\
                                           saveFolder    = saveFolder   ,\
                                           varMax        = varMax       ,\
                                           varMin        = varMin       ,\
                                           varyMaxMin    = varyMaxMin   ,\
                                           timeFolder    = timeFolder   ,\
                                           useSubProcess = useSubProcess,\
                                          )
        #}}}

    return timeFolder
#}}}

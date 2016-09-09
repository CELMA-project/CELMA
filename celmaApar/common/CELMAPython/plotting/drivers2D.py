#!/usr/bin/env python

"""
Contains drivers for plotting 2D plots
"""

from .plotters import Plot2D
from .getStrings import getSaveString
import matplotlib.pyplot as plt
import numpy as np
from multiprocessing import Process
from boutdata import collect

#{{{allMainFields2DDriver
def allMainFields2DDriver(path                ,\
                          useSubProcess = True,\
                          saveFolder    = None,\
                          timeFolder    = None,\
                          **kwargs):
    #{{{docstring
    """
    Driver for all main fields.

    Specific allMainFields input:
    useSubProcess - Each plot will be made by a new sub process
    timeFolder    - The name of the timefolder
                    Enables several plots to be put into same timeFolder

    Specific allMainFields output:
    timeFolder - The name of the timefolder
                 Enables several plots to be put into same timeFolder

    For more details, see single2DDriver
    """
    #}}}

    functions = [\
                 lnN2DDriver        ,\
                 uEPar2DDriver      ,\
                 uIPar2DDriver      ,\
                 momDens2DDriver    ,\
                 JMPar2DDriver,\
                 vortD2DDriver      ,\
                 vort2DDriver       ,\
                 phi2DDriver        ,\
                 APar2DDriver       ,\
                 jPar2DDriver       ,\
                 n2DDriver          ,\
                ]

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
        for func in functions:
            processes.append(\
                Process(\
                        target = func                    ,\
                        args   = (path,)                 ,\
                        kwargs = {'saveFolder':saveFolder,\
                                  'timeFolder':timeFolder,\
                                  **kwargs
                                 }
                       )
            )

        for process in processes:
            process.start()
        #}}}
    else:
        #{{{ Normal function call
        # Do the plotting
        for func in functions:
            timeFolder = func(path                   ,\
                              timeFolder = timeFolder,\
                              saveFolder = saveFolder,\
                              **kwargs)
        #}}}

    return timeFolder
#}}}

#{{{lnN2DDriver
def lnN2DDriver(path             ,\
                timeFolder = None,\
                **kwargs):
    #{{{docstring
    """
    Driver for the lnN field.

    Specific lnN input:
    timeFolder - The name of the timefolder
                 Enables several plots to be put into same timeFolder

    Specific lnN output:
    timeFolder - The name of the timefolder
                 Enables several plots to be put into same timeFolder

    For more details, see single2DDriver
    """
    #}}}

    varName = 'lnN'
    pltName = r'\ln(n)'

    # Do the plot
    timeFolder = single2DDriver(path                   ,\
                                varName                ,\
                                pltName    = pltName   ,\
                                timeFolder = timeFolder,\
                                **kwargs)

    return timeFolder
#}}}

#{{{uEPar2DDriver
def uEPar2DDriver(path             ,\
                  timeFolder = None,\
                  **kwargs):
    #{{{docstring
    """
    Driver for the uEPar field.

    Specific uEPar input:
    timeFolder - The name of the timefolder
                 Enables several plots to be put into same timeFolder

    Specific uEPar output:
    timeFolder - The name of the timefolder
                 Enables several plots to be put into same timeFolder

    For more details, see single2DDriver
    """
    #}}}

    varName = 'uEPar'
    pltName = r'u_{e,\parallel}'

    # Do the plot
    timeFolder = single2DDriver(path                   ,\
                                varName                ,\
                                pltName    = pltName   ,\
                                timeFolder = timeFolder,\
                                **kwargs)

    return timeFolder
#}}}

#{{{uIPar2DDriver
def uIPar2DDriver(path             ,\
                   timeFolder = None,\
                   **kwargs):
    #{{{docstring
    """
    Driver for the uIPar field.

    Specific uIPar input:
    timeFolder - The name of the timefolder
                 Enables several plots to be put into same timeFolder

    Specific uIPar output:
    timeFolder - The name of the timefolder
                 Enables several plots to be put into same timeFolder

    For more details, see single2DDriver
    """
    #}}}

    varName = 'uIPar'
    pltName = r'u_{i,\parallel}'

    # Do the plot
    timeFolder = single2DDriver(path                   ,\
                                varName                ,\
                                pltName    = pltName   ,\
                                timeFolder = timeFolder,\
                                **kwargs)

    return timeFolder
#}}}

#{{{vortD2DDriver
def vortD2DDriver(path             ,\
                   timeFolder = None,\
                   **kwargs):
    #{{{docstring
    """
    Driver for the vortD field.

    Specific vortD input:
    timeFolder - The name of the timefolder
                 Enables several plots to be put into same timeFolder

    Specific vortD output:
    timeFolder - The name of the timefolder
                 Enables several plots to be put into same timeFolder

    For more details, see single2DDriver
    """
    #}}}

    varName = 'vortD'
    pltName = r'\Omega^D'

    # Do the plot
    timeFolder = single2DDriver(path                   ,\
                                varName                ,\
                                pltName    = pltName   ,\
                                timeFolder = timeFolder,\
                                **kwargs)

    return timeFolder
#}}}

#{{{vort2DDriver
def vort2DDriver(path             ,\
                 timeFolder = None,\
                 **kwargs):
    #{{{docstring
    """
    Driver for the vort field.

    Specific vort input:
    timeFolder - The name of the timefolder
                 Enables several plots to be put into same timeFolder

    Specific vort output:
    timeFolder - The name of the timefolder
                 Enables several plots to be put into same timeFolder

    For more details, see single2DDriver
    """
    #}}}

    varName = 'vort'
    pltName = r'\Omega'

    # Do the plot
    timeFolder = single2DDriver(path                   ,\
                                varName                ,\
                                pltName    = pltName   ,\
                                timeFolder = timeFolder,\
                                **kwargs)

    return timeFolder
#}}}

#{{{JMPar2DDriver
def JMPar2DDriver(path             ,\
                        timeFolder = None,\
                        **kwargs):
    #{{{docstring
    """
    Driver for the jMPar field.

    Specific jMPar input:
    timeFolder - The name of the timefolder
                 Enables several plots to be put into same timeFolder

    Specific jMPar output:
    timeFolder - The name of the timefolder
                 Enables several plots to be put into same timeFolder

    For more details, see single2DDriver
    """
    #}}}

    varName = 'jMPar'
    pltName = r'j^M_\parallel'

    # Do the plot
    timeFolder = single2DDriver(path                   ,\
                                varName                ,\
                                pltName    = pltName   ,\
                                timeFolder = timeFolder,\
                                **kwargs)

    return timeFolder

#}}}

#{{{momDens2DDriver
def momDens2DDriver(path             ,\
                    timeFolder = None,\
                    **kwargs):
    #{{{docstring
    """
    Driver for the jPar field.

    Specific jPar input:
    timeFolder - The name of the timefolder
                 Enables several plots to be put into same timeFolder

    Specific jPar output:
    timeFolder - The name of the timefolder
                 Enables several plots to be put into same timeFolder

    For more details, see single2DDriver
    """
    #}}}

    varName = 'momDensPar'
    pltName = r'nu_{i,\parallel}'

    # Do the plot
    timeFolder = single2DDriver(path                   ,\
                                varName                ,\
                                pltName    = pltName   ,\
                                timeFolder = timeFolder,\
                                **kwargs)

    return timeFolder

#}}}

#{{{phi2DDriver
def phi2DDriver(path             ,\
                timeFolder = None,\
                **kwargs):
    #{{{docstring
    """
    Driver for the phi field.

    Specific phi input:
    timeFolder - The name of the timefolder
                 Enables several plots to be put into same timeFolder

    Specific phi output:
    timeFolder - The name of the timefolder
                 Enables several plots to be put into same timeFolder

    For more details, see single2DDriver
    """
    #}}}

    varName = 'phi'
    pltName = r'\phi'

    # Do the plot
    timeFolder = single2DDriver(path                   ,\
                                varName                ,\
                                pltName    = pltName   ,\
                                timeFolder = timeFolder,\
                                **kwargs)

    return timeFolder
#}}}

#{{{APar2DDriver
def APar2DDriver(path             ,\
                timeFolder = None,\
                **kwargs):
    #{{{docstring
    """
    Driver for the APar field.

    Specific APar input:
    timeFolder - The name of the timefolder
                 Enables several plots to be put into same timeFolder

    Specific APar output:
    timeFolder - The name of the timefolder
                 Enables several plots to be put into same timeFolder

    For more details, see single2DDriver
    """
    #}}}

    varName = 'APar'
    pltName = r'A_\parallel'

    # Do the plot
    timeFolder = single2DDriver(path                   ,\
                                varName                ,\
                                pltName    = pltName   ,\
                                timeFolder = timeFolder,\
                                **kwargs)

    return timeFolder

#}}}

#{{{jPar2DDriver
def jPar2DDriver(path             ,\
                timeFolder = None,\
                **kwargs):
    #{{{docstring
    """
    Driver for the jPar field.

    Specific jPar input:
    timeFolder - The name of the timefolder
                 Enables several plots to be put into same timeFolder

    Specific jPar output:
    timeFolder - The name of the timefolder
                 Enables several plots to be put into same timeFolder

    For more details, see single2DDriver
    """
    #}}}

    varName = 'jPar'
    pltName = r'j_\parallel'

    # Do the plot
    timeFolder = single2DDriver(path                   ,\
                                varName                ,\
                                pltName    = pltName   ,\
                                timeFolder = timeFolder,\
                                **kwargs)

    return timeFolder

#}}}

#{{{n2DDriver
def n2DDriver(path             ,\
              timeFolder = None,\
              **kwargs):
    #{{{docstring
    """
    Driver for the n field.

    Specific n input:
    timeFolder - The name of the timefolder
                 Enables several plots to be put into same timeFolder

    Specific n output:
    timeFolder - The name of the timefolder
                 Enables several plots to be put into same timeFolder

    For more details, see single2DDriver
    """
    #}}}

    varName = None
    pltName = r'n'

    #{{{varFunc
    def varFunc(path    = None,\
                yguards = None,\
                xguards = None,\
                tind    = None,\
                **kwargs):
        """
        Function which returns the parallel current variable.

        This function will be called in plotters.Plot2D
        """

        lnN = collect("lnN"               ,\
                      path    = path      ,\
                      yguards = yguards   ,\
                      xguards = xguards   ,\
                      tind    = tind      ,\
                      info    = False     ,\
                      )

        n = np.exp(lnN)

        return n
    #}}}

    # Do the plot
    timeFolder = single2DDriver(path                   ,\
                                varName                ,\
                                pltName    = pltName   ,\
                                timeFolder = timeFolder,\
                                varFunc    = varFunc   ,\
                                **kwargs)

    return timeFolder
#}}}

#{{{single2DDriver
def single2DDriver(path                      ,\
                   varName                   ,\
                   xguards    = False        ,\
                   yguards    = False        ,\
                   xSlice     = slice(0,None),\
                   ySlice     = slice(0,None),\
                   zSlice     = slice(0,None),\
                   tSlice     = None         ,\
                   subPolAvg  = False        ,\
                   showPlot   = False        ,\
                   savePlot   = True         ,\
                   saveFolder = None         ,\
                   pltName    = None         ,\
                   varMax     = None         ,\
                   varMin     = None         ,\
                   varyMaxMin = None         ,\
                   timeFolder = None         ,\
                   varFunc    = None         ,\
                   ):
    #{{{docstring
    """
    Driver for a single plot.

    Input:
    path       - The simulation folder
    varName    - Name of the field which is going to be plotted
    xguards    - If the xguards are to be collected
    yguards    - If the yguards are to be collected
    xSlice     - How to slice in x
    ySlice     - How to slice in y
    zSlice     - How to slice in z
    tSlice     - How to slice in t
    subPolAvg  - Whether or not the poloidal average should be
                 subtracted from the data
    showPlot   - If the plot is to be displayed
    savePlot   - If the plot is to be saved
    saveFolder - Name of save folder
    pltName    - Name of plot to make
    varMax     - Setting a hard upper limit z-axis in the plot
    varMin     - Setting a hard lower limit z-axis in the plot
    varyMaxMin - Whether or not the limits of the z-axis should be set
                 to the max/min of the current timestep or not
    timeFolder - The name of the timefolder
                 Enables several plots to be put into same timeFolder
    varFunc    - Function which returns the variable (used if variables
                 is not collectable)

    Output:
    timeFolder - The name of the timefolder
                 Enables several plots to be put into same timeFolder
    """
    #}}}

    # Make the plotter object
    plotter = Plot2D(path                   ,\
                     varName                ,\
                     xguards    = xguards   ,\
                     yguards    = yguards   ,\
                     xSlice     = xSlice    ,\
                     ySlice     = ySlice    ,\
                     zSlice     = zSlice    ,\
                     tSlice     = tSlice    ,\
                     subPolAvg  = subPolAvg ,\
                     showPlot   = showPlot  ,\
                     savePlot   = savePlot  ,\
                     saveFolder = saveFolder,\
                     varMax     = varMax    ,\
                     varMin     = varMin    ,\
                     varyMaxMin = varyMaxMin,\
                     varFunc    = varFunc   ,\
                    )

    # Do the plot
    timeFolder = plotter.plotDriver(pltName, timeFolder=timeFolder)

    return timeFolder
#}}}

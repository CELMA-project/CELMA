#!/usr/bin/env python

"""
Contains drivers for plotting 1D plots
"""

from .organizer import Organizer
from .line import Line
from .plotters import Plot1D
from .getStrings import getSaveString
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
            additional = ['visualization', saveFolder]
            saveString, timeFolder = getSaveString(''                     ,\
                                                   path                   ,\
                                                   timeFolder = timeFolder,\
                                                   additional = additional,\
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
                  "uEParFields",\
                  "uIParFields",\
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
            additional = ['visualization', saveFolder]
            saveString, timeFolder = getSaveString(''                     ,\
                                                   path                   ,\
                                                   timeFolder = timeFolder,\
                                                   additional = additional,\
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
                  "uEParFields",\
                  "uIParFields",\
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
            additional = ['visualization', saveFolder]
            saveString, timeFolder = getSaveString(''                     ,\
                                                   path                   ,\
                                                   timeFolder = timeFolder,\
                                                   additional = additional,\
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
                   polAvg     = False        ,\
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
    polAvg     - Whether or not to perform a poloidal average of
                 the data
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
                     polAvg     = polAvg    ,\
                     showPlot   = showPlot  ,\
                     savePlot   = savePlot  ,\
                     saveFolder = saveFolder,\
                    )

    # Select plt type
    if pltName == None:
        pltName = saveFolder

    if pltName == "mainFields":
        orgObj = getMainFields()
    elif pltName == "lnNFields":
        orgObj = getLnNFields()
    elif pltName == "uEParFields":
        orgObj = getUEParFields()
    elif pltName == "uIParFields":
        orgObj = getUIParFields()
    elif pltName == "vortDFields":
        orgObj = getVortDFields()
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
        # Create the jPar line (the first extra line)
        orgObj.extraLines['jPar'].field = np.exp(lnN)*(uIPar - uEPar)
        if orgObj.extraLines['jPar'].plotPos:
            orgObj.lines.insert(orgObj.extraLines['jPar'].plotPos,\
                                orgObj.extraLines['jPar'])
        else:
            orgObj.lines.append(orgObj.extraLines['jPar'])
        # Create the n line (the second extra line)
        orgObj.extraLines['n'].field = np.exp(lnN)
        if orgObj.extraLines['n'].plotPos:
            orgObj.lines.insert(orgObj.extraLines['n'].plotPos,\
                                orgObj.extraLines['n'])
        else:
            orgObj.lines.append(orgObj.extraLines['n'])

    # Do the plot
    timeFolder = plotter.plotDriver(fig, orgObj, timeFolder=timeFolder)

    return timeFolder
#}}}
#}}}

#{{{Field getters
#{{{getMainFields
def getMainFields():
    """
    Prepares the main fields
    """

    # Making the orgObj instance
    mainFields = Organizer("mainFields")
    # Making lines in the pattern name, lable, plotPos
    # Evolved fields
    mainFields.lines.append(Line('lnN'  , r'\ln(n)'         , plotPos=0))
    mainFields.lines.append(Line('uIPar', r'u_{i,\parallel}', plotPos=6))
    mainFields.lines.append(Line('uEPar', r'u_{e,\parallel}', plotPos=4))
    mainFields.lines.append(Line('vortD', r'\Omega^D'       , plotPos=3))
    # Helping field
    mainFields.lines.append(Line('phi'  , r'\phi'           , plotPos=1))
    mainFields.lines.append(Line('S'    , r'S'              , plotPos=7))
    mainFields.lines.append(Line('vort' , r'\Omega'         , plotPos=5))
    # Extra lines
    # FIXME: This API is not so intuitive, consider change
    mainFields.extraLines['jPar'] = Line('jPar', r'j_\parallel', plotPos=8)
    mainFields.extraLines['n'   ] = Line('n'   , r'n'          , plotPos=2)

    return mainFields
#}}}

#{{{getLnNFields
def getLnNFields():
    """
    Prepares the lnN fields
    """

    # Making the orgObj instance
    lnN = Organizer(r"\ln(n)", useCombinedPlot=True)
    # Making lines in the pattern name, lable, plotPos
    # Evolved fields
    lnN.lines.append(Line('lnNAdv'                             ,\
        r'-J\{\phi,\ln(n)\}'                                   ))
    lnN.lines.append(Line('lnNRes'                             ,\
        r'\frac{0.51\nu_{ei}}{\mu}\left(\nabla^2_\perp\ln(n) +'+\
        r'\left[\nabla_\perp\ln(n)\right]^2\right)'            ))
    lnN.lines.append(Line('gradUEPar'                          ,\
        r'-\partial_{\parallel}u_{e,\parallel}'                ))
    lnN.lines.append(Line('lnNUeAdv'                           ,\
        r'-u_{e,\parallel}\cdot\nabla_\parallel\ln(n)'         ))
    lnN.lines.append(Line('srcN'                               ,\
        r'\frac{S}{n}'                                         ))
    lnN.lines.append(Line('lnNPerpArtVisc'                     ,\
        r'-D_{n,\perp} \partial^2_{\perp}\ln(n)'               ))
    lnN.lines.append(Line('lnNParArtVisc'                      ,\
        r'-D_{n,\parallel} \partial^2_{\parallel}\ln(n)'       ))

    return lnN
#}}}

#{{{getUEParFields
def getUEParFields():
    """
    Prepares the uEPar fields
    """

    # Making the orgObj instance
    uEPar = Organizer(r"u_{e,\parallel}", useCombinedPlot=True)
    # Making lines in the pattern name, lable, plotPos
    uEPar.lines.append(Line('uEParAdv'                                      ,\
    r'-r\{\phi,u_{e,\parallel}\}'                                           ))
    uEPar.lines.append(Line('uEParParAdv'                                   ,\
    r'- u_{e,\parallel}\partial_{\parallel}u_{e,\parallel}'                 ))
    uEPar.lines.append(Line('gradPhiLnN'                                    ,\
    r'\mu\nabla_\parallel\left( \phi - \ln(n)\right)'                       ))
    uEPar.lines.append(Line('uEParRes'                                      ,\
    r'-0.51\nu_{ei}\left(u_{e,\parallel}-u_{i,\parallel} \right)'           ))
    uEPar.lines.append(Line('ueSrc'                                         ,\
    r'-\frac{S u_{e,\parallel}}{n}'                                         ))
    uEPar.lines.append(Line('ueNeutral'                                     ,\
    r'-\nu_{en}u_{e,\parallel}'                                             ))
    uEPar.lines.append(Line('uEParPerpArtVisc'                              ,\
    r'-D_{u_{e,\parallel}, \perp}\nabla^2_\perp u_{e,\parallel}'            ))
    uEPar.lines.append(Line('uEParParArtVisc'                               ,\
    r'-D_{u_{e,\parallel}, \parallel} \partial^2_{\parallel}u_{e,\parallel}'))

    return uEPar
#}}}

#{{{getUIParFields
def getUIParFields():
    """
    Prepares the uIPar fields
    """

    # Making the orgObj instance
    uIPar = Organizer(r"u_{i,\parallel}", useCombinedPlot=True)
    # Making lines in the pattern name, lable, plotPos
    uIPar.lines.append(Line('uIParAdv'                                    ,\
             r'-J\{\phi,u_{i,\parallel}\}'                                ))
    uIPar.lines.append(Line('uIParParAdv'                                 ,\
             r'-u_{i,\parallel}\partial_{\parallel}u_{i,\parallel}'       ))
    uIPar.lines.append(Line('gradPhi'                                     ,\
             r'-\nabla\left(\phi\right)'                                  ))
    uIPar.lines.append(Line('uIParRes'                                    ,\
             r'-0.51\nu_{ei}\left(u_{i,\parallel}-u_{e,\parallel} \right)'))
    uIPar.lines.append(Line('uIParPerpArtVisc'                            ,\
             r'D_{u_{i,\parallel}, \perp}\nabla^2_\perp u_{i,\parallel}'  ))
    uIPar.lines.append(Line('uiNeutral'                                   ,\
             r'-\nu_{ei}u_{i,\parallel}'                                  ))
    uIPar.lines.append(Line('uiSrc'                                       ,\
             r'-\frac{S u_{i,\parallel}}{n}'                              ))
    uIPar.lines.append(Line('uIParParArtVisc'                             ,\
             r'-D_{u_i,\parallel} \partial^2_{\parallel}u_{i,\parallel}'  ))

    return uIPar
#}}}

#{{{getVortDFields
def getVortDFields():
    """
    Prepares the vortD fields
    """

    # Making the orgObj instance
    vortD = Organizer(r"\Omega^D", useCombinedPlot=True)
    # Making lines in the pattern name, lable, plotPos
    vortD.lines.append(Line('vortNeutral'                        ,\
         r'-\nu_{in}n\Omega'                                     ))
    vortD.lines.append(Line('potNeutral'                         ,\
         r'-\nu_{in}\nabla_\perp \phi \cdot \nabla_\perp n'      ))
    vortD.lines.append(Line('divExBAdvGradPerpPhiN'              ,\
         r'-\nabla\cdot\left('                                   +\
         r'\mathbf{u}_E \cdot\nabla'                             +\
         r'\left[ n \nabla_\perp \phi \right]'                   +\
         r'\right)'                                              ))
    vortD.lines.append(Line('parDerDivUIParNGradPerpPhi'         ,\
         r'-\partial_\parallel \nabla\cdot\left('                +\
         r'u_{i,\parallel} '                                     +\
         r'n \nabla_\perp \phi'                                  +\
         r'\right)'                                              ))
    vortD.lines.append(Line('nGradUiUe'                          ,\
         r'n\partial_\parallel (u_{i,\parallel}-u_{e,\parallel})'))
    vortD.lines.append(Line('uiUeGradN'                          ,\
        r'(u_{i,\parallel}-u_{e,\parallel})\partial_{\parallel}n'))
    vortD.lines.append(Line('vortDParArtVisc'                    ,\
         r'-D_{\Omega^D} \partial^2_{\parallel}\Omega^D'         ))
    vortD.lines.append(Line('vortDPerpArtVisc'                   ,\
         r'-D_{\Omega^D, \perp} \nabla_\perp^2\Omega^D'          ))

    return vortD
#}}}
#}}}

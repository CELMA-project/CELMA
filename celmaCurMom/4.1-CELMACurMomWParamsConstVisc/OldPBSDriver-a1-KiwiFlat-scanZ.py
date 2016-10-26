#!/usr/bin/env python

"""Driver which runs using PBS."""

import os, sys
# If we add to sys.path, then it must be an absolute path
commonDir = os.path.abspath('./../common')
# Sys path is a list of system paths
sys.path.append(commonDir)

import re
import numpy as np
from boutdata import collect
from bout_runners import PBS_runner
from CELMAPython.drivers import postBoutRunner

#{{{restartFromFunc
def restartFromFunc(dmp_folder     = None,\
                    aScanPath      = None,\
                    scanParameters = None,\
                    **kwargs):
    """
    Function which converts a path belonging to one paths in a scan
    to the path belonging to the current scan.

    The function obtains the current scan parameters from dmp_folder
    (the dmp_folder given from bout_runners), and inserts the
    current scan parameters into aScanPath (the function input which is
    one of the paths belonging to the scan).

    NOTE: This will not work if the values of one of the scan parameters
          contains an underscore.

    Parameters
    ----------
    dmp_folder : str
        Given by the bout_runners. Used to find the current scan
        values.
    aScanPath : str
        One of the restart paths from a previously run simulation.
    scanParameters : list
        List of strings of the names of the scan paramemters.
    **kwargs : keyword dictionary
        Dictionary with additional keyword arguments, given by
        bout_runners.

    Returns
    -------
    scanPath : str
        aScanPath converted to the scan parameters of the current run.
    """

    # Make a template string of aScanPath
    scanPathTemplate = aScanPath
    for scanParameter in scanParameters:
        hits = [m.start() for m in \
                re.finditer(scanParameter, scanPathTemplate)]
        # FIXME: If initial hit is in root folder this algorithm will fail
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

    # Get the values from the current dmp_folder
    values = {}
    for scanParameter in scanParameters:
        hits = [m.start() for m in \
                re.finditer(scanParameter, dmp_folder)]
        # Choose the first hit to get the value from (again we assume
        # that the value does not contain a _)
        value_start = hits[0] + len(scanParameter) + 1
        # Here we assume that the value is not separated by an
        # underscore
        values[scanParameter] = dmp_folder[value_start:].split("_")[0]

    # Insert the values
    restartFrom = scanPathTemplate.format(values)

    return restartFrom
#}}}

# If you just want to post-process
justPostProcess = True
# Normal post-processors
postProcessInit = True
postProcessExp  = True
postProcessLin  = True
postProcessTurb = True
# Extra post-processors
postProcessLinProfiles     = False
postProcessTurbProfiles    = False
postProcessProbesAndEnergy = True
postProcessGrowthRates     = True

#{{{Main options
#{{{The scan
# NOTE: Calling this len will overshadow the len() function
length = [1       , 2       , 4       , 6       , 8       , 10       ]
Ly     = [102.2235, 204.4469, 408.8938, 613.3408, 817.7877, 1022.2346]
scanParameters  = ["len", "Ly"]
series_add = [\
              ('input', 'len', length),\
              ('geom' , 'Ly' , Ly),\
             ]
#}}}
#{{{The options for the post processing function
saveFolderFunc         = "scanWTagSaveFunc"
convertToPhysical      = False
showPlot               = False
savePlot               = True
xguards                = False
yguards                = False
xSlice                 = 0
ySlice                 = 8*2
zSlice                 = 0
axisEqualParallel      = False
varName                = "n"
pltName                = "n"
nProbes                = 5
maxMode                = 10
yInd                   = ySlice
var                    = "n"
useSteadyStatePathFunc = True
extension              = "png"
useSubProcess          = True
#}}}
#{{{File handeling options
remove_old = False
directory  = "a1-KiwiFlatZ"
make       = False
cpy_source = True
#}}}
#{{{The PBS options
nproc                 = 48
BOUT_nodes            = 3
BOUT_ppn              = 16
post_process_nproc    = 1
post_process_nodes    = 1
post_process_ppn      = 20
#}}}
#}}}

#{{{Abbrevations
#{{{Dictionaries with common runner options
commonRunnerKwargs =\
        {\
         "directory"          : directory         ,\
         "make"               : make              ,\
         "nproc"              : nproc             ,\
         "cpy_source"         : cpy_source        ,\
         "BOUT_nodes"         : BOUT_nodes        ,\
         "BOUT_ppn"           : BOUT_ppn          ,\
         "post_process_nproc" : post_process_nproc,\
         "post_process_nodes" : post_process_nodes,\
         "post_process_ppn"   : post_process_ppn  ,\
        }
#}}}
#{{{Dictionaries with common post-processing options
commonPlotterKwargs =\
        {\
         "saveFolderFunc"   : saveFolderFunc   ,\
         "convertToPhysical": convertToPhysical,\
         "showPlot"         : showPlot         ,\
         "savePlot"         : savePlot         ,\
         "extension"        : extension        ,\
         "useSubProcess"    : useSubProcess    ,\
        }
fieldPlotterKwargs =\
        {\
         "xguards"          : xguards          ,\
         "yguards"          : yguards          ,\
         "xSlice"           : xSlice           ,\
         "ySlice"           : ySlice           ,\
         "zSlice"           : zSlice           ,\
         "axisEqualParallel": axisEqualParallel,\
         **commonPlotterKwargs                 ,\
        }
#}}}
#}}}

#{{{Init runner
if postProcessInit:
    curPostProcessor = postBoutRunner
else:
    curPostProcessor = None
#{{{Init options
# Name
theRunName = "a1-KiwiFlatZ-0-initialize"
# Set the spatial domain
nz = 1
# Set the temporal domain
restart    = None
timestep   = [2e3]
nout       = [2]
# Filter
ownFilterType = "none"
#Switches
useHyperViscAzVortD = [False]
# Specify the numbers used for the BOUT runs
BOUT_walltime         = '08:00:00'
BOUT_run_name         = theRunName
post_process_run_name = 'post' + theRunName.capitalize()
post_process_walltime = '0:29:00'
post_process_queue    = 'xpresq'
# Post processing option
tSlice     = slice(-2, None)
#}}}
#{{{Run and post processing
initRunner = PBS_runner(\
                # Set spatial domain
                nz         = nz      ,\
                # Set temporal domain
                nout       = nout    ,\
                timestep   = timestep,\
                # Set the restart option
                restart    = restart ,\
                # Set additional option
                additional = [
                              ('tag',theRunName,0),\
                              ('ownFilters'  , 'type', ownFilterType),\
                              ('switch'      , 'useHyperViscAzVortD', useHyperViscAzVortD),\
                             ],\
                series_add = series_add                      ,\
                # PBS options
                BOUT_walltime         = BOUT_walltime        ,\
                BOUT_run_name         = BOUT_run_name        ,\
                post_process_walltime = post_process_walltime,\
                post_process_queue    = post_process_queue   ,\
                post_process_run_name = post_process_run_name,\
                # Common options
                **commonRunnerKwargs                         ,\
                )

init_dmp_folders, PBS_ids = initRunner.execute_runs(\
                             remove_old               = remove_old,\
                             post_processing_function = curPostProcessor,\
                             # This function will be called every time after
                             # performing a run
                             post_process_after_every_run = True,\
                             # Below are the kwargs arguments being passed to
                             # the post processing function
                             driverName        = "plot1DAnd2DDriver",\
                             tSlice            = tSlice             ,\
                             theRunName        = theRunName         ,\
                             # Common kwargs
                             **fieldPlotterKwargs                   ,\
                            )
#}}}
#}}}

if justPostProcess:
    restart = None
else:
    restart = "overwrite"

#{{{The expand runner
if postProcessExp:
    curPostProcessor = postBoutRunner
else:
    curPostProcessor = None
#{{{Expand runner options
# Set the spatial domain
nz = 256
# Set the temporal domain
timestep   = [50]
nout       = [2]
# Filter
ownFilterType = "none"
#Switches
useHyperViscAzVortD = [False]
# From previous outputs
aScanPath = init_dmp_folders[0]
# Name
theRunName = "a1-KiwiFlatZ-1-expand"
# PBS options
BOUT_walltime         = '08:00:00'
BOUT_run_name         = theRunName
post_process_run_name = 'post' + theRunName.capitalize()
post_process_walltime = '0:29:00'
post_process_queue    = 'xpresq'
# Post processing option
tSlice     = slice(-2, None)
#}}}
#{{{Run and post processing
expandRunner = PBS_runner(\
                # Set spatial domain
                nz           = nz               ,\
                # Set temporal domain
                nout         = nout             ,\
                timestep     = timestep         ,\
                # Set restart options
                restart      = restart          ,\
                restart_from = restartFromFunc  ,\
                # Set additional options
                additional = [
                              ('tag',theRunName,0),\
                              ('ownFilters'  , 'type', ownFilterType),\
                              ('switch'      , 'useHyperViscAzVortD', useHyperViscAzVortD),\
                             ],\
                series_add = series_add                      ,\
                # PBS options
                BOUT_walltime         = BOUT_walltime        ,\
                BOUT_run_name         = BOUT_run_name        ,\
                post_process_walltime = post_process_walltime,\
                post_process_queue    = post_process_queue   ,\
                post_process_run_name = post_process_run_name,\
                # Common options
                **commonRunnerKwargs                         ,\
                )

expand_dmp_folders, PBS_ids = expandRunner.execute_runs(\
                               remove_old               = remove_old  ,\
                               post_processing_function = curPostProcessor,\
                               # Declare dependencies
                               job_dependencies = PBS_ids,\
                               # This function will be called every time after
                               # performing a run
                               post_process_after_every_run = True,\
                               # Below are the kwargs arguments being passed to
                               # the post processing function
                               driverName        = "plot1DAnd2DDriver",\
                               tSlice            = tSlice             ,\
                               theRunName        = theRunName         ,\
                               # Below are the kwargs given to the
                               # restartFromFunc
                               aScanPath      = aScanPath     ,\
                               scanParameters = scanParameters,\
                               # Common kwargs
                               **fieldPlotterKwargs           ,\
                              )
#}}}
#}}}

#{{{ If profiles are to be plotted
if postProcessLinProfiles or postProcessTurbProfiles:
    noutProfile                 = 3
    timestepProfile             = 10
    restartProfile              = "overwrite"
    useHyperViscAzVortDProfile  = True
    saveTermsProfile            = True

    # Create the options for the runners
    # Notice that we would like to save all the fields here
    profileRunOptions = {\
                  # Set temporal domain
                  "nout"         : noutProfile    ,\
                  "timestep"     : timestepProfile,\
                  # Set restart options
                  "restart"      : restartProfile ,\
                  "restart_from" : restartFromFunc,\
                  # Set additional options
                  "additional" : [
                                ('switch'      , 'useHyperViscAzVortD',useHyperViscAzVortDProfile),\
                                ('switch'      , 'saveTerms'          ,saveTermsProfile),\
                               ],\
                  "series_add" : series_add       ,\
                  # Common options
                  **commonRunnerKwargs            ,\
                        }
#}}}

#{{{The linear runner
if postProcessLin:
    curPostProcessor = postBoutRunner
else:
    curPostProcessor = None
#{{{ Linear options
#Switches
saveTerms           = False
useHyperViscAzVortD = [True]
includeNoise     = True
forceAddNoise    = True
# As this is scan dependent, the driver finds the correct folder
maxGradRhoFolder = expand_dmp_folders[0]
# From previous outputs
aScanPath = expand_dmp_folders[0]
# Set the temporal domain
timestep  = [1]
nout     = [500]
# Name
theRunName = "a1-KiwiFlatZ-2-linearPhase1"
# PBS options
BOUT_run_name         = theRunName
BOUT_walltime         = '100:00:00'
post_process_run_name = 'post' + theRunName.capitalize()
post_process_walltime = '03:00:00'
post_process_queue    = 'workq'
# Post processing options
tSlice     = slice(0, None, 2)
varyMaxMin = True
subPolAvg  = True
mode       = "perpAndPol"
#}}}
#{{{Run and post processing
linearRun = PBS_runner(\
            # Set temporal domain
            nout         = nout           ,\
            timestep     = timestep       ,\
            # Set restart options
            restart      = restart        ,\
            restart_from = restartFromFunc,\
            # Set additional options
            additional = [
                          ('tag',theRunName,0),\
                          ('switch'      , 'includeNoise'       , includeNoise ),\
                          ('switch'      , 'forceAddNoise'      ,forceAddNoise),\
                          ('switch'      , 'useHyperViscAzVortD',useHyperViscAzVortD),\
                          ('switch'      , 'saveTerms'          ,saveTerms),\
                         ],\
            series_add = series_add                      ,\
            # PBS options
            BOUT_walltime         = BOUT_walltime        ,\
            BOUT_run_name         = BOUT_run_name        ,\
            post_process_walltime = post_process_walltime,\
            post_process_queue    = post_process_queue   ,\
            post_process_run_name = post_process_run_name,\
            # Common options
            **commonRunnerKwargs                         ,\
            )

linear_dmp_folders, PBS_ids = linearRun.execute_runs(\
                                 remove_old               = remove_old,\
                                 post_processing_function = curPostProcessor,\
                                 # Declare dependencies
                                 job_dependencies = PBS_ids,\
                                 # This function will be called every time after
                                 # performing a run
                                 post_process_after_every_run = True,\
                                 # Below are the kwargs arguments being passed to
                                 # the post processing function
                                 # Switches
                                 driverName       = "single2DDriver",\
                                 theRunName       = theRunName      ,\
                                 tSlice           = tSlice          ,\
                                 subPolAvg        = subPolAvg       ,\
                                 varName          = varName         ,\
                                 pltName          = pltName         ,\
                                 varyMaxMin       = varyMaxMin      ,\
                                 mode             = mode            ,\
                                 maxGradRhoFolder = maxGradRhoFolder,\
                                 # Below are the kwargs given to the
                                 # restartFromFunc
                                 aScanPath      = aScanPath     ,\
                                 scanParameters = scanParameters,\
                                 # Common kwargs
                                 **fieldPlotterKwargs           ,\
                                )
#}}}
#{{{ If linear profiles are to be plotted
if postProcessLinProfiles:
    curPostProcessor = postBoutRunner
    theRunName = "a1-KiwiFlatZ-2-linearPhaseParProfiles"
    aScanPathProfiles = linear_dmp_folders[0]
    tSlice = slice(-30, None, 10)

    _, _ = linearRun.execute_runs(\
                                 remove_old               = remove_old,\
                                 post_processing_function = curPostProcessor,\
                                 # Declare dependencies
                                 job_dependencies = PBS_ids,\
                                 # This function will be called every time after
                                 # performing a run
                                 post_process_after_every_run = True,\
                                 # Below are the kwargs arguments being passed to
                                 # the post processing function
                                 # Switches
                                 driverName     = "parDriver"   ,\
                                 tSlice         = tSlice        ,\
                                 theRunName     = theRunName    ,\
                                 # Below are the kwargs given to the
                                 # restartFromFunc
                                 aScanPath      = aScanPath     ,\
                                 scanParameters = scanParameters,\
                                 # Common kwargs
                                 **fieldPlotterKwargs           ,\
                                )
#}}}
#}}}

#{{{Turbulence runner
if postProcessTurb:
    curPostProcessor = postBoutRunner
else:
    curPostProcessor = None
#{{{Turbulence options
# Switches
saveTerms           = False
useHyperViscAzVortD = [True]
# Set the temporal domain
nout     = [5000]
timestep = [1]
# Name
theRunName = "a1-KiwiFlatZ-3-turbulentPhase1"
# PBS options
BOUT_run_name         = theRunName
BOUT_walltime         = '100:00:00'
post_process_run_name = 'post' + theRunName.capitalize()
post_process_walltime = '03:00:00'
post_process_queue    = 'workq'
# Post processing options
tSlice    = slice(-5000, None, 10)
aScanPath = linear_dmp_folders[0]
#}}}
#{{{Run and post processing
turboRun = PBS_runner(\
                # Set temporal domain
                nout       = nout               ,\
                timestep   = timestep           ,\
                # Set restart options
                restart      = restart          ,\
                restart_from = restartFromFunc,\
                additional = [
                              ('tag',theRunName,0),\
                              ('switch'      , 'useHyperViscAzVortD',useHyperViscAzVortD),\
                              ('switch'      , 'saveTerms'          ,saveTerms),\
                             ],\
                series_add = series_add                      ,\
                # PBS options
                BOUT_walltime         = BOUT_walltime        ,\
                BOUT_run_name         = BOUT_run_name        ,\
                post_process_walltime = post_process_walltime,\
                post_process_queue    = post_process_queue   ,\
                post_process_run_name = post_process_run_name,\
                # Common options
                **commonRunnerKwargs                         ,\
                )

turbo_dmp_folders, PBS_ids = turboRun.execute_runs(\
                                 remove_old               = remove_old,\
                                 post_processing_function = curPostProcessor,\
                                 # Declare dependencies
                                 job_dependencies = PBS_ids,\
                                 # This function will be called every time after
                                 # performing a run
                                 post_process_after_every_run = True,\
                                 # Below are the kwargs arguments being passed to
                                 # the post processing function
                                 # Switches
                                 driverName     = "single2DDriver",\
                                 tSlice         = tSlice          ,\
                                 theRunName     = theRunName      ,\
                                 varName        = varName         ,\
                                 pltName        = pltName         ,\
                                 # Below are the kwargs given to the
                                 # restartFromFunc
                                 aScanPath      = aScanPath     ,\
                                 scanParameters = scanParameters,\
                                 # Common kwargs
                                 **fieldPlotterKwargs           ,\
                                )
#}}}
#{{{ If linear profiles are to be plotted
if postProcessTurbProfiles:
    curPostProcessor = postBoutRunner
    theRunName = "a1-KiwiFlatElTemp-3-turbulentPhase1ParProfiles"
    aScanPathProfiles = turbo_dmp_folders[0]
    tSlice = None

    # Add the tag and the run name
    if postProcessLinProfiles:
        # Tag is already present in the dict:
        _ = profileRunOptions["additional"].pop()

    profileRunOptions["additional"].append(('tag',theRunName,0))
    profileRunOptions["BOUT_run_name"] = theRunName
    # Create the runner
    profileRun = PBS_runner(**profileRunOptions)
    # Execute
    _, _ = profileRun.execute_runs(\
                                 remove_old               = remove_old,\
                                 post_processing_function = curPostProcessor,\
                                 # Declare dependencies
                                 job_dependencies = PBS_ids,\
                                 # This function will be called every time after
                                 # performing a run
                                 post_process_after_every_run = True,\
                                 # Below are the kwargs arguments being passed to
                                 # the post processing function
                                 # Switches
                                 driverName     = "parDriver"   ,\
                                 tSlice         = tSlice        ,\
                                 theRunName     = theRunName    ,\
                                 # Below are the kwargs given to the
                                 # restartFromFunc
                                 aScanPath      = aScanPath     ,\
                                 scanParameters = scanParameters,\
                                 # Common kwargs
                                 **fieldPlotterKwargs           ,\
                                )
#}}}
#}}}

#{{{Growth rates (run this driver after all, as we need the collectionFolders)
if postProcessGrowthRates:
    scanParam  = scanParameters[0]
    theRunName = "a1-KiwiFlatZ-growthRates"
    curPostProcessor = postBoutRunner

    # Make a list of list, where each sublist will be used as the paths
    # in collectiveCollect
    collectionFolders = list(zip(linear_dmp_folders, turbo_dmp_folders))

    _, _ = linearRun.execute_runs(\
                                 remove_old               = remove_old,\
                                 post_processing_function = curPostProcessor,\
                                 # This function will be called every time after
                                 # performing a run
                                 post_process_after_every_run = False,\
                                 # Below are the kwargs arguments being passed to
                                 # the post processing function
                                 # Switches
                                 driverName       = "plotGrowthRates"  ,\
                                 # PostProcessDriver input
                                 **commonPlotterKwargs                 ,\
                                 theRunName       = theRunName         ,\
                                 # StatsAndSignalsDrivers input
                                 paths            = collectionFolders  ,\
                                 # DriversProbes input
                                 var              = var                  ,\
                                 scanParam        = scanParam            ,\
                                 yInd             = yInd                 ,\
                                 nProbes          = nProbes              ,\
                                 steadyStatePaths = expand_dmp_folders   ,\
                                 maxMode          = maxMode              ,\
                                 # Below are the kwargs given to the
                                 # restartFromFunc
                                 aScanPath      = aScanPath     ,\
                                 scanParameters = scanParameters,\
                                )
#}}}

#{{{Probes and energy (run this driver after all, as we need the collectionFolders)
if postProcessProbesAndEnergy:
    collectionFolders = [linear_dmp_folders[0],\
                         turbo_dmp_folders[0]]
    theRunName = "a1-KiwiFlatZ-all-energyProbesPlot"
    curPostProcessor = postBoutRunner

    # Found from the overshoot at the energy plot
    # Overshoot happening around timestep 4400, timestep 4700 looks ok
    # Init + expand = 4100 => 4700 - 4100 = 600
    tIndSaturatedTurb = 600

    _, _ = turboRun.execute_runs(\
                                 remove_old               = remove_old,\
                                 post_processing_function = curPostProcessor,\
                                 # Declare dependencies
                                 job_dependencies = PBS_ids,\
                                 # This function will be called every time after
                                 # performing a run
                                 post_process_after_every_run = True,\
                                 # Below are the kwargs arguments being passed to
                                 # the post processing function
                                 # postBoutRunner option
                                 driverName = "plotEnergyAndProbes"    ,\
                                 # PostProcessDriver input
                                 **commonPlotterKwargs                 ,\
                                 theRunName        = theRunName        ,\
                                 # StatsAndSignalsDrivers input
                                 paths             = collectionFolders,\
                                 # DriversProbes input
                                 var               = var              ,\
                                 yInd              = yInd             ,\
                                 nProbes           = nProbes          ,\
                                 maxMode           = maxMode          ,\
                                 tIndSaturatedTurb = tIndSaturatedTurb,\
                                 # The steady state path will be
                                 # converted using convertToCurrentScanParameters
                                 steadyStatePath = expand_dmp_folders[0],\
                                 # Below are the kwargs given to the
                                 # restartFromFunc and convertToCurrentScanParameters
                                 aScanPath      = aScanPath     ,\
                                 scanParameters = scanParameters,\
                                )
#}}}

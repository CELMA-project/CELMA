#!/usr/bin/env python

"""Tests that all drivers are working."""

import os, sys
# If we add to sys.path, then it must be an absolute path
commonDir = os.path.abspath('./../common')
# Sys path is a list of system paths
sys.path.append(commonDir)

import re
import numpy as np
from boutdata import collect
from bout_runners import basic_runner
from CELMAPython.drivers import postBoutRunner
from CELMAPython.plotHelpers import findLargestRadialGrad

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
postProcessInit = False
postProcessExp  = False
postProcessLin  = False
postProcessTrub = False
# Extra post-processors
postProcessLinProfiles     = False
postProcessGrowthRates     = True
postProcessTurbProfiles    = False
postProcessProbesAndEnergy = False

#{{{Main options
#{{{The scan
B0 = [1.0e-1  , 9.0e-2  ]
Lx = [4.8296  , 4.3466  ]
Ly = [270.4579, 243.4121]
scanParameters  = ["B0", "Lx", "Ly"]
series_add = [\
              ('input', 'B0', B0),\
              ('geom' , 'Lx', Lx),\
              ('geom' , 'Ly', Ly),\
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
useSubProcess          = False
#}}}
#{{{File handeling options
remove_old = False
directory  = "0-scanTest"
make       = False
cpy_source = True
#}}}
nproc = 4
#}}}

#{{{Abbrevations
#{{{Dictionaries with common runner options
commonRunnerKwargs =\
        {\
         "directory"          : directory         ,\
         "make"               : make              ,\
         "nproc"              : nproc             ,\
         "cpy_source"         : cpy_source        ,\
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
theRunName = "0-TestScan-0-initialize"
# Set the spatial domain
nz = 1
# Set the temporal domain
restart    = None
timestep   = [1e-10]
nout       = [2]
# Filter
ownFilterType = "none"
# Post processing option
tSlice     = slice(-2, None)
#}}}
#{{{Run and post processing
initRunner = basic_runner(\
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
                             ],\
                series_add = series_add                      ,\
                # Common options
                **commonRunnerKwargs                         ,\
                )

init_dmp_folders, _ = initRunner.execute_runs(\
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
timestep   = [1e-10]
nout       = [2]
# Filter
ownFilterType = "none"
# From previous outputs
aScanPath = init_dmp_folders[0]
# Name
theRunName = "0-TestScan-1-expand"
# Post processing option
tSlice     = slice(-2, None)
#}}}
#{{{Run and post processing
expandRunner = basic_runner(\
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
                             ],\
                series_add = series_add                      ,\
                # Common options
                **commonRunnerKwargs                         ,\
                )

expand_dmp_folders, _ = expandRunner.execute_runs(\
                               remove_old               = remove_old  ,\
                               post_processing_function = curPostProcessor,\
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
    noutProfile                 = 0
    timestepProfile             = 1e-10
    useHyperViscAzVortDProfile  = True
    saveTermsProfile            = True

    # Create a new runner as we would like to save all the fields
    profileRun = basic_runner(\
                  # Set temporal domain
                  nout         = noutProfile    ,\
                  timestep     = timestepProfile,\
                  # Set restart options
                  restart      = restart        ,\
                  restart_from = restartFromFunc,\
                  # Set additional options
                  additional = [
                                ('tag',theRunName,0),\
                                ('switch'      , 'useHyperViscAzVortD',useHyperViscAzVortDProfile),\
                                ('switch'      , 'saveTerms'          ,saveTermsProfile),\
                               ],\
                  series_add = series_add                      ,\
                  # Common options
                  **commonRunnerKwargs                         ,\
                  )
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
timestep  = [1e-10]
nout     = [10]
# Name
theRunName = "0-TestScan-2-linearPhase1"
# Post processing options
tSlice           = slice(-2, None)
varyMaxMin       = True
subPolAvg        = True
mode             = "perpAndPol"
#}}}
#{{{Run and post processing
linearRun = basic_runner(\
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
            # Common options
            **commonRunnerKwargs                         ,\
            )

linear_dmp_folders, _ = linearRun.execute_runs(\
                                 remove_old               = remove_old,\
                                 post_processing_function = curPostProcessor,\
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
#{{{ If growth rates are to be plotted
if postProcessGrowthRates:
    scanParam  = "B0"
    theRunName = "0-TestScan-all-growthRates"
    curPostProcessor = postBoutRunner

    _, _ = linearRun.execute_runs(\
                                 remove_old               = remove_old,\
                                 post_processing_function = curPostProcessor,\
                                 # This function will be called every time after
                                 # performing a run
                                 post_process_after_every_run = False,\
                                 # Below are the kwargs arguments being passed to
                                 # the post processing function
                                 # Switches
                                 driverName       = "plotGrowthRates",\
                                 # PostProcessDriver input
                                 **commonPlotterKwargs                 ,\
                                 theRunName        = theRunName        ,\
                                 # StatsAndSignalsDrivers input
                                 paths            = linear_dmp_folders,\
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
#{{{ If linear profiles are to be plotted
if postProcessLinProfiles:
    curPostProcessor = postBoutRunner
    theRunName = "0-TestScan-linearPhaseParProfiles"
    aScanPathProfiles = linear_dmp_folders[0]

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
if postProcessTrub:
    curPostProcessor = postBoutRunner
else:
    curPostProcessor = None
#{{{Turbulence options
# Switches
saveTerms           = False
useHyperViscAzVortD = [True]
# Set the temporal domain
nout     = [2]
timestep = [1e-10]
# Name
theRunName = "0-TestScan-3-turbulentPhase1"
# Post processing options
tSlice    = slice(-2, None)
aScanPath = linear_dmp_folders[0]
#}}}
#{{{Run and post processing
turboRun = basic_runner(\
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
                # Common options
                **commonRunnerKwargs                         ,\
                )

turbo_dmp_folders, _ = turboRun.execute_runs(\
                                 remove_old               = remove_old,\
                                 post_processing_function = curPostProcessor,\
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
    theRunName = "0-TestScan-turboPhaseParProfiles"
    aScanPathProfiles = turbo_dmp_folders[0]

    _, _ = turboRun.execute_runs(\
                                 remove_old               = remove_old      ,\
                                 post_processing_function = curPostProcessor,\
                                 # Declare dependencies
                                 job_dependencies = PBS_ids         ,\
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

#{{{ Probes and energy (run this driver after all, as we need the collectionFolders)
if postProcessProbesAndEnergy:
    collectionFolders = [linear_dmp_folders[0],\
                         turbo_dmp_folders[0]]
    theRunName = "0-TestScan-all-energyProbesPlot"
    curPostProcessor = postBoutRunner

    # Found from the overshoot at the energy plot
    # Here just set to None
    tIndSaturatedTurb = None

    _, _ = turboRun.execute_runs(\
                                 remove_old               = remove_old,\
                                 post_processing_function = curPostProcessor,\
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
                                 var             = var                  ,\
                                 yInd            = yInd                 ,\
                                 nProbes         = nProbes              ,\
                                 maxMode         = maxMode              ,\
                                 # The steady state path will be
                                 # converted using convertToCurrentScanParameters
                                 steadyStatePath = expand_dmp_folders[0],\
                                 # Below are the kwargs given to the
                                 # restartFromFunc and convertToCurrentScanParameters
                                 aScanPath      = aScanPath     ,\
                                 scanParameters = scanParameters,\
                                )
#}}}

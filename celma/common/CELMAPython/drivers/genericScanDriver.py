#!/usr/bin/env python

"""
Contains the restartFromFunc and GenericScanDriver
"""

from .postBoutRunner import postBoutRunner
from bout_runners import PBS_runner
import re
import inspect

# FIXME: Could make a setOption for the post processing of each type as
#        well
# FIXME: Could also add kwargs in runner options

#{{{restartFromFunc
def restartFromFunc(dmp_folder     = None,\
                    aScanPath      = None,\
                    scanParameters = None,\
                    **kwargs):
    #{{{docstring
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
    scanParameters : sequence except str
        Sequence of strings of the names of the scan paramemters.
    **kwargs : keyword dictionary
        Dictionary with additional keyword arguments, given by
        bout_runners.

    Returns
    -------
    scanPath : str
        aScanPath converted to the scan parameters of the current run.
    """
    #}}}

    # Make a template string of aScanPath
    scanPathTemplate = aScanPath
    if type(scanParameters) == str:
        message = ("restartFromFunc was given the string '{}' as an "
                   "input parameter when a non-string sequence was required").\
                    format(scanParameters)
        raise ValueError(message)
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

#{{{GenericScanDriver
class GenericScanDriver(object):
    #{{{docstring
    """
    Generic scan driver class.

    Constructor sets default options, which can be changed through the
    class' setters.
    """
    #}}}

    #{{{Constructor
    def __init__(self):
        #{{{docstring
        """
        Sets the default options.
        """
        #}}}
        # Set warning flags
        self._calledFunctions = {
                    "mainOptions"          : False,\
                    "postProcessingFlag"   : False,\
                    "commonRunnerOptions"  : False,\
                    "commonPlotterOptions" : False,\
                    "fieldPlotterOptions"  : False,\
                    "probePlotterOptions"  : False,\
                    "initOptions"          : False,\
                    "expandOptions"        : False,\
                    "linearOptions"        : False,\
                    "turbulenceOptions"    : False,\
                }

        # Set default parameters
        self._postProcessingFlags = {\
                "justPostProcess"            : False,\
                "postProcessInit"            : False,\
                "postProcessExp"             : False,\
                "postProcessLin"             : False,\
                "postProcessTurb"            : False,\
                "postProcessLinProfiles"     : False,\
                "postProcessTurbProfiles"    : False,\
                "postProcessProbesAndEnergy" : False,\
                "postProcessGrowthRates"     : False,\
                }

        self._commonRunnerOptions =\
                {\
                 "nproc"              : 48   ,\
                 "cpy_source"         : True ,\
                 "BOUT_nodes"         : 3    ,\
                 "BOUT_ppn"           : 16   ,\
                 "post_process_nproc" : 1    ,\
                 "post_process_nodes" : 1    ,\
                 "post_process_ppn"   : 20   ,\
                }

        self._commonPlotterOptions =\
                {\
                 "saveFolderFunc"   : "scanWTagSaveFunc",\
                 "convertToPhysical": False             ,\
                 "showPlot"         : False             ,\
                 "savePlot"         : True              ,\
                 "extension"        : "png"             ,\
                 "useSubProcess"    : True              ,\
                }

        self._fieldPlotterOptions =\
                {\
                 "xguards"          : False,\
                 "yguards"          : False,\
                 "xSlice"           : 0    ,\
                 "ySlice"           : 16   ,\
                 "zSlice"           : 0    ,\
                 "axisEqualParallel": False,\
                }

        self._probesPlotterOptions =\
                {\
                 "nProbes" : 5 ,\
                 "maxMode" : 10,\
                 "yInd"    : 16,\
                }
        #}}}

    #{{{setMainOptions
    def setMainOptions(self                         ,\
                       directory                    ,\
                       scanParameters               ,\
                       series_add                   ,\
                       theRunName                   ,\
                       make                  = False,\
                       varName               = "n"  ,\
                       pltName               = "n"  ,\
                       timeStepMultiplicator = 1    ,\
                       boutRunnersNoise      = None ,\
                      ):
        #{{{docstring
        """
        Sets the main options for the scan.

        Parameters
        ----------
        directory : str
            Path to BOUT.inp
        scanParameters : sequence of strings
            Sequence of all quantities that will change during a scan
        series_add : sequence of sequence which is not string
            The series_add for the scan
        theRunName : str
            Name of the run
        make : bool
            Shall the program be made or not
        varName : str
            Name to be collected
        pltName : str
            Name to be plotted
        timeStepMultiplicator : int
            How much the default time step should be multiplied with
        boutRunnersNoise : [None|float]
            Should the noise generator from bout runners be used in the
            linear runs. If this is None, the noise will be generated
            from the noise generator in noiseGenerator.cxx
        """
        #}}}

        if boutRunnersNoise is None:
            boutRunnersNoise = False

        self._calledFunctions["mainOptions"] = True

        self._directory             = directory
        self._scanParameters        = scanParameters
        self._series_add            = series_add
        self._theRunName            = theRunName
        self._varName               = varName
        self._make                  = make
        self._var                   = varName
        self._pltName               = pltName
        self._timeStepMultiplicator = timeStepMultiplicator
        self._boutRunnersNoise      = boutRunnersNoise
    #}}}

    #{{{setPostProcessingFlags
    def setPostProcessingFlags(self                              ,\
                               justPostProcess            = False,\
                               postProcessInit            = False,\
                               postProcessExp             = False,\
                               postProcessLin             = False,\
                               postProcessTurb            = False,\
                               postProcessLinProfiles     = False,\
                               postProcessTurbProfiles    = False,\
                               postProcessProbesAndEnergy = False,\
                               postProcessGrowthRates     = False,\
                               tIndSaturatedTurb          = None ,\
            ):
        #{{{docstring
        """
        Sets the postProcessing flags.

        Parameters
        ----------
        justPostProcess : bool
            If the driver should only do the post-processing.
            NOTE: justPostProcess sets the restart parameter. Do **NOT**
            set this to True the first run, as this will prevent proper
            restarting.
        postProcessInit : bool
            If plots for the initialization should be run
        postProcessExp : bool
            If plots for the expand should be run
        postProcessLin : bool
            If plots for the linear phase should be run
        postProcessTurb : bool
            If plots for the turbulent phase should be run
        postProcessLinProfiles : bool
            If plots for all fields should be plotted after the linear phase
        postProcessTurbProfiles : bool
            If plots for all fields should be plotted after the turbulent phase
        postProcessProbesAndEnergy : bool
            If the probe data and the energy should be plotted
        postProcessGrowthRates : bool
            If the growth rates should be plotted
        tIndSaturatedTurb : int
            The index in the turbulence run where the turbulence is
            saturated (usually taken to be after the overshoot in the
            kinetic energy)
        """
        #}}}

        self._calledFunctions["postProcessingFlag"] = True

        self._postProcessingFlags = {\
                "justPostProcess"            : justPostProcess           ,\
                "postProcessInit"            : postProcessInit           ,\
                "postProcessExp"             : postProcessExp            ,\
                "postProcessLin"             : postProcessLin            ,\
                "postProcessTurb"            : postProcessTurb           ,\
                "postProcessLinProfiles"     : postProcessLinProfiles    ,\
                "postProcessTurbProfiles"    : postProcessTurbProfiles   ,\
                "postProcessProbesAndEnergy" : postProcessProbesAndEnergy,\
                "postProcessGrowthRates"     : postProcessGrowthRates    ,\
                "tIndSaturatedTurb"          : tIndSaturatedTurb         ,\
                }
    #}}}

    #{{{setCommonRunnerOptions
    def setCommonRunnerOptions(self                     ,\
                               nproc              = 48  ,\
                               cpy_source         = True,\
                               BOUT_nodes         = 3   ,\
                               BOUT_ppn           = 16  ,\
                               post_process_nproc = 1   ,\
                               post_process_nodes = 1   ,\
                               post_process_ppn   = 20  ,\
            ):
        #{{{docstring
        """
        Sets the kwargs for the common runner.

        Parameters
        ----------
        nproc : int
            Number of processors
        cpy_source : bool
            If the source should be copied
        BOUT_nodes : int
            How many nodes to run on
        BOUT_ppn : int
            How many processors per node
        post_process_nproc : int
            How many processors for post processor
        post_process_nodes : int
            How many nodes for post processor
        post_process_ppn : int
            How many processors per node for post processor
        """
        #}}}

        self._calledFunctions["commonRunnerOptions"] = True

        self._commonRunnerOptions =\
                {\
                 "nproc"              : nproc             ,\
                 "cpy_source"         : cpy_source        ,\
                 "BOUT_nodes"         : BOUT_nodes        ,\
                 "BOUT_ppn"           : BOUT_ppn          ,\
                 "post_process_nproc" : post_process_nproc,\
                 "post_process_nodes" : post_process_nodes,\
                 "post_process_ppn"   : post_process_ppn  ,\
                }
    #}}}

    #{{{setCommonPlotterOptions
    def setCommonPlotterOptions(self                                 ,\
                               saveFolderFunc    = "scanWTagSaveFunc",\
                               convertToPhysical = False             ,\
                               showPlot          = False             ,\
                               savePlot          = True              ,\
                               extension         = "png"             ,\
                               useSubProcess     = True              ,\
            ):
        #{{{docstring
        """
        Sets the kwargs for common plots.

        Parameters
        ----------
        saveFolderFunc : str
           What svae folder function to use
        convertToPhysical : bool
            Whether or not to convert to physical units
        showPlot : bool
            If the plot should be displayed
        savePlot : bool
            If the plot should be saved
        extension : str
            Extension of the file when not animating
        useSubProcess : bool
            If sub process should be used for the plotting
        """
        #}}}

        self._calledFunctions["commonPlotterOptions"] = True

        self._commonPlotterOptions =\
                {\
                 "saveFolderFunc"   : saveFolderFunc   ,\
                 "convertToPhysical": convertToPhysical,\
                 "showPlot"         : showPlot         ,\
                 "savePlot"         : savePlot         ,\
                 "extension"        : extension        ,\
                 "useSubProcess"    : useSubProcess    ,\
                }
        #}}}

    #{{{setFieldPlottersOptions
    def setFieldPlottersOptions(self                    ,\
                               xguards           = False,\
                               yguards           = False,\
                               xSlice            = 0    ,\
                               ySlice            = 16   ,\
                               zSlice            = 0    ,\
                               axisEqualParallel = False,\
            ):
        #{{{docstring
        """
        Sets the kwargs for field plotters.

        Parameters
        ----------
        xguards : bool
            Should the xguards be collected and plotted?
        yguards : bool
            Should the xguards be collected and plotted?
        xSlice : int
            Constant of x when slicing is not taking place
        ySlice : int
            Constant of y when slicing is not taking place
        zSlice : int
            Constant of z when slicing is not taking place
        axisEqualParallel : bool
            If axis equal should be used on the parallel plots
        """
        #}}}

        self._calledFunctions["fieldPlotterOptions"] = True

        self._fieldPlotterOptions =\
                {\
                 "xguards"          : xguards          ,\
                 "yguards"          : yguards          ,\
                 "xSlice"           : xSlice           ,\
                 "ySlice"           : ySlice           ,\
                 "zSlice"           : zSlice           ,\
                 "axisEqualParallel": axisEqualParallel,\
                }
    #}}}

    #{{{setProbePlottersOptions
    def setProbePlottersOptions(self        ,\
                               nProbes = 5  ,\
                               maxMode = 10 ,\
                               yInd    = 16 ,\
            ):
        #{{{docstring
        """
        Sets the kwargs for probe plots.

        Parameters
        ----------
        nProbes : int
            Number of probes to be used
        maxMode : int
            Maximum mode number to investigate
        yInd : int
            y index to investigate
        """
        #}}}

        self._calledFunctions["probePlotterOptions"] = True

        self._probesPlotterOptions =\
                {\
                 "nProbes" : nProbes,\
                 "maxMode" : maxMode,\
                 "yInd"    : yInd   ,\
                }
    #}}}

    #{{{setInitOptions
    def setInitOptions(self                              ,\
                       timestep              = 2e3       ,\
                       nout                  = 2         ,\
                       BOUT_walltime         = "05:00:00",\
                       post_process_walltime = "0:29:00" ,\
                       post_process_queue    = "xpresq"  ,\
            ):
        #{{{docstring
        """
        Sets the kwargs for the init run.

        Parameters
        ----------
        NOTE: timestep will be multiplied with timeStepMultiplicator
        For details, see bout_runners input
        """
        #}}}
        if self._calledFunctions["mainOptions"] == False:
            message =\
                "self.setMainOptions must be called prior to setInitOptions"
            raise ValueError(message)

        # Check where this function is called from
        curframe = inspect.currentframe()
        calframe = inspect.getouterframes(curframe, 2)
        if calframe[1][3] != "__init__":
            self._calledFunctions["setInitOptions"] = True

        self._initOptions =\
                {\
                 "timestep" : (timestep*self._timeStepMultiplicator,),\
                 "nout"     : (nout,),\
                }

        self._initPBSOptions =\
                {\
                 "BOUT_walltime"         : BOUT_walltime        ,\
                 "post_process_walltime" : post_process_walltime,\
                 "post_process_queue"    : post_process_queue   ,\
                }

    #}}}

    #{{{setExpandOptions
    def setExpandOptions(self                              ,\
                         nz                    = 256       ,\
                         timestep              = 50        ,\
                         nout                  = 2         ,\
                         BOUT_walltime         = "24:00:00",\
                         post_process_walltime = "0:29:00" ,\
                         post_process_queue    = "xpresq"  ,\
            ):
        #{{{docstring
        """
        Sets the kwargs for the expand run.

        Parameters
        ----------
        NOTE: timestep will be multiplied with timeStepMultiplicator
        For details, see bout_runners input
        """
        #}}}
        if self._calledFunctions["mainOptions"] == False:
            message =\
                "self.setMainOptions must be called prior to setExpandOptions"
            raise ValueError(message)

        # Check where this function is called from
        curframe = inspect.currentframe()
        calframe = inspect.getouterframes(curframe, 2)
        if calframe[1][3] != "__init__":
            self._calledFunctions["setExpandOptions"] = True

        self._expandOptions =\
                {\
                 "nz"       : nz,\
                 "timestep" : (timestep*self._timeStepMultiplicator,),\
                 "nout"     : (nout,),\
                }

        self._expandPBSOptions =\
                {\
                 "BOUT_walltime"         : BOUT_walltime        ,\
                 "post_process_walltime" : post_process_walltime,\
                 "post_process_queue"    : post_process_queue   ,\
                }

    #}}}

    #{{{setLinearOptions
    def setLinearOptions(self                              ,\
                         timestep              = 1         ,\
                         nout                  = 1000      ,\
                         BOUT_walltime         = '72:00:00',\
                         post_process_walltime = '03:00:00',\
                         post_process_queue    = 'workq'   ,\
                        ):
        #{{{docstring
        """
        Sets the kwargs for the linear run.

        Parameters
        ----------
        NOTE: timestep will be multiplied with timeStepMultiplicator
        For details, see bout_runners input
        """
        #}}}
        if self._calledFunctions["mainOptions"] == False:
            message =\
                "self.setMainOptions must be called prior to setLinearOptions"
            raise ValueError(message)

        # Check where this function is called from
        curframe = inspect.currentframe()
        calframe = inspect.getouterframes(curframe, 2)
        if calframe[1][3] != "__init__":
            self._calledFunctions["setLinearOptions"] = True

        self._linearOptions =\
                {\
                 "timestep" : (timestep*self._timeStepMultiplicator,),\
                 "nout"     : (nout,),\
                }

        self._linearPBSOptions =\
                {\
                 "BOUT_walltime"         : BOUT_walltime        ,\
                 "post_process_walltime" : post_process_walltime,\
                 "post_process_queue"    : post_process_queue   ,\
                }

    #}}}

    #{{{setTurbulenceOptions
    def setTurbulenceOptions(self                              ,\
                             timestep              = 1         ,\
                             nout                  = 5000      ,\
                             BOUT_walltime         = "72:00:00",\
                             post_process_walltime = "03:00:00",\
                             post_process_queue    = "workq"   ,\
            ):
        #{{{docstring
        """
        Sets the kwargs for the turbulence run.

        Parameters
        ----------
        NOTE: timestep will be multiplied with timeStepMultiplicator
        For details, see bout_runners input
        """
        #}}}
        if self._calledFunctions["mainOptions"] == False:
            message =\
                "self.setMainOptions must be called prior to setTurbulenceOptions"
            raise ValueError(message)

        # Check where this function is called from
        curframe = inspect.currentframe()
        calframe = inspect.getouterframes(curframe, 2)
        if calframe[1][3] != "__init__":
            self._calledFunctions["setTurbulenceOptions"] = True

        self._turbulenceOptions =\
                {\
                 "timestep" : (timestep*self._timeStepMultiplicator,),\
                 "nout"     : (nout,),\
                }

        self._turbulencePBSOptions =\
                {\
                 "BOUT_walltime"         : BOUT_walltime        ,\
                 "post_process_walltime" : post_process_walltime,\
                 "post_process_queue"    : post_process_queue   ,\
                }
    #}}}

    #{{{runScan
    def runScan(self):
        """
        Calls the drivers used for running scans
        """

        if not(self._calledFunctions["mainOptions"]):
            message = "self.setMainOptions must be called prior to a run"
            raise ValueError(message)

        keysToBeCalled = []

        for key, val in self._calledFunctions.items():
            if val:
                print("{:<25} has been called".format(key))
            else:
                keysToBeCalled.append(key)
                print("{:<25} has NOT been called (using defaults)".format(key))

        # Set default runner options if not set
        for key in keysToBeCalled:
            if key == "initOptions":
                self.setInitOptions()
            elif key == "expandOptions":
                self.setExpandOptions()
            elif key == "linearOptions":
                self.setLinearOptions()
            elif key == "turbulenceOptions":
                self.setTurbulenceOptions()


        # Make dictionary to variables
        for (flag, value) in self._postProcessingFlags.items():
            setattr(self, flag, value)

        # Update dicts
        self._fieldPlotterOptions.update(self._commonPlotterOptions)
        self._commonRunnerOptions["directory"] = self._directory

        #{{{Init runner
        if self.postProcessInit:
            curPostProcessor = postBoutRunner
        else:
            curPostProcessor = None
        #{{{Init options
        # Name
        theRunName = self._theRunName + "-0-initialize"
        # Set the spatial domain
        nz = 1
        # Set the temporal domain
        restart    = None
        # Filter
        ownFilterType = "none"
        #Switches
        useHyperViscAzVortD = (False,)
        # Post processing option
        tSlice = None
        # PBS options
        BOUT_run_name         = theRunName
        post_process_run_name = 'post' + theRunName.capitalize()
        #}}}

        #{{{Run and post processing
        initRunner = PBS_runner(\
            # Shall we make
            make       = self._make,\
            # Set spatial domain
            nz         = nz        ,\
            # Set the restart option
            restart    = restart   ,\
            # Init options
            **self._initOptions    ,\
            # Set additional option
            additional = (
                ('tag',theRunName,0),\
                ('ownFilters'  , 'type', ownFilterType),\
                ('switch'      , 'useHyperViscAzVortD', useHyperViscAzVortD),\
                         ),\
            series_add = self._series_add,\
            # PBS options
            BOUT_run_name         = BOUT_run_name        ,\
            post_process_run_name = post_process_run_name,\
            **self._initPBSOptions                       ,\
            # Common options
            **self._commonRunnerOptions                  ,\
                        )

        init_dmp_folderss, PBS_ids = initRunner.execute_runs(\
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
            **self._fieldPlotterOptions            ,\
                                    )
        #}}}
        #}}}

        if self.justPostProcess:
            restart = None
        else:
            restart = "overwrite"

        #{{{The expand runner
        if self.postProcessExp:
            curPostProcessor = postBoutRunner
        else:
            curPostProcessor = None
        #{{{Expand runner options
        # Set the spatial domain
        # Filter
        ownFilterType = "none"
        #Switches
        useHyperViscAzVortD = (False,)
        # From previous outputs
        aScanPath = init_dmp_folderss[0]
        # Name
        theRunName = self._theRunName + "-1-expand"
        # PBS options
        BOUT_run_name         = theRunName
        post_process_run_name = 'post' + theRunName.capitalize()
        # Post processing option
        tSlice = None
        #}}}

        #{{{Run and post processing
        expandRunner = PBS_runner(\
            # Expand options
            **self._expandOptions         ,\
            # Set restart options
            restart      = restart        ,\
            restart_from = restartFromFunc,\
            # Set additional options
            additional = (
                ('tag',theRunName,0),\
                ('ownFilters'  , 'type', ownFilterType),\
                ('switch'      , 'useHyperViscAzVortD', useHyperViscAzVortD),\
                         ),\
            series_add = self._series_add                ,\
            # PBS options
            BOUT_run_name         = BOUT_run_name        ,\
            post_process_run_name = post_process_run_name,\
            **self._expandPBSOptions                     ,\
            # Common options
            **self._commonRunnerOptions                  ,\
                        )

        expand_dmp_folderss, PBS_ids = expandRunner.execute_runs(\
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
            aScanPath      = aScanPath            ,\
            scanParameters = self._scanParameters ,\
            # Common kwargs
            **self._fieldPlotterOptions           ,\
                                      )
        #}}}
        #}}}

        #{{{ If profiles are to be plotted
        if self.postProcessLinProfiles or self.postProcessTurbProfiles:
            noutProfile                 = 3
            timestepProfile             = 10*self._timeStepMultiplicator
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
                "additional" : (
                 ('switch', 'useHyperViscAzVortD',useHyperViscAzVortDProfile),\
                 ('switch', 'saveTerms'          ,saveTermsProfile),\
                             ),\
                "series_add" : self._series_add,\
                # Common options
                **self._commonRunnerOptions    ,\
                                }
        #}}}

        #{{{The linear runner
        if self.postProcessLin:
            curPostProcessor = postBoutRunner
        else:
            curPostProcessor = None
        #{{{ Linear options
        #Switches
        saveTerms           = False
        useHyperViscAzVortD = (True,)
        if not(self._boutRunnersNoise):
            includeNoise  = True
            forceAddNoise = True
            add_noise     = None
        else:
            includeNoise  = False
            forceAddNoise = False
            add_noise     = {"lnN":self._boutRunnersNoise}
        # As this is scan dependent, the driver finds the correct folder
        maxGradRhoFolder = expand_dmp_folderss[0]
        # From previous outputs
        aScanPath = expand_dmp_folderss[0]
        # Name
        theRunName = self._theRunName + "-2-linearPhase1"
        # PBS options
        BOUT_run_name         = theRunName
        post_process_run_name = 'post' + theRunName.capitalize()
        # Post processing options
        tSlice     = slice(0, None, 2)
        varyMaxMin = True
        subPolAvg  = True
        mode       = "perpAndPol"
        #}}}
        #{{{Run and post processing
        linearRun = PBS_runner(\
            # Linear options
            **self._linearOptions         ,\
            # Set restart options
            restart      = restart        ,\
            restart_from = restartFromFunc,\
            # Set additional options
            additional = (
                ('tag'   , theRunName           ,0),\
                ('switch', 'includeNoise'       , includeNoise ),\
                ('switch', 'forceAddNoise'      ,forceAddNoise),\
                ('switch', 'useHyperViscAzVortD',useHyperViscAzVortD),\
                ('switch', 'saveTerms'          ,saveTerms),\
                         ),\
            series_add = self._series_add                ,\
            # Set eventual noise
            add_noise             = add_noise            ,\
            # PBS options
            BOUT_run_name         = BOUT_run_name        ,\
            post_process_run_name = post_process_run_name,\
            **self._linearPBSOptions                     ,\
            # Common options
            **self._commonRunnerOptions                  ,\
                    )

        linear_dmp_folderss, PBS_ids = linearRun.execute_runs(\
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
                varName          = self._varName   ,\
                pltName          = self._pltName   ,\
                varyMaxMin       = varyMaxMin      ,\
                mode             = mode            ,\
                maxGradRhoFolder = maxGradRhoFolder,\
                # Below are the kwargs given to the
                # restartFromFunc
                aScanPath      = aScanPath         ,\
                scanParameters = self._scanParameters,\
                # Common kwargs
                **self._fieldPlotterOptions          ,\
                                        )
        #}}}
        #{{{ If linear profiles are to be plotted
        if self.postProcessLinProfiles:
            curPostProcessor = postBoutRunner
            theRunName = self._theRunName + "-2-linearPhaseParProfiles"
            tSlice = None

            # Add the tag and the run name
            profileRunOptions["additional"] =\
                    list(profileRunOptions["additional"])
            profileRunOptions["additional"].append(('tag',theRunName,0))
            profileRunOptions["additional"] =\
                    tuple(profileRunOptions["additional"])
            profileRunOptions["BOUT_run_name"] = theRunName

            # Create the runner
            profileRun = PBS_runner(**profileRunOptions)
            # Execute
            _, _ = profileRun.execute_runs(\
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
                scanParameters = self._scanParameters,\
                # Common kwargs
                **self._fieldPlotterOptions     ,\
                                        )
        #}}}
        #}}}

        #{{{Turbulence runner
        if self.postProcessTurb:
            curPostProcessor = postBoutRunner
        else:
            curPostProcessor = None
        #{{{Turbulence options
        # Switches
        saveTerms           = False
        useHyperViscAzVortD = (True,)
        # Name
        theRunName = self._theRunName + "-3-turbulentPhase1"
        # PBS options
        BOUT_run_name         = theRunName
        post_process_run_name = 'post' + theRunName.capitalize()
        # Post processing options
        tSlice    = slice(0, None, 10)
        aScanPath = linear_dmp_folderss[0]
        #}}}
        #{{{Run and post processing
        turboRun = PBS_runner(\
            # Set turbulence options
            **self._turbulenceOptions       ,\
            # Set restart options
            restart      = restart          ,\
            restart_from = restartFromFunc  ,\
            additional = (
                ('tag',theRunName,0),\
                ('switch'      , 'useHyperViscAzVortD',useHyperViscAzVortD),\
                ('switch'      , 'saveTerms'          ,saveTerms),\
                         ),\
            series_add = self._series_add                ,\
            # PBS options
            BOUT_run_name         = BOUT_run_name        ,\
            post_process_run_name = post_process_run_name,\
            **self._turbulencePBSOptions                 ,\
            # Common options
            **self._commonRunnerOptions                  ,\
                        )

        turbo_dmp_folderss, PBS_ids = turboRun.execute_runs(\
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
            varName        = self._varName   ,\
            pltName        = self._pltName   ,\
            # Below are the kwargs given to the
            # restartFromFunc
            aScanPath      = aScanPath           ,\
            scanParameters = self._scanParameters,\
            # Common kwargs
            **self._fieldPlotterOptions          ,\
                                        )
        #}}}
        #{{{ If linear profiles are to be plotted
        if self.postProcessTurbProfiles:
            curPostProcessor = postBoutRunner
            theRunName = self._theRunName + "-3-turbulentPhase1ParProfiles"
            tSlice = None

            # Add the tag and the run name
            if self.postProcessLinProfiles:
                # Tag is already present in the dict:
                profileRunOptions["additional"] =\
                        list(profileRunOptions["additional"])
                _ = profileRunOptions["additional"].pop()
                profileRunOptions["additional"] =\
                        tuple(profileRunOptions["additional"])

            profileRunOptions["additional"] = list(profileRunOptions["additional"])
            profileRunOptions["additional"].append(('tag',theRunName,0))
            profileRunOptions["additional"] = tuple(profileRunOptions["additional"])
            profileRunOptions["BOUT_run_name"] = theRunName
            # Create the runner
            profileRun = PBS_runner(**profileRunOptions)
            # Execute
            _, _ = profileRun.execute_runs(\
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
                scanParameters = self._scanParameters,\
                # Common kwargs
                **self._fieldPlotterOptions    ,\
                                        )
        #}}}
        #}}}

        #{{{Growth rates (run after all others, as collectionFolders is needed)
        if self.postProcessGrowthRates:
            scanParam  = self._scanParameters[0]
            theRunName = self._theRunName + "-growthRates"
            curPostProcessor = postBoutRunner

            # Make a tuple of tuples, where each subtuple will be used as the
            # paths in collectiveCollect
            collectionFolders =\
                tuple(zip(linear_dmp_folderss, turbo_dmp_folderss))

            _, _ = linearRun.execute_runs(\
                post_processing_function = curPostProcessor,\
                # This function will be called every time after
                # performing a run
                post_process_after_every_run = False,\
                # Below are the kwargs arguments being passed to
                # the post processing function
                # Switches
                driverName       = "plotGrowthRates"  ,\
                # PostProcessDriver input
                **self._commonPlotterOptions          ,\
                theRunName       = theRunName         ,\
                # StatsAndSignalsDrivers input
                paths            = collectionFolders  ,\
                # DriversProbes input
                var              = self._var          ,\
                scanParam        = scanParam          ,\
                steadyStatePaths = expand_dmp_folderss ,\
                **self._probesPlotterOptions          ,\
                # Below are the kwargs given to the
                # restartFromFunc
                aScanPath      = aScanPath            ,\
                scanParameters = self._scanParameters ,\
                                        )
        #}}}

        #{{{Probes and energy (run this driver after all, as we need the collectionFolders)
        if self.postProcessProbesAndEnergy:
            collectionFolders = (linear_dmp_folderss[0],\
                                 turbo_dmp_folderss[0])

            theRunName = self._theRunName + "-energyProbesPlot"
            curPostProcessor = postBoutRunner

            _, _ = turboRun.execute_runs(\
                post_processing_function = curPostProcessor,\
                # Declare dependencies
                job_dependencies = PBS_ids,\
                # This function will be called every time after
                # performing a run
                post_process_after_every_run = True,\
                # Below are the kwargs arguments being passed to
                # the post processing function
                # postBoutRunner option
                driverName = "plotEnergyAndProbes"         ,\
                # PostProcessDriver input
                **self._commonPlotterOptions               ,\
                theRunName        = theRunName             ,\
                # StatsAndSignalsDrivers input
                paths             = collectionFolders      ,\
                # DriversProbes input
                var               = self._var              ,\
                **self._probesPlotterOptions               ,\
                tIndSaturatedTurb = self.tIndSaturatedTurb ,\
                # The steady state path will be
                # converted using convertToCurrentScanParameters
                steadyStatePath = expand_dmp_folderss[0],\
                # Below are the kwargs given to the
                # restartFromFunc and convertToCurrentScanParameters
                aScanPath      = aScanPath           ,\
                scanParameters = self._scanParameters,\
                                        )
        #}}}
    #}}}
#}}}

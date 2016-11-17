#!/usr/bin/env python

"""
Contains the restartFromFunc and GenericScanDriver
"""

from .postBoutRunner import postBoutRunner
from bout_runners import PBS_runner
import re
import inspect
import os

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
                    "runOptions"           : False,\
                    "postProcessingFlag"   : False,\
                    "commonRunnerOptions"  : False,\
                    "commonPlotterOptions" : False,\
                    "fieldPlotterOptions"  : None ,\
                    "probePlotterOptions"  : None ,\
                    "initOptions"          : False,\
                    "initPostOptions"      : False,\
                    "expandOptions"        : False,\
                    "expandPostOptions"    : False,\
                    "linearOptions"        : False,\
                    "linearPostOptions"    : False,\
                    "turbulenceOptions"    : False,\
                    "turbulencePostOptions": False,\
                }
        #}}}

    #{{{setMainOptions
    def setMainOptions(self                         ,\
                       directory                    ,\
                       scanParameters               ,\
                       series_add                   ,\
                       theRunName                   ,\
                       make                  = False,\
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
        timeStepMultiplicator : int
            How much the default time step should be multiplied with
        boutRunnersNoise : [None|dict]
            If this is None, the noise will be generated
            from the noise generator in noiseGenerator.cxx
            If this is set as a dict, the dict must have the key of one
            of the variables evolved in time, and the value must be the
            amplitude of the noise.
        """
        #}}}

        if boutRunnersNoise is None:
            boutRunnersNoise = False

        self._calledFunctions["mainOptions"] = True

        self._directory             = directory
        self._scanParameters        = scanParameters
        self._series_add            = series_add
        self._theRunName            = theRunName
        self._make                  = make
        self._timeStepMultiplicator = timeStepMultiplicator
        self._boutRunnersNoise      = boutRunnersNoise
    #}}}

    #{{{setRunOptions
    def setRunOptions(self            ,\
                      runInit   = True,\
                      runExpand = True,\
                      runLin    = True,\
                      runTurb   = True,\
            ):
        #{{{docstring
        """
        Set which runs which will be simulated.

        Parameters
        ----------
        runInit : bool
            Inital run will be performed if True
        runExpand : bool
            Expand run will be performed if True
        runLin : bool
            Linear run will be performed if True
        runTurb : bool
            Turbulence run will be performed if True
        """
        #}}}
        # Check where this function is called from
        curframe = inspect.currentframe()
        calframe = inspect.getouterframes(curframe, 2)
        if calframe[1][3] != "__init__":
            self._calledFunctions["runOptions"] = True

        self._runOptions =\
                {\
                  "runInit"   : runInit  ,\
                  "runExpand" : runExpand,\
                  "runLin"    : runLin   ,\
                  "runTurb"   : runTurb  ,\
                }
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

        # Update the post options
        if hasattr(self, "_initPostOptions"):
            if self._initUsesFieldPlotter:
                self._initPostOptions.update(self._fieldPlotterOptions)
        if hasattr(self, "_expandPostOptions"):
            if self._expandUsesFieldPlotter:
                self._expandPostOptions.update(self._fieldPlotterOptions)
        if hasattr(self, "_linearPostOptions"):
            if self._linearUsesFieldPlotter:
                self._linearPostOptions.update(self._fieldPlotterOptions)
        if hasattr(self, "_turbulencePostOptions"):
            if self._turbulenceUsesFieldPlotter:
                self._turbulencePostOptions.update(self._fieldPlotterOptions)

    #}}}

    #{{{setProbePlottersOptions
    def setProbePlottersOptions(self        ,\
                               nProbes = 5  ,\
                               maxMode = 10 ,\
                               yInd    = 16 ,\
                               var     = "n",\
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
        var : str
            The variable to investigate
        """
        #}}}

        self._calledFunctions["probePlotterOptions"] = True

        self._probesPlotterOptions =\
                {\
                 "nProbes" : nProbes,\
                 "maxMode" : maxMode,\
                 "yInd"    : yInd   ,\
                 "var"     : var    ,\
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
            self._calledFunctions["initOptions"] = True

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

    #{{{setInitPostOptions
    def setInitPostOptions(self                         ,\
                           useDefault             = True,\
                           useFieldPlotterOptions = True,\
                           **kwargs):
        #{{{docstring
        """
        Sets the kwargs for the post processing of the init runs.

        Parameters
        ----------
        useDefault : bool
        useFieldPlotterOptions : bool
            Whether or not the fieldPlotterOptions should be used.
        kwargs : keyword arguments
            Do not use keyword arguments if useDefault is True.
            keywords to the post processors.
            Must include "driverName" if useDefault is False.
            For details more details, see the info about each post
            processing driver in common/CELMAPython/drivers
        """
        #}}}

        if useDefault and len(kwargs) !=0:
            message = ("No extra keywords should be given when "
                       "useDefault is set.\n"
                       "These keywords were found {}".format(kwargs))
            raise ValueError(message)

        if useDefault:
            # Make the field plotter
            self._initUsesFieldPlotter = True

            if self._calledFunctions["fieldPlotterOptions"] == None:
                self.setFieldPlottersOptions()
                self._calledFunctions["fieldPlotterOptions"] == False

            self._initPostOptions =\
                    {"driverName" : "plot1DAnd2DDriver",\
                     "tSlice"     : None      ,\
                     **self._fieldPlotterOptions
                    }

        else:
            if not("driverName" in kwargs.keys()):
                message = "'driverName' missing from keyword arguments"
                raise ValueError(message)

            # Check that field plotter is called first (if set)
            if useFieldPlotterOptions:
                self._initUsesFieldPlotter = True
                if self._calledFunctions["fieldPlotterOptions"] == None:
                    self.setFieldPlottersOptions()
                    self._calledFunctions["fieldPlotterOptions"] == False
            else:
                self._initUsesFieldPlotter = False

            self._initPostOptions = {**kwargs}

            if useFieldPlotterOptions:
                self._initPostOptions.update(self._fieldPlotterOptions)
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
            self._calledFunctions["expandOptions"] = True

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

    #{{{setExpandPostOptions
    def setExpandPostOptions(self,
                             useDefault             = True,\
                             useFieldPlotterOptions = True,\
                             **kwargs):
        #{{{docstring
        """
        Sets the kwargs for the post processing of the expand runs.

        Parameters
        ----------
        useDefault : bool
        useFieldPlotterOptions : bool
            Whether or not the fieldPlotterOptions should be used.
        kwargs : keyword arguments
            Do not use keyword arguments if useDefault is True.
            keywords to the post processors.
            Must include "driverName" if useDefault is False.
            For details more details, see the info about each post
            processing driver in common/CELMAPython/drivers
        """
        #}}}

        if useDefault and len(kwargs) !=0:
            message = ("No extra keywords should be given when "
                       "useDefault is set.\n"
                       "These keywords were found {}".format(kwargs))
            raise ValueError(message)

        if useDefault:
            # Make the field plotter
            self._expandUsesFieldPlotter = True

            if self._calledFunctions["fieldPlotterOptions"] == None:
                self.setFieldPlottersOptions()
                self._calledFunctions["fieldPlotterOptions"] == False

            self._expandPostOptions =\
                    {"driverName" : "plot1DAnd2DDriver",\
                     "tSlice"     : None      ,\
                     **self._fieldPlotterOptions
                    }

        else:
            if not("driverName" in kwargs.keys()):
                message = "'driverName' missing from keyword arguments"
                raise ValueError(message)

            # Check that field plotter is called first (if set)
            if useFieldPlotterOptions:
                self._expandUsesFieldPlotter = True
                if self._calledFunctions["fieldPlotterOptions"] == None:
                    self.setFieldPlottersOptions()
                    self._calledFunctions["fieldPlotterOptions"] == False
            else:
                self._expandUsesFieldPlotter = False

            self._expandPostOptions = {**kwargs}

            if useFieldPlotterOptions:
                self._expandPostOptions.update(self._fieldPlotterOptions)
    #}}}

    #{{{setLinearOptions
    def setLinearOptions(self                              ,\
                         timestep              = 1         ,\
                         nout                  = 1000      ,\
                         BOUT_walltime         = "72:00:00",\
                         post_process_walltime = "03:00:00",\
                         post_process_queue    = "workq"   ,\
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
            self._calledFunctions["linearOptions"] = True

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

    #{{{setLinearPostOptions
    def setLinearPostOptions(self,
                             useDefault             = True,\
                             useFieldPlotterOptions = True,\
                             useMaxGrad             = True,\
                             **kwargs):
        #{{{docstring
        """
        Sets the kwargs for the post processing of the linear runs.

        Parameters
        ----------
        useDefault : bool
        useFieldPlotterOptions : bool
            Whether or not the fieldPlotterOptions should be used.
        useMaxGrad : bool
            Whether or not to use the max gradient calculated from the
            expand runs.
        kwargs : keyword arguments
            Do not use keyword arguments if useDefault is True.
            keywords to the post processors.
            Must include "driverName" if useDefault is False.
            For details more details, see the info about each post
            processing driver in common/CELMAPython/drivers
        """
        #}}}

        # Check that field plotter is called first (if set)
        if useDefault and len(kwargs) !=0:
            message = ("No extra keywords should be given when "
                       "useDefault is set.\n"
                       "These keywords were found {}".format(kwargs))
            raise ValueError(message)

        if useDefault:
            self._linearUsesFieldPlotter       = True
            self._linearPostOptionsUsesMaxGrad = True

            # Make the field plotter
            if self._calledFunctions["fieldPlotterOptions"] == None:
                self.setFieldPlottersOptions()
                self._calledFunctions["fieldPlotterOptions"] == False

            self._linearPostOptions =\
                    {"driverName" : "single2DDriver" ,\
                     "varName"    : "n"              ,\
                     "pltName"    : "n"              ,\
                     "tSlice"     : slice(0, None, 2),\
                     "varyMaxMin" : True             ,\
                     "subPolAvg"  : True             ,\
                     "mode"       : "perpAndPol"     ,\
                     **self._fieldPlotterOptions     ,\
                    }

        else:
            if not("driverName" in kwargs.keys()):
                message = "'driverName' missing from keyword arguments"
                raise ValueError(message)

            # Check that field plotter is called first (if set)
            if useFieldPlotterOptions:
                self._linearUsesFieldPlotter = True
                if self._calledFunctions["fieldPlotterOptions"] == None:
                    self.setFieldPlottersOptions()
                    self._calledFunctions["fieldPlotterOptions"] == False
            else:
                self._linearUsesFieldPlotter = False

            self._linearPostOptionsUsesMaxGrad = True if useMaxGrad else False

            self._linearPostOptions = {**kwargs}

            if useFieldPlotterOptions:
                self._linearPostOptions.update(self._fieldPlotterOptions)
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
            self._calledFunctions["turbulenceOptions"] = True

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

    #{{{setTurbulencePostOptions
    def setTurbulencePostOptions(self,
                                 useDefault             = True,\
                                 useFieldPlotterOptions = True,\
                                 **kwargs):
        #{{{docstring
        """
        Sets the kwargs for the post processing of the turbulence runs.

        Parameters
        ----------
        useDefault : bool
        useFieldPlotterOptions : bool
            Whether or not the fieldPlotterOptions should be used.
        kwargs : keyword arguments
            Do not use keyword arguments if useDefault is True.
            keywords to the post processors.
            Must include "driverName" if useDefault is False.
            For details more details, see the info about each post
            processing driver in common/CELMAPython/drivers
        """
        #}}}

        if useDefault and len(kwargs) !=0:
            message = ("No extra keywords should be given when "
                       "useDefault is set.\n"
                       "These keywords were found {}".format(kwargs))
            raise ValueError(message)

        if useDefault:
            # Make the field plotter
            self._turbulenceUsesFieldPlotter = True

            if self._calledFunctions["fieldPlotterOptions"] == None:
                self.setFieldPlottersOptions()
                self._calledFunctions["fieldPlotterOptions"] == False

            self._turbulencePostOptions =\
                    {"driverName" : "single2DDriver"  ,\
                     "tSlice"     : slice(0, None, 10),\
                     "varName"    : "n"               ,\
                     "pltName"    : "n"               ,\
                     **self._fieldPlotterOptions
                    }

        else:
            if not("driverName" in kwargs.keys()):
                message = "'driverName' missing from keyword arguments"
                raise ValueError(message)

            # Check that field plotter is called first (if set)
            if useFieldPlotterOptions:
                self._turbulenceUsesFieldPlotter = True
                if self._calledFunctions["fieldPlotterOptions"] == None:
                    self.setFieldPlottersOptions()
                    self._calledFunctions["fieldPlotterOptions"] == False
            else:
                self._turbulenceUsesFieldPlotter = False

            self._turbulencePostOptions = {**kwargs}

            if useFieldPlotterOptions:
                self._turbulencePostOptions.update(self._fieldPlotterOptions)
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
            if val is not(None):
                if val:
                    print("{:<25} has been called".format(key))
                else:
                    keysToBeCalled.append(key)
                    print("{:<25} has NOT been called (using defaults)".format(key))

        # Set default runner options if not set
        for key in keysToBeCalled:
            # Call the function from their names
            getattr(self, "set{}".format(key[0].upper()+key[1:]))()

        # Make dictionary to variables
        for (flag, value) in self._postProcessingFlags.items():
            setattr(self, flag, value)
        for (flag, value) in self._runOptions.items():
            setattr(self, flag, value)

        # Update dicts
        self._fieldPlotterOptions  .update(self._commonPlotterOptions)
        self._initPostOptions      .update(self._commonPlotterOptions)
        self._expandPostOptions    .update(self._commonPlotterOptions)
        self._linearPostOptions    .update(self._commonPlotterOptions)
        self._turbulencePostOptions.update(self._commonPlotterOptions)
        self._commonRunnerOptions["directory"] = self._directory

        # Call the runners
        if self.runInit:
            self._callInitRunner()

        if self.justPostProcess:
            self._restart = None
        else:
            self._restart = "overwrite"

        if self.runExpand:
            self._callExpandRunner()

        if self.runLin:
            self._callLinearRunner()

        if self.runTurb:
            self._callTurboRunner()

        if self.postProcessGrowthRates:
            self._callPostProcessGrowthRates()

        if self.postProcessProbesAndEnergy:
            self._callPostProcessProbesAndEnergy()

        #{{{ If profiles are to be plotted
        # FIXME: Consider making this into a function
        if self.postProcessLinProfiles or self.postProcessTurbProfiles:
            noutProfile                 = 3
            timestepProfile             = 10*self._timeStepMultiplicator
            restartProfile              = "overwrite"
            useHyperViscAzVortDProfile  = True
            saveTermsProfile            = True

            # Create the options for the runners
            # Notice that we would like to save all the fields here
            self._profileRunOptions = {\
                # Set temporal domain
                "nout"         : noutProfile    ,\
                "timestep"     : timestepProfile,\
                # Set restart options
                "restart"      : restartProfile ,\
                "restart_from" : restartFromFunc,\
                # Set additional options
                "additional" : (
                 ("switch", "useHyperViscAzVortD",useHyperViscAzVortDProfile),\
                 ("switch", "saveTerms"          ,saveTermsProfile),\
                             ),\
                "series_add" : self._series_add,\
                # Common options
                **self._commonRunnerOptions    ,\
                                }
        #}}}

        if self.postProcessLinProfiles:
            self._callpostProcessLinProfiles()

        if self.postProcessTurbProfiles:
            self._callPostProcessTurbProfiles()
    #}}}

    #{{{_callInitRunner
    def _callInitRunner(self):
        """Executes the init run and the post processing"""

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
        self._restart = None
        # Filter
        ownFilterType = "none"
        #Switches
        useHyperViscAzVortD = (False,)
        # PBS options
        BOUT_run_name         = theRunName
        post_process_run_name = "post" + theRunName.capitalize()
        #}}}

        #{{{Run and post processing
        initRunner = PBS_runner(\
            # Shall we make
            make       = self._make,\
            # Set spatial domain
            nz         = nz        ,\
            # Set the restart option
            restart    = self._restart,\
            # Init options
            **self._initOptions       ,\
            # Set additional option
            additional = (
                ("tag",theRunName,0),\
                ("ownFilters"  , "type", ownFilterType),\
                ("switch"      , "useHyperViscAzVortD", useHyperViscAzVortD),\
                         ),\
            series_add = self._series_add,\
            # PBS options
            BOUT_run_name         = BOUT_run_name        ,\
            post_process_run_name = post_process_run_name,\
            **self._initPBSOptions                       ,\
            # Common options
            **self._commonRunnerOptions                  ,\
                        )

        self._init_dmp_folders, self._init_PBS_ids = initRunner.execute_runs(\
            post_processing_function = curPostProcessor,\
            # This function will be called every time after
            # performing a run
            post_process_after_every_run = True,\
            # Below are the kwargs arguments being passed to
            # the post processing function
            theRunName  = theRunName,\
            **self._initPostOptions ,\
            )
        #}}}
    #}}}

    #{{{_callExpandRunner
    def _callExpandRunner(self):
        """Executes the expand run and the post processing"""

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
        self._initAScanPath = self._init_dmp_folders[0]
        # Name
        theRunName = self._theRunName + "-1-expand"
        # PBS options
        BOUT_run_name         = theRunName
        post_process_run_name = "post" + theRunName.capitalize()
        #}}}

        #{{{Run and post processing
        expandRunner = PBS_runner(\
            # Expand options
            **self._expandOptions         ,\
            # Set restart options
            restart      = self._restart  ,\
            restart_from = restartFromFunc,\
            # Set additional options
            additional = (
                ("tag",theRunName,0),\
                ("ownFilters"  , "type", ownFilterType),\
                ("switch"      , "useHyperViscAzVortD", useHyperViscAzVortD),\
                         ),\
            series_add = self._series_add                ,\
            # PBS options
            BOUT_run_name         = BOUT_run_name        ,\
            post_process_run_name = post_process_run_name,\
            **self._expandPBSOptions                     ,\
            # Common options
            **self._commonRunnerOptions                  ,\
                        )

        self._expand_dmp_folders, self._expand_PBS_ids =\
        expandRunner.execute_runs(\
            post_processing_function = curPostProcessor,\
            # Declare dependencies
            job_dependencies = self._init_PBS_ids,\
            # This function will be called every time after
            # performing a run
            post_process_after_every_run = True,\
            # Below are the kwargs arguments being passed to
            # the post processing function
            theRunName        = theRunName        ,\
            **self._expandPostOptions             ,\
            # Below are the kwargs given to the
            # restartFromFunc
            aScanPath      = self._initAScanPath  ,\
            scanParameters = self._scanParameters ,\
                                      )
        #}}}
    #}}}

    #{{{_callLinearRunner
    def _callLinearRunner(self):
        """Executes the linear run and the post processing"""
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
            add_noise     = self._boutRunnersNoise
        if self._linearPostOptionsUsesMaxGrad:
            # As this is scan dependent, the driver finds the correct folder
            maxGradRhoFolder = self._expand_dmp_folders[0]
            self._linearPostOptions["maxGradRhoFolder"] =\
                                            maxGradRhoFolder
        # From previous outputs
        self._expandAScanPath = self._expand_dmp_folders[0]
        # Name
        theRunName = self._theRunName + "-2-linearPhase1"
        # PBS options
        BOUT_run_name         = theRunName
        post_process_run_name = "post" + theRunName.capitalize()
        #}}}
        #{{{Run and post processing
        self._linearRun = PBS_runner(\
            # Linear options
            **self._linearOptions         ,\
            # Set restart options
            restart      = self._restart  ,\
            restart_from = restartFromFunc,\
            # Set additional options
            additional = (
                ("tag"   , theRunName           ,0),\
                ("switch", "includeNoise"       , includeNoise ),\
                ("switch", "forceAddNoise"      ,forceAddNoise),\
                ("switch", "useHyperViscAzVortD",useHyperViscAzVortD),\
                ("switch", "saveTerms"          ,saveTerms),\
                         ),\
            series_add = self._series_add                ,\
            # Set eventual noise
            add_noise             = add_noise            ,\
            # PBS options
            BOUT_run_name         = BOUT_run_name        ,\
            post_process_run_name = post_process_run_name,\
            # Common options
            **self._linearPBSOptions                     ,\
            **self._commonRunnerOptions                  ,\
                    )

        self._linear_dmp_folders, self._linear_PBS_ids = \
        self._linearRun.execute_runs(\
            post_processing_function = curPostProcessor,\
            # Declare dependencies
            job_dependencies = self._expand_PBS_ids,\
            # This function will be called every time after
            # performing a run
            post_process_after_every_run = True,\
            # Below are the kwargs arguments being passed to
            # the post processing function
            theRunName        = theRunName        ,\
            **self._linearPostOptions             ,\
            # Below are the kwargs given to the
            # restartFromFunc
            aScanPath      = self._expandAScanPath,\
            scanParameters = self._scanParameters ,\
                )

        if self._boutRunnersNoise:
            # Add file which states what noise is used
            # If ownNoise is used, params are written to the BOUT.log file
            for dmp in self._linear_dmp_folders:
                # Create the file
                with open(os.path.join(dmp, "addnoise.log"), "w") as f:
                    f.write("{}".format(self._boutRunnersNoise))
        #}}}
    #}}}

    #{{{_callTurboRunner
    def _callTurboRunner(self):
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
        post_process_run_name = "post" + theRunName.capitalize()
        # Set aScanPath
        self._turboAScanPath  = self._linear_dmp_folders[0]
        #}}}
        #{{{Run and post processing
        self._turboRun = PBS_runner(\
            # Set turbulence options
            **self._turbulenceOptions       ,\
            # Set restart options
            restart      = self._restart    ,\
            restart_from = restartFromFunc  ,\
            additional = (
                ("tag",theRunName,0),\
                ("switch"      , "useHyperViscAzVortD",useHyperViscAzVortD),\
                ("switch"      , "saveTerms"          ,saveTerms),\
                         ),\
            series_add = self._series_add                ,\
            # PBS options
            BOUT_run_name         = BOUT_run_name        ,\
            post_process_run_name = post_process_run_name,\
            **self._turbulencePBSOptions                 ,\
            # Common options
            **self._commonRunnerOptions                  ,\
                        )

        self._turbo_dmp_folders, self._turbo_PBS_ids =\
        self._turboRun.execute_runs(\
            post_processing_function = curPostProcessor,\
            # Declare dependencies
            job_dependencies = self._linear_PBS_ids,\
            # This function will be called every time after
            # performing a run
            post_process_after_every_run = True,\
            # Below are the kwargs arguments being passed to
            # the post processing function
            theRunName        = theRunName        ,\
            **self._turbulencePostOptions         ,\
            # Below are the kwargs given to the
            # restartFromFunc
            aScanPath      = self._initAScanPath  ,\
            scanParameters = self._scanParameters ,\
                                        )
        #}}}
     #}}}

    #{{{_callPostProcessGrowthRates
    def _callPostProcessGrowthRates(self):
        """Calls the post processor for the growth rates"""

        scanParam  = self._scanParameters[0]
        theRunName = self._theRunName + "-growthRates"
        curPostProcessor = postBoutRunner

        # Make a tuple of tuples, where each subtuple will be used as the
        # paths in collectiveCollect
        collectionFolders =\
            tuple(zip(self._linear_dmp_folders, self._turbo_dmp_folders))

        _, _ = self._linearRun.execute_runs(\
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
            scanParam        = scanParam               ,\
            steadyStatePaths = self._expand_dmp_folders,\
            **self._probesPlotterOptions               ,\
            # Below are the kwargs given to the
            # restartFromFunc
            aScanPath      = self._turboAScanPath ,\
            scanParameters = self._scanParameters ,\
                                    )
    #}}}

    #{{{_callPostProcessProbesAndEnergy
    def _callPostProcessProbesAndEnergy(self):
        """Calls the post processor for the growth rates"""

        collectionFolders = (self._linear_dmp_folders[0],\
                             self._turbo_dmp_folders[0])

        theRunName = self._theRunName + "-energyProbesPlot"
        curPostProcessor = postBoutRunner

        _, _ = self._turboRun.execute_runs(\
            post_processing_function = curPostProcessor,\
            # Declare dependencies
            job_dependencies = self._turbo_PBS_ids,\
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
            **self._probesPlotterOptions               ,\
            tIndSaturatedTurb = self.tIndSaturatedTurb ,\
            # The steady state path will be
            # converted using convertToCurrentScanParameters
            steadyStatePath = self._expand_dmp_folders[0],\
            # Below are the kwargs given to the
            # restartFromFunc and convertToCurrentScanParameters
            aScanPath      = self._tuboAScanPath         ,\
            scanParameters = self._scanParameters)
    #}}}

    #{{{_callpostProcessLinProfiles
    def _callpostProcessLinProfiles(self):
        """Calls the postProcessor for linear profiles"""

        curPostProcessor = postBoutRunner
        theRunName = self._theRunName + "-2-linearPhaseParProfiles"
        tSlice = None

        # Add the tag and the run name
        self._profileRunOptions["additional"] =\
                list(self._profileRunOptions["additional"])
        self._profileRunOptions["additional"].append(("tag",theRunName,0))
        self._profileRunOptions["additional"] =\
                tuple(self._profileRunOptions["additional"])
        self._profileRunOptions["BOUT_run_name"] = theRunName

        # Create the runner
        profileRun = PBS_runner(**self._profileRunOptions)
        # Execute
        _, _ = profileRun.execute_runs(\
            post_processing_function = curPostProcessor,\
            # Declare dependencies
            job_dependencies = self._linear_PBS_ids,\
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
            aScanPath      = self._linearAScanPath,\
            scanParameters = self._scanParameters ,\
            # Common kwargs
            **self._fieldPlotterOptions)
    #}}}

    #{{{_callPostProcessTurbProfiles
    def _callPostProcessTurbProfiles(self):
        curPostProcessor = postBoutRunner
        theRunName = self._theRunName + "-3-turbulentPhase1ParProfiles"
        tSlice = None

        # Add the tag and the run name
        if self.postProcessLinProfiles:
            # Tag is already present in the dict:
            self._profileRunOptions["additional"] =\
                    list(self._profileRunOptions["additional"])
            _ = self._profileRunOptions["additional"].pop()
            self._profileRunOptions["additional"] =\
                    tuple(self._profileRunOptions["additional"])

        self._profileRunOptions["additional"] =\
             list(self._profileRunOptions["additional"])
        self._profileRunOptions["additional"].append(("tag",theRunName,0))
        self._profileRunOptions["additional"] =\
             tuple(self._profileRunOptions["additional"])
        self._profileRunOptions["BOUT_run_name"] = theRunName
        # Create the runner
        profileRun = PBS_runner(**self._profileRunOptions)
        # Execute
        _, _ = profileRun.execute_runs(\
            post_processing_function = curPostProcessor,\
            # Declare dependencies
            job_dependencies = self._turbo_PBS_ids,\
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
            aScanPath      = self._linearAScanPath,\
            scanParameters = self._scanParameters ,\
            # Common kwargs
            **self._fieldPlotterOptions           ,\
                                    )
    #}}}
#}}}

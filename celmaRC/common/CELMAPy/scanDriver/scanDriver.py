#!/usr/bin/env python

"""
Contains the restartFromFunc and ScanDriver
"""

from bout_runners import basic_runner, PBS_runner
import inspect
import os
import pickle
import re

# NOTE: Smells of code duplication in the "call" functions

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

#{{{ScanDriver
class ScanDriver(object):
    #{{{docstring
    """
    Generic scan driver class.

    Constructor sets default options, which can be changed through the
    class' setters.
    """
    #}}}

    #{{{Constructor
    def __init__(self, directory, runner = PBS_runner):
        #{{{docstring
        """
        Sets the default options and selects the runner.

        Parameters
        ----------
        directory : str
            Path to BOUT.inp
        runner : [basic_runner | PBS_runner]
            The runner to use
        """
        #}}}
        # Set warning flags
        self._calledFunctions = {
                    "mainOptions"         : False,\
                    "runTypeOptions"      : False,\
                    "commonRunnerOptions" : False,\
                    "initOptions"         : False,\
                    "expandOptions"       : False,\
                    "linearOptions"       : False,\
                    "turbulenceOptions"   : False,\
                }

        self._directory = directory
        self._runner    = runner
        #}}}

    #{{{setMainOptions
    def setMainOptions(self                         ,\
                       scanParameters               ,\
                       series_add                   ,\
                       theRunName                   ,\
                       make                  = False,\
                       boutRunnersNoise      = None ,\
                       restartFrom           = None ,\
                      ):
        #{{{docstring
        """
        Sets the main options for the scan.

        Parameters
        ----------
        scanParameters : sequence of strings
            Sequence of all quantities that will change during a scan
        series_add : sequence of sequence which is not string
            The series_add for the scan
        theRunName : str
            Name of the run
        make : bool
            Shall the program be made or not
        boutRunnersNoise : [None|dict]
            If this is None, the noise will be generated
            from the noise generator in noiseGenerator.cxx
            If this is set as a dict, the dict must have the key of one
            of the variables evolved in time, and the value must be the
            amplitude of the noise.
        restartFrom : [None|str]
            Set this to a path if you would like to restart the run from
            one particular file.
        """
        #}}}

        if boutRunnersNoise is None:
            boutRunnersNoise = False

        self._calledFunctions["mainOptions"] = True

        self._scanParameters        = scanParameters
        self._series_add            = series_add
        self._theRunName            = theRunName
        self._make                  = make
        self._boutRunnersNoise      = boutRunnersNoise
        self._restartFrom           = restartFrom
    #}}}

    #{{{setRunTypeOptions
    def setRunTypeOptions(self            ,\
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
            self._calledFunctions["runTypeOptions"] = True

        self._runOptions =\
                {\
                  "runInit"   : runInit  ,\
                  "runExpand" : runExpand,\
                  "runLin"    : runLin   ,\
                  "runTurb"   : runTurb  ,\
                }
    #}}}

    #{{{setCommonRunnerOptions
    def setCommonRunnerOptions(self               ,\
                               nproc        = 48  ,\
                               cpy_source   = True,\
                               BOUT_nodes   = 3   ,\
                               BOUT_ppn     = 16  ,\
                               BOUT_queue   = None,\
                               BOUT_account = None,\
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
            How many nodes to run on (only used when the
            runner is PBS_runner)
        BOUT_ppn : int
            How many processors per node (only used when the
            runner is PBS_runner)
        BOUT_queue : [None|str]
            What queue to submit the jobs to (only used when the
            runner is PBS_runner)
        BOUT_account : [None|str]
            Account number to use for the runs (only used when the
            runner is PBS_runner)
        """
        #}}}

        self._calledFunctions["commonRunnerOptions"] = True

        self._commonRunnerOptions =\
                {\
                 "nproc"               : nproc       ,\
                 "cpy_source"          : cpy_source  ,\
                }

        if self._runner == PBS_runner:
            self._commonRunnerOptions["BOUT_nodes"]   = BOUT_nodes
            self._commonRunnerOptions["BOUT_ppn"  ]   = BOUT_ppn
            self._commonRunnerOptions["BOUT_queue"]   = BOUT_queue
            self._commonRunnerOptions["BOUT_account"] = BOUT_account
    #}}}

    #{{{setInitOptions
    def setInitOptions(self                              ,\
                       timestep              = 2e3       ,\
                       nout                  = 2         ,\
                       BOUT_walltime         = "24:00:00"):
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
                 "timestep" : (timestep,),\
                 "nout"     : (nout,),\
                }

        self._initPBSOptions =\
                {\
                 "BOUT_walltime" : BOUT_walltime,\
                }

    #}}}

    #{{{setExpandOptions
    def setExpandOptions(self                              ,\
                         nz                    = 256       ,\
                         timestep              = 50        ,\
                         nout                  = 2         ,\
                         BOUT_walltime         = "24:00:00",\
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
                 "timestep" : (timestep,),\
                 "nout"     : (nout,),\
                }

        self._expandPBSOptions =\
                {\
                 "BOUT_walltime" : BOUT_walltime,\
                }

    #}}}

    #{{{setLinearOptions
    def setLinearOptions(self                              ,\
                         timestep              = 1         ,\
                         nout                  = 1000      ,\
                         BOUT_walltime         = "72:00:00",\
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
                 "timestep" : (timestep,),\
                 "nout"     : (nout,),\
                }

        self._linearPBSOptions =\
                {\
                 "BOUT_walltime" : BOUT_walltime ,\
                }
    #}}}

    #{{{setTurbulenceOptions
    def setTurbulenceOptions(self                      ,\
                             timestep      = 1         ,\
                             nout          = 5000      ,\
                             BOUT_walltime = "72:00:00",\
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
                 "timestep" : (timestep,),\
                 "nout"     : (nout,),\
                }

        self._turbulencePBSOptions =\
                {\
                 "BOUT_walltime" : BOUT_walltime,\
                }
    #}}}

    #{{{runScan
    def runScan(self, boussinesq=False, restartTurb=None):
        """
        Calls the drivers used for running scans

        Parameters
        ----------
        boussinesq : bool
            Whether or not boussinesq approximation is used
        restartTurb : [None|int]
            Number of times to restart the trubulence runs
        """

        # Set boussinesq and restart turb
        self._boussinesq  = boussinesq
        self._restartTurb = restartTurb

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

        if self._runner == basic_runner:
            # Remove PBS option
            members = tuple(member for member in dir(self) if "PBS" in member)
            for member in members:
                setattr(self, member, {})

        # Make dictionary to variables
        for (flag, value) in self._runOptions.items():
            setattr(self, flag, value)

        # Update dicts
        self._commonRunnerOptions["directory"] = self._directory

        # Check that no troubles running with restart from
        if self._restartFrom is not None:
            messageTemplate = "Cannot run {} when restart from is set\n"
            tests = {"runInit" : self.runInit}
            problem = False
            message = ""
            for key in tests.keys():
                if tests[key]:
                    problem = True
                    message += messageTemplate.format(key)

            if problem:
                raise ValueError(message)

        # Set the path to the dmp folders
        self._dmpFoldersDictPath =\
                os.path.join(self._directory, "dmpFoldersDict.pickle")

        # Call the runners
        if self.runInit:
            self._callInitRunner()
            # Load the dmpFolders pickle, update it and save it
            dmpFoldersDict = self._getDmpFolderDict()
            dmpFoldersDict["init"] = self._init_dmp_folders
            self._pickleDmpFoldersDict(dmpFoldersDict)

        self._restart = "overwrite"

        if self.runExpand:
            self._callExpandRunner()
            # Load the dmpFolders pickle, update it and save it
            dmpFoldersDict = self._getDmpFolderDict()
            dmpFoldersDict["expand"] = self._expand_dmp_folders
            self._pickleDmpFoldersDict(dmpFoldersDict)

        if self.runLin:
            self._callLinearRunner()
            # Load the dmpFolders pickle, update it and save it
            dmpFoldersDict = self._getDmpFolderDict()
            dmpFoldersDict["linear"] = self._linear_dmp_folders
            self._pickleDmpFoldersDict(dmpFoldersDict)

        if self.runTurb:
            self._callTurboRunner()
            # Load the dmpFolders pickle, update it and save it
            dmpFoldersDict = self._getDmpFolderDict()
            dmpFoldersDict["turbulence"] = self._turbo_dmp_folders
            self._pickleDmpFoldersDict(dmpFoldersDict)


        if self._restartTurb is not None:
            self._callExtraTurboRunner()
            # Load the dmpFolders pickle, update it and save it
            dmpFoldersDict = self._getDmpFolderDict()
            dmpFoldersDict["extraTurbulence"] = list(dmpFoldersDict["extraTurbulence"])
            dmpFoldersDict["extraTurbulence"].append(self._extra_turbo_folders)
            dmpFoldersDict["extraTurbulence"] = tuple(dmpFoldersDict["extraTurbulence"])
            self._pickleDmpFoldersDict(dmpFoldersDict)
    #}}}

    #{{{_getDmpFolderDict
    def _getDmpFolderDict(self):
        """
        Returns the dmpFolderDict
        """

        if os.path.exists(self._dmpFoldersDictPath):
            with open(self._dmpFoldersDictPath, "rb") as f:
                dmpFoldersDict = pickle.load(f)
        else:
            dmpFoldersDict = {\
                              "init"            : None,\
                              "expand"          : None,\
                              "linear"          : None,\
                              "turbulence"      : None,\
                              "extraTurbulence" : ()  ,\
                             }

        return dmpFoldersDict
    #}}}

    #{{{_pickleDmpFoldersDict
    def _pickleDmpFoldersDict(self, dmpFoldersDict):
        """
        Pickles the dmpFolderDict
        """

        with open(self._dmpFoldersDictPath, "wb") as f:
            pickle.dump(dmpFoldersDict, f)
    #}}}

    #{{{_callInitRunner
    def _callInitRunner(self):
        """Executes the init run"""

        #{{{Init options
        # Name
        theRunName = self._theRunName + "-0-initialize"
        # Set the spatial domain
        nz = 1
        # Set the temporal domain
        self._restart = None
        # Filter
        ownFilterType = "none"
        if not(self._boussinesq):
            hyper = ("switch", "useHyperViscAzVortD",False)
        else:
            hyper = ("switch", "useHyperViscAzVort",False)
        # PBS options
        if self._runner == PBS_runner:
            self._initPBSOptions.update({"BOUT_run_name": theRunName })
        #}}}
        #{{{Run
        initRunner = self._runner(\
            make       = self._make   ,\
            nz         = nz           ,\
            restart    = self._restart,\
            **self._initOptions       ,\
            # Set additional option
            additional = (
                ("tag",        theRunName, 0            ),\
                ("ownFilters", "type"    , ownFilterType),\
                hyper                                    ,\
                         )               ,\
            series_add = self._series_add,\
            # PBS options
            **self._initPBSOptions,\
            # Common options
            **self._commonRunnerOptions,\
                        )

        self._init_dmp_folders, self._init_PBS_ids =\
            initRunner.execute_runs()
        #}}}
    #}}}

    #{{{_callExpandRunner
    def _callExpandRunner(self):
        """Executes the expand run"""

        #{{{Expand runner options
        # Set the spatial domain
        # Filter
        ownFilterType = "none"
        if not(self._boussinesq):
            hyper = ("switch", "useHyperViscAzVortD",False)
        else:
            hyper = ("switch", "useHyperViscAzVort",False)
        # From previous outputs
        try:
            self._initAScanPath = self._init_dmp_folders[0]
        except AttributeError as er:
            if "has no attribute" in er.args[0]:
                raise RuntimeError("Init must be run prior to expand")
            else:
                raise er

        # Name
        theRunName = self._theRunName + "-1-expand"
        # PBS options
        if self._runner == PBS_runner:
            self._expandPBSOptions.update({"BOUT_run_name" : theRunName})

        restart_from =\
                self._restartFrom if self._restartFrom else restartFromFunc
        #}}}
        #{{{Run
        expandRunner = self._runner(\
            restart      = self._restart,\
            restart_from = restart_from ,\
            **self._expandOptions       ,\
            # Set additional options
            additional = (
                ("tag"       ,theRunName, 0            ),\
                ("ownFilters", "type"   , ownFilterType),\
                hyper,\
                         )               ,\
            series_add = self._series_add,\
            # PBS options
            **self._expandPBSOptions   ,\
            # Common options
            **self._commonRunnerOptions,\
                        )

        self._expand_dmp_folders, self._expand_PBS_ids =\
        expandRunner.execute_runs(\
            # Declare dependencies
            job_dependencies = self._init_PBS_ids,\
            # Below are the kwargs given to the
            # restartFromFunc
            aScanPath      = self._initAScanPath  ,\
            scanParameters = self._scanParameters ,\
                                      )
        #}}}
    #}}}

    #{{{_callLinearRunner
    def _callLinearRunner(self):
        """Executes the linear run"""

        #{{{ Linear options
        #Switches
        saveTerms = False
        if not(self._boussinesq):
            hyper = ("switch", "useHyperViscAzVortD",True)
        else:
            hyper = ("switch", "useHyperViscAzVort",True)
        if not(self._boutRunnersNoise):
            includeNoise  = True
            forceAddNoise = True
            add_noise     = None
        else:
            includeNoise  = False
            forceAddNoise = False
            add_noise     = self._boutRunnersNoise
        try:
            # From previous outputs
            self._expandAScanPath = self._expand_dmp_folders[0]
        except AttributeError as er:
            if "has no attribute" in er.args[0]:
                raise RuntimeError("Expand must be run prior to linear")
            else:
                raise er

        # Name
        theRunName = self._theRunName + "-2-linearPhase1"
        # PBS options
        if self._runner == PBS_runner:
            self._linearPBSOptions.update({"BOUT_run_name" : theRunName})

        restart_from =\
                self._restartFrom if self._restartFrom else restartFromFunc
        #}}}
        #{{{Run
        self._linearRun = self._runner(\
            restart      = self._restart,\
            restart_from = restart_from ,\
            **self._linearOptions       ,\
            # Set additional options
            additional = (
                ("tag"   , theRunName      ,0),\
                ("switch", "includeNoise"  , includeNoise ),\
                ("switch", "forceAddNoise" ,forceAddNoise) ,\
                ("switch", "saveTerms"     ,saveTerms)     ,\
                hyper,\
                         )               ,\
            series_add = self._series_add,\
            # Set noise
            add_noise = add_noise,\
            # PBS options
            **self._linearPBSOptions,\
            # Common options
            **self._commonRunnerOptions,\
                    )

        self._linear_dmp_folders, self._linear_PBS_ids = \
        self._linearRun.execute_runs(\
            # Declare dependencies
            job_dependencies = self._expand_PBS_ids,\
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
        """Executes the turbulence run"""

        #{{{Turbulence options
        # Switches
        saveTerms           = False
        if not(self._boussinesq):
            hyper = ("switch", "useHyperViscAzVortD",True)
        else:
            hyper = ("switch", "useHyperViscAzVort",True)

        # Name
        theRunName = self._theRunName + "-3-turbulentPhase1"
        # PBS options
        if self._runner == PBS_runner:
            self._turbulencePBSOptions.update({"BOUT_run_name" : theRunName })

        # Set aScanPath
        try:
            self._linearAScanPath  = self._linear_dmp_folders[0]
        except AttributeError as er:
            raise RuntimeError("Linear must be run prior to turbulence")

        restart_from =\
                self._restartFrom if self._restartFrom else restartFromFunc
        #}}}
        #{{{Run
        self._turboRun = self._runner(\
            # Set turbulence options
            **self._turbulenceOptions   ,\
            # Set restart options
            restart      = self._restart,\
            restart_from = restart_from ,\
            additional = (
                ("tag"   , theRunName , 0        ),\
                ("switch", "saveTerms", saveTerms),\
                hyper,\
                         )               ,\
            series_add = self._series_add,\
            # PBS options
            **self._turbulencePBSOptions,\
            # Common options
            **self._commonRunnerOptions ,\
                        )

        self._turbo_dmp_folders, self._turbo_PBS_ids =\
        self._turboRun.execute_runs(\
            # Declare dependencies
            job_dependencies = self._linear_PBS_ids,\
            # Below are the kwargs given to the
            # restartFromFunc
            aScanPath      = self._linearAScanPath ,\
            scanParameters = self._scanParameters ,\
                                        )
        #}}}
    #}}}

    #{{{_callExtraTurboRunner
    def _callExtraTurboRunner(self):
        """
        Calls the turbulence runner in a loop, so that it can be restarted.
        """

        saveTerms   = False

        # Name
        theRunName = self._theRunName + "-3-turbulentPhase1"

        # Switches
        if not(self._boussinesq):
            hyper = ("switch", "useHyperViscAzVortD",True)
        else:
            hyper = ("switch", "useHyperViscAzVort",True)
        # Make extra turbulence folder
        self._extra_turbo_folders = []
        # Create new runner
        self._turboRun = self._runner(\
            **self._turbulenceOptions,\
            restart    = "overwrite" ,\
            additional = (
                ("tag"   , theRunName , 0        ),\
                ("switch", "saveTerms", saveTerms),\
                hyper,\
                         )               ,\
            series_add = self._series_add,\
            # PBS options
            **self._turbulencePBSOptions ,\
            # Common options
            **self._commonRunnerOptions  ,\
                        )

        for nr in range(self._restartTurb):
            turbo_dmp_folders, self._turbo_PBS_ids =\
                self._turboRun.execute_runs(\
                    # Declare dependencies
                    job_dependencies = self._turbo_PBS_ids,\
                    # Below are the kwargs given to the
                    # restartFromFunc
                    aScanPath      = self._linearAScanPath,\
                    scanParameters = self._scanParameters ,\
                                                )
            # Update the list
            self._extra_turbo_folders.append(turbo_dmp_folders)

        # Recast to tuple
        self._extra_turbo_folders = tuple(self._extra_turbo_folders)
    #}}}
#}}}

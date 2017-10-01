#!/usr/bin/env python

"""
Contains the restartFromFunc and ScanDriver
"""

from bout_runners import basic_runner, PBS_runner
import inspect, os, pickle, re, pathlib, shutil, difflib

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

        self._scanParameters   = scanParameters
        self._series_add       = series_add
        self._theRunName       = theRunName
        self._make             = make
        self._boutRunnersNoise = boutRunnersNoise
        self._restartFrom      = restartFrom
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
                               BOUT_mail    = None,\
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
        BOUT_mail : [None|str]
            Mail address to notify when a BOUT job has finished
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
            self._commonRunnerOptions["BOUT_mail"]    = BOUT_mail
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
    def runScan(self,\
                boussinesq=False,\
                restartTurb=None,\
                checkForEmptyRestarts=True):
        """
        Calls the drivers used for running scans

        Parameters
        ----------
        boussinesq : bool
            Whether or not boussinesq approximation is used
        restartTurb : [None|int]
            Number of times to restart the trubulence runs
        checkForEmptyRestarts : bool
            If True, the fuction will search for folders which contains
            *.restart.*-files without *.dmp.*-files (with exception of
            the restart_0 folder). The function will start running the
            simulation of these folders without any queue dependecy.
        """

        # Set boussinesq and restart turb
        self._boussinesq  = boussinesq
        self._restartTurb = restartTurb

        if not(self._calledFunctions["mainOptions"]):
            message = "self.setMainOptions must be called prior to a run"
            raise ValueError(message)

        # Set default options for unset values
        self._setNotCalledToDefault()

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

        # Find the empty restarts, set to "set()" if switched off
        if checkForEmptyRestarts:
            self._emptyRestarts = self._searchForEmptyRestarts()
        else:
            self._emptyRestarts = set()

        # YOU ARE HERE:
        # TESTS: - One emptyRestart in expand
        #        - Two in lin
        #        - One expand, one lin
        #        - Just turb
        #        - The whole range [with one restart] (one)
        #        - One of the restarts missing (with 2 restarts)
        if len(self._emptyRestarts) == 0:
            self._normalRun()
        else:
            # Find the project root
            self._projectRoot = list(list(self._emptyRestarts)[0].parents)[-2]
            self._runOnlyFromRestart()
    #}}}

    #{{{_setNotCalledToDefault
    def _setNotCalledToDefault(self):
        """
        Sets non-called options to default values.

        Will also print the status of the options.
        """

        keysToBeCalled = []

        for key, val in self._calledFunctions.items():
            if val is not(None):
                if val:
                    print("{:<25} has been called".format(key))
                else:
                    keysToBeCalled.append(key)
                    print("{:<25} has NOT been called (using defaults)".\
                          format(key))

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
    #}}}

    #{{{_normalRun
    def _normalRun(self):
        """
        Normal run procedure (as opposed to _runOnlyFromRestart)

        * The simulations are performed in chronological order
        * The simulations are performed with dependencies
        * The dmpFoldersDict is updated normally
        """

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
# FIXME: END


# Idea: Could move restart_0 to rst_from_linear and restart files in
#       root to restart_last_restart_files
#       In this way: Could start running in the restart folder without
#       hasle
# FIXME: Need the self._turbo_dmp_folders in order do run this
        if self._restartTurb is not None:
            # NOTE: dmpFolders are treated internally in this function
            self._callExtraTurboRunner()
    #}}}

    #{{{_searchForEmptyRestarts
    def _searchForEmptyRestarts(self):
        """
        Searches for folder containing *.restart.* files without *.dmp.* files.

        Also moves restart_0 folders to rst_BAK_*

        Returns
        -------
        onlyRestart : set
            A set of all folders containing *.restart.*, but no *.dmp.* files.
        """

        restartFiles = list(pathlib.Path(self._directory).glob("**/*.restart.*"))
        dmpFiles = list(pathlib.Path(self._directory).glob("**/*.dmp.*"))

        # .parents[0] is the directory holding the file
        # https://stackoverflow.com/questions/35490148/how-to-get-folder-name-in-which-given-file-resides-from-pathlib-path
        restartFolders = set(f.parents[0] for f in restartFiles)
        dmpFolders = set(f.parents[0] for f in dmpFiles)

        # Non-symmetric difference
        # https://stackoverflow.com/questions/3462143/get-difference-between-two-lists
        onlyRestart = restartFolders - dmpFolders

        # Extraxt restart_0, rst_BAK_* and root_rst_files folders
        restart0Folders = set(e for e in onlyRestart if "restart_0" in str(e))
        rstBAKFolders = set(e for e in onlyRestart if "rst_BAK" in str(e))
        rootRstFiles = set(e for e in onlyRestart if "root_rst_files" in str(e))
        onlyRestart -= restart0Folders
        onlyRestart -= rstBAKFolders
        onlyRestart -= rootRstFiles

        # Move restart_0 folders to rst_BAK_*
        for f in restart0Folders:
            # Check the number of rst_BAK_* folders already present
            curRstBak = list(pathlib.Path(f.parents[0]).glob("rst_BAK*"))
            curRstBak = [e for e in curRstBak if e.is_dir()]
            newRstBakFolder =\
                f.parents[0].joinpath("rst_BAK_{}".format(len(curRstBak)))
            os.makedirs(str(newRstBakFolder))
            f.rename(newRstBakFolder)

        return onlyRestart
    #}}}

    #{{{_runOnlyFromRestart
    def _runOnlyFromRestart(self):
        """
        Run procedure which runs from 'empty restart' (as opposed to _normalRun)

        * The simulations are performed in reverse order
        * The simulations are performed without dependencies
        * The dmpFoldersDict is not updated
        """

        print(("Found the following folders with *.restart.* files"
               " without *.dmp.* files:"))
        for f in self._emptyRestarts:
            print("   * {}".format(f))
        print("\n'runScan' called with 'checkForEmptyRestarts', rerun "\
              "the driver in order to fix 'dmpFoldersDict.pickle'")
        if self.runTurb:
            # Move the root restart files to root_rst_files
            turboRoots = [e for e in self._emptyRestarts if \
                         "turbulentPhase1" in str(e)]
            if len(turboRoots) != 0:
                self._linear_dmp_folders = self._findPreviousFolders("linear")
                self._moveRootRestart(turboRoot)
                self._linear_PBS_ids = None
                # Run the simulation
                self._callTurboRunner()
        if self.runLin:
            linearRoots = [e for e in self._emptyRestarts if \
                          "linearPhase1" in str(e)]
            if len(linearRoots) != 0:
                self._expand_dmp_folders = self._findPreviousFolders("expand")
                self._moveRootRestart(linearRoots)
                self._expand_PBS_ids = None
                self._callLinearRunner()
        if self.runExpand:
            expandRoots = [e for e in self._emptyRestarts if \
                          "expand" in str(e)]
            if len(expandRoots) != 0:
                # self._init_dmp_folders found in runScan
                self._moveRootRestart(expandRoots)
                self._init_PBS_ids = None
                self._callExpandRunner()

# FIXME: Need own logic here
# FIXME: Question? Would need own logic on init run, or all subsequent
#        runs?
# FIXME: Fails as no self._turbo_dmp_folders
        import pdb; pdb.set_trace()
        if self._restartTurb is not None:
            # NOTE: dmpFolders are treated internally in this function
            self._callExtraTurboRunner()
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

    #{{{_moveRootRestart
    def _moveRootRestart(self, rootFolders):
        """
        Moves the restart files in a root folder to root_rst_files

        Parameters
        ----------
        rootFolders : list
            List containing the folders to move the *.restart.* files from
        """

        for root in rootFolders:
            # Find all files
            content = root.glob("*")
            files = tuple(e for e in content if e.is_file())

            # Create new folder
            rootRstFolder = root.joinpath("root_rst_files")
            os.makedirs(str(rootRstFolder))

            for f in files:
                shutil.move(str(f), str(rootRstFolder))
    #}}}

    #{{{_pickleDmpFoldersDict
    def _pickleDmpFoldersDict(self, dmpFoldersDict):
        """
        Pickles the dmpFolderDict
        """

        with open(self._dmpFoldersDictPath, "wb") as f:
            pickle.dump(dmpFoldersDict, f)
    #}}}

    #{{{_findPreviousFolders
    def _findPreviousFolders(self, folderName):
        """
        Finds the dmp folders from where we are restarting from

        Parameters
        ----------
        folderName : str
            Folder name to search for, like for example 'expand'

        Returns
        -------
        previousFolders : tuple
            Tuple of the previous dmp folders
        """
        # Search for all different "expand" folders
        expandFolders = tuple(projectRoot.glob("**/*expand*/*"))
        expandFolders = tuple(set(e.parents[0] for e in expandFolders))
        # Ensure that only scanParameter is changing
        for ef in expandFolders:
            self._checkOnlyScanParametersVaries(str(expandFolders[0]), str(ef))

        previousFolders = tuple(str(e) for e in expandFolders)

        return previousFolders
    #}}}

    #{{{_checkOnlyScanParametersVaries
    def _checkOnlyScanParametersVaries(self, a, b):
        """
        Asserts that only the scan parameter varies in a folder.

        Throws an error if other elements, such as switches varies.

        Parameters
        ----------
        a : str
            One of the paths to be compared
        b : str
            The other path to be compared
        """

        # Find the difference between two strings
        # a = "foo_bar" b = "foo_baz" => ['  f  o  o  ', '  b  a- r+ z']
        # https://stackoverflow.com/questions/17904097/python-difference-between-two-strings
        splittedText = "".join(difflib.ndiff(a, b)).split("_")
        # Make a diff-like split of the elements (nested list comprehension)
        # => (('  f', '  o', '  o', '  '), ('  b', '  a', '- r', '+ z'))
        # https://stackoverflow.com/questions/9475241/split-string-every-nth-character
        chunkedComparison = tuple(tuple(e[i:i+3] for i in range(0, len(e), 3))\
                                  for e in splittedText)

        onlyDiffs = []
        for e in chunkedComparison:
            # If there is no difference, combine the subelement
            # => "foo"
            # else, keep to subelement
            # ('  b', '  a', '- r', '+ z')
            if all(char[0] == " " for char in e):
                onlyDiffs.append("".join(e).replace(" ", ""))
            else:
                onlyDiffs.append(e)

        # If there is a diff, and the previous element was not in scanParameters
        # success will be false
        success = True
        for i in range(len(onlyDiffs)):
        if type(onlyDiffs[i]) is not str:
          if i == 0:
              success = False
          if onlyDiffs[i-1] not in self._scanParameters:
              success = False

        if not success:
            message =\
                ("Could not run in empty restart mode as folders belonging to"
                "different scans were found in the same folder, for example"
                "\n{}and\n{}").format(a, b)
            raise RuntimeError(message)
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

        add_noise = self._boutRunnersNoise

        try:
            # From previous outputs
            import pdb; pdb.set_trace()
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
                ("tag"   , theRunName , 0),\
                ("switch", "saveTerms", saveTerms),\
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
        # Create new runner
        self._turboRun = self._runner(\
            **self._turbulenceOptions,\
            restart    = "overwrite" ,\
            additional = (
# FIXME: tag
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

        # Check if runs have already been performed
        restartNumber = self._getNumberOfPerformedRestarts()

        # Update the restart number
        additionalRestartsNeeded = self._restartTurb - restartNumber

        for nr in range(additionalRestartsNeeded):
# FIXME: Actually: Maybe the best would be to crawl through the folders,
# and check for folders with only restart (no dmp), can read from those
# if linear, turb or similar
# NOT restart_0 folders!
# !!! In init, expand, linear and turb: restart_0 NOT superfluous!
# !!! In init, expand and linear: restart files in root ARE superfluous? Well, files are first copied to restart_0
# ? Will restart start if restart_0 is populated? If so, no extra logic needed (at least not in expand and linear)...but problem if starts the preceding first, and make dependencies from that
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! If root restart files does not exist: will copy from preceding root folder => must run the restart files in reverse!!!
# Latest turbulence is the latest, usually: Will not need to restart this, as this will append the total amount of runs
# TODO: Up an including linear: Only extra logic needed: reverse order


# NOTE: Hypothesis: numberOfRestartFolders could be 0 or neg, in that case, no new run
#       would be performed, must be checked

# NOTE: nr will start on 0, will refer to restart_0? check
# NOTE: Potential problem: dmp_folders are given after the run has been
#       given...however, do not believe this would be a problem here, as
#       turbo_dmp_folders are input
# FIXME: appropriate to add check-logic here?
#        NOTE: calling bout_runners, but created above
            turbo_dmp_folders, self._turbo_PBS_ids =\
                self._turboRun.execute_runs(\
                    # Declare dependencies
                    job_dependencies = self._turbo_PBS_ids,\
                    # Below are the kwargs given to the
                    # restartFromFunc through bout_runners
                    aScanPath      = self._linearAScanPath,\
                    scanParameters = self._scanParameters ,\
                                            )

            # Load the dmpFolders pickle, update it and save it
            dmpFoldersDict = self._getDmpFolderDict()
            dmpFoldersDict["extraTurbulence"] =\
                list(dmpFoldersDict["extraTurbulence"])
            for dmpNr in range(len(turbo_dmp_folders)):
                dmpFoldersDict["extraTurbulence"][dmpNr] =\
                    list(dmpFoldersDict["extraTurbulence"][dmpNr])
                dmpFoldersDict["extraTurbulence"][dmpNr].\
                    append(turbo_dmp_folders[dmpNr])
                dmpFoldersDict["extraTurbulence"][dmpNr] =\
                    tuple(dmpFoldersDict["extraTurbulence"][dmpNr])

            dmpFoldersDict["extraTurbulence"] =\
                tuple(dmpFoldersDict["extraTurbulence"])
            self._pickleDmpFoldersDict(dmpFoldersDict)
    #}}}

    #{{{_getNumberOfPerformedRestarts
    def _getNumberOfPerformedRestarts(self):
        """
        Will return the number of performed restarts.

        Also checks that the scans has an equal amount of restarts

        Returns
        -------
        restartNumber : int
            The number of already performed restarts.
        """

        numberOfRestartFolders = []
        # Placeholder for the dmp folders
        extraTrubulenceRuns = []
        for turboDmp in self._turbo_dmp_folders:
            # Find the folders
            folders = tuple(folder for folder in os.listdir(turboDmp) if
                            os.path.isdir(os.path.join(turboDmp, folder)))
            # NOTE: bout_runners is making a restart 0 folder for
            #       copied restart files, as a hack, the 0th restart
            #       folder is excluded

            # Remove un-needed folders
            folders = tuple(folder for folder in folders if\
                            ("restart" in folder) and (folder != "restart_0"))

            # Sort
            # http://stackoverflow.com/questions/4836710/does-python-have-a-built-in-function-for-string-natural-sort
            # NOTE: The lambda function returns a tuple for each element
            #       The sorting is made from the first element in the
            #       tuple
            #       If the first element is equal, then the sorting
            #       happens for in the next element, and so on.
            digit   = re.compile(r"(\d+)")
            folders =\
                sorted(folders,\
                    key=lambda el: [int(s) if s.isdigit() else s.lower()\
                        for s in re.split(digit, el)])

            extraTrubulenceRuns.append(\
                   tuple(os.path.join(turboDmp, folder) for folder in folders))

            numberOfRestartFolders.append(len(folders))

        # Load the dmpFolders pickle, update it and save it
        dmpFoldersDict = self._getDmpFolderDict()
        dmpFoldersDict["extraTurbulence"] = tuple(extraTrubulenceRuns)
        self._pickleDmpFoldersDict(dmpFoldersDict)

        if numberOfRestartFolders.count(numberOfRestartFolders[0]) != len(numberOfRestartFolders):
            message = ("The dmp folders did not have an equal number of "
                       "restarts. The restart counter remains ambigous. "
                       "Fix by manually restarting missing restart "
                       "simulations.")
            raise RuntimeError(message)

        restartNumber = numberOfRestartFolders[0]

        return restartNumber
    #}}}
#}}}

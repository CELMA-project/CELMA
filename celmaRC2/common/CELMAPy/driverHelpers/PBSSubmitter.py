#!/usr/bin/env python

"""
Class for submitting functions using PBS
"""

from ..driverHelpers import getTime
from subprocess import run, PIPE
import inspect
import os

#{{{PBSSubmitter
class PBSSubmitter(object):
    """
    Class which can be used to submit functions (from a file) using PBS.

    Member functions of classes is not yet supported.
    Functions using "from .." or "from ." is not yet supported.
    """

    #{{{constructor
    def __init__(self):
        """
        Sets functions not called and call setMisc
        """

        self._notCalled = [\
                            "setJobName" ,\
                            "setNodes"   ,\
                            "setQueue"   ,\
                            "setWalltime",\
                          ]

        self._submitWithPBS = True
    #}}}

    #{{{toggleSubmitOrRun
    def toggleSubmitOrRun(self):
        """
        Toggles submit or run.
        """

        # Use xor
        self._submitWithPBS = bool(self._submitWithPBS^1)

        print("Submission to PBS queue is: {}".format(self._submitWithPBS))
    #}}}

    #{{{setMisc
    def setMisc(self               ,\
                mail    = None     ,\
                account = None     ,\
                logPath = "postLogs"):
        #{{{docstring
        """
        Sets miscellaneous options, creates logPath if not set.

        Parameters
        ----------
        queue : [None|str]
            Name of the queue to run on
        mail : [None|str]
            Mail address to send status to
        account : [None|str]
            Account number to use
        logPath : str
            Path where the submission log and submission errors are
            posted to.
        """
        #}}}

        self._mail    = mail
        self._account = account
        self._logPath = logPath

        # Make dir if not exists
        if not os.path.exists(logPath):
            os.makedirs(logPath)

        self._miscCalled = True
    #}}}

    #{{{setJobName
    def setJobName(self, jobName):
        #{{{docstring
        """
        Sets the job name, and gets the time

        Parameters
        ----------
        jobName : str
            Name of the job
        """
        #}}}

        self._jobName = jobName
        self._time    = getTime(depth = "microsecond")

        if "setJobName" in self._notCalled:
            self._notCalled.remove("setJobName")
    #}}}

    #{{{setNodes
    def setNodes(self, nodes, ppn):
        #{{{docstring
        """
        Sets the job name, and gets the time

        Parameters
        ----------
        nodes : int
            Number of nodes to use
        ppn : int
            Processors per node to use
        """
        #}}}

        self._nodes = nodes
        self._ppn   = ppn

        if "setNodes" in self._notCalled:
            self._notCalled.remove("setNodes")
    #}}}

    #{{{setQueue
    def setQueue(self, queue):
        #{{{docstring
        """
        Sets the queue to use

        Parameters
        ----------
        queue : [None|str]
            Name of the queue to run on
        """
        #}}}

        self._queue = queue

        if "setQueue" in self._notCalled:
            self._notCalled.remove("setQueue")
    #}}}

    #{{{setWalltime
    def setWalltime(self, walltime):
        #{{{docstring
        """
        Sets the walltime

        Parameters
        ----------
        walltime : [None|str]
            Walltime on the format "HH:MM:SS"
        """
        #}}}

        self._walltime = walltime

        if "setWalltime" in self._notCalled:
            self._notCalled.remove("setWalltime")
    #}}}

    #{{{submitFunction
    def submitFunction(self               ,\
                       function           ,\
                       args         = ()  ,\
                       kwargs       = {}  ,\
                       dependencies = None):
        #{{{docstring
        """
        Function which submits a function using PBS

        This is done by making a self deleting temporary python file
        that will be called by a self deleting PBS script.

        Parameters
        ----------
        args : [None|tuple]
            Tuple of the positional arguments to use
        kwargs : [None|dict]
            Dictionary of the keyword arguments to use
        dependencies : [None|tuple]
            The job will be set on hold until the dependencies are
            finished.
        """
        #}}}

        # Guard
        if len(self._notCalled) > 0:
            message = "The following functions were not called:\n{}".\
                        format("\n".join(self._notCalled))
            raise RuntimeError(message)

        if not(self._miscCalled):
            self.setMisc()

        if self._submitWithPBS:
            # The name of the file
            fileName = "tmp_{}_{}.py".format(function.__name__, self._time)

            # Make the script
            functionPath = inspect.getmodule(function).__file__
            with open(functionPath, "r") as f:
                tmpFile = f.read()

            # Add arguments
            tmpFile += "\n#Added by PBSSubmitter\n"
            tmpFile += "args={}\nkwargs={}\n".format(args, kwargs)

            # Add function call in the end of the script
            tmpFile += "{}({},{})\n".\
                    format(function.__name__, "*args", "**kwargs")

            # When the script has run, it will delete itself
            tmpFile += "import os\nos.remove('{}')\n".format(fileName)

            # Write the python script
            with open(fileName, "w") as f:
                f.write(tmpFile)

            # Get core of the job string
            jobString = self._createPBSCoreString()

            # Call the python script in the submission
            jobString += "python {}\n".format(fileName)
            jobString += "exit"

            # Create the dependencies
            if dependencies is not None:
                dependencies = ":".join(dependencies)

            # Submit the job
            print("\nSubmitting '{}'\n".format(self._jobName))

            self._submit(jobString, dependentJob = dependencies)
        else:
            # Running the job
            print("\nRunning '{}'\n".format(self._jobName))
            function(*args, **kwargs)
    #}}}

    #{{{_createPBSCoreString
    def _createPBSCoreString(self):
        """
        Creates the core of a PBS script as a string
        """

        # Shebang line
        jobString = "#!/bin/bash\n"
        # The job name
        jobString += "#PBS -N {}\n".format(self._jobName)
        jobString += "#PBS -l nodes={}:ppn={}\n".format(self._nodes, self._ppn)
        # If walltime is set
        if self._walltime is not None:
            # Wall time, must be in format HOURS:MINUTES:SECONDS
            jobString += "#PBS -l walltime={}\n".format(self._walltime)
        # If submitting to a specific queue
        if self._queue is not None:
            jobString += "#PBS -q {}\n".format(self._queue)
        jobString += "#PBS -o {}.log\n".\
                        format(os.path.join(self._logPath, self._jobName))
        jobString += "#PBS -e {}.err\n".\
                        format(os.path.join(self._logPath, self._jobName))
        if self._account is not None:
            jobString += "#PBS -A {}\n".format(self._account)
        # If we want to be notified by mail
        if self._mail is not None:
            jobString += "#PBS -M {}\n".format(self._mail)
        # #PBS -m abe
        # a=aborted b=begin e=ended
        jobString += "#PBS -m e\n"
        # cd to the folder you are sending the qsub from
        jobString += "cd $PBS_O_WORKDIR\n"

        return jobString
    #}}}

    #{{{_submit
    def _submit(self, jobString, dependentJob=None):
        """
        Saves the jobString as a shell script, submits it and deletes
        it. Returns the output from PBS as a string
        """

        # Create the name of the temporary shell script
        scriptName = "tmp_{}.sh".format(self._time)

        # Save the string as a script
        with open(scriptName, "w") as shell_script:
                shell_script.write(jobString)

        # Submit the jobs
        if dependentJob is None:
            command = "qsub ./{0}".format(scriptName).split(" ")
            completedProcess = run(command, stdout=PIPE, stderr=PIPE)
        else:
            # If the length of the depend job is 0, then all the jobs
            # have completed, and we can carry on as usual without
            # dependencies
            if len(dependentJob) == 0:
                command = "qsub ./"+scriptName
                completedProcess = run(command, stdout=PIPE, stderr=PIPE)
            else:
                # With dependencies
                command = "qsub -W depend=afterok:{} ./{}".\
                          format(dependentJob, scriptName)
                completedProcess = run(command, stdout=PIPE, stderr=PIPE)

        # Check for success
        if completedProcess.returncode != 0:
            print("\nSubmission failed, printing output\n")
            print(completedProcess.stdout)
            print(completedProcess.stderr)
            message = ("The submission failed with exit code {}"
                       ", see the output above").\
                               format(completedProcess.returncode)
            raise RuntimeError(message)

        # Delete the shell script
        try:
            os.remove(scriptName)
        except FileNotFoundError:
            # Do not raise an error
            pass
    #}}}
#}}}

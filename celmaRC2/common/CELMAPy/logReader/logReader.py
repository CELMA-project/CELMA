#!/usr/bin/env python

"""Contains function which reads the log files."""

from collections import OrderedDict
from glob import glob
import numpy as np
import os

#{{{getLogNumbers
def getLogNumbers(path):
    #{{{docstring
    """
    Get the average simulation numbers from the BOUT.log.* files.

    Parameters
    ----------
    path : str
        The path to the log files.

    Returns
    -------
    logNumbers : dict
        Although the keys may vary depending on the settings, the
        standard keys are:
            * SimTime  - The normalized time in the simulation
            * RHSevals - Number of right hand side evaluations before a
                         time step
            * WallTime - Physical time used on the step
            * Calc     - Percentage of time used on arithmetical calculation
            * Inv      - Percentage of time used on laplace inversion
            * Comm     - Percentage of time used on communication
            * I/O      - Percentage of time used on inupt/output
            * SOLVER   - Percentage of time used in the solver
        The values are stored in numpy arrays.
    """
    #}}}

    fileNames = glob(os.path.join(path, "BOUT.log.*"))

    nrFiles = len(fileNames)
    if nrFiles == 0:
        raise RuntimeError("No log files found in {}".format(path))

    # Create the empty dict
    data = createEmptyDict(fileNames)

    for fileName in fileNames:
        # Create a temporary dataDict
        tmpData = OrderedDict({header[0]:[]})
        with open(fileName,"r") as f:
            # Read until the start of the log
            for line in f:
                if "Sim Time  |  RHS evals  | Wall Time |" in line:
                   break

            newlineCount = 0
            for line in f:
                if line == "\n":
                    if newlineCount > 0:
                        # No more time log after this (either error or summary)
                        break
                    else:
                        # First line after header is a newline
                        newlineCount += 1
                        continue

                # Each line ends with a newline character
                line = line.replace("\n", "")
                # Columns separated by spaces
                line = line.split(" ")
                # Get the numbers, not the excess spaces
                line = [column for column in line if column != ""]
                for key, value in zip(tmpData.keys(), line):
                    # Cast to float
                    tmpData[key].append(float(value))

        # Cast to numpy array and add to data dict
        for key in data.keys():
            if data[key] is None:
                data[key] = np.array(tmpData[key])
            else:
                try:
                    data[key] += np.array(tmpData[key])
                except ValueError as e:
                    if "not be broadcast together" in e.args[0]:
                        message = ("The data is corrupted as the number of "
                                   "outputs varies with processor number. "
                                   "Repair by running 'repairBrokenExit'."
                                )
                        raise ValueError(message)
                    else:
                        raise e

    # Divide by number of processors and cast to non-writeable array
    for key in data.keys():
        data[key] = data[key]/nrFiles
        data[key].setflags(write=False)

    return dict(data)
#}}}

#{{{createEmptyDict
def createEmptyDict(fileNames):
    """
    Creates an empty dictionary which will be filled by the getLogNumbers
    routine.

    Parameters
    ----------
    fileNames : list
        List of the file names found with glob.

    Returns
    -------
    data : OrderedDict
        The empty dict to be filled
    """
    with open(fileNames[0],"r") as f:
        # Get the header
        for line in f:
            # Find the start of the timestep log
            if "Sim Time  |  RHS evals  | Wall Time |" in line:
             line = line.replace("\n", "")
             # First part is split by a | character...
             header = line.split("|")
             # ...whilst the last part is split by " "
             headerWithSpaces = header.pop().split(" ")
             header = [*header, *headerWithSpaces]
             header = [head.replace(" ","") for head in header if head != ""]
             data   = OrderedDict({header[0]:[]})
             for head in header[0:]:
                 data[head] = None

             break

    return data
#}}}

#{{{collectiveGetLogNumbers
def collectiveGetLogNumbers(paths):
    #{{{docstring
    """
    Get the merges the simulation numbers for several BOUT.log.0 files.

    NOTE: The paths should be ordered in ascending temporal order.
    NOTE: The 0th iteration of each path will be removed.

    Parameters
    ----------
    paths : tuple
        The path to the log files.

    Returns
    -------
    logNumbers : dict
        Although the keys may vary depending on the settings, the
        standard keys are:
            * SimTime  - The normalized time in the simulation
            * RHSevals - Number of right hand side evaluations before a
                         time step
            * WallTime - Physical time used on the step
            * Calc     - Percentage of time used on arithmetical calculation
            * Inv      - Percentage of time used on laplace inversion
            * Comm     - Percentage of time used on communication
            * I/O      - Percentage of time used on inupt/output
            * SOLVER   - Percentage of time used in the solver
        The values are stored in numpy arrays.
    """
    #}}}

    data = None

    for path in paths:
        curData = getLogNumbers(path)

        if data is None:
            data = curData
            # Remove the first point as this is not counted as a time
            # step
            for key in data.keys():
                data[key]=data[key][1:]
        else:
            for key in data.keys():
                # Remove first point in time in the current time as this
                # is the same as the last of the previous
                data[key]=np.concatenate((data[key], curData[key][1:]), axis=0)

    return data
#}}}

#!/usr/bin/env python

""" Contains repairBrokenExit """

from boututils.datafile import DataFile
import glob, os, shutil
import numpy as np

#{{{repairBrokenExit
def repairBrokenExit(path, mode="rmLastTime"):
    #{{{docstring
    """
    Repairs run where the run was interrupted.

    Will also do a crude check for corruption of the restart files

    Parameters
    ----------
    path : str
        Path to the dmp files to fix.
    mode : ["rmLastTime"|"repair"]
        Which mode to run.
        "rmLastTime" simply removes the last time point for all fields,
        whereas "repair" will try to repair the files by copying from
        the restart files.
    """
    #}}}

    newFiles = glob.glob(os.path.join(path, "BOUT.dmp.*"))
    newFiles.sort()
    nfiles = len(newFiles)

    if nfiles == 0:
        raise ValueError("No data found in {}".format(path))

    print("Number of files found: " + str(nfiles))

    # Make a backup
    bakPath = path + "BAK"
    shutil.move(path, bakPath)
    print("Moved backup to" + bakPath)
    os.mkdir(path)

    # Copy all the files, except the dmp files
    cpyFiles = glob.glob(os.path.join(bakPath, "*"))
    for f in cpyFiles:
        if not(".dmp." in f):
            shutil.copy2(f, f.replace(bakPath, path))

    dmpFiles = glob.glob(os.path.join(bakPath, "BOUT.dmp.*"))
    dmpFiles.sort()

    restartFiles = glob.glob(os.path.join(bakPath, "BOUT.restart.*"))
    restartFiles.sort()
    checkForCorruption(restartFiles)

    if mode == "repair":
        repair(newFiles, dmpFiles, restartFiles)
    elif mode == "rmLastTime":
        rmLastTime(newFiles, dmpFiles)
    else:
        raise NotImplementedError("mode = '{}' is not implemented".format(mode))
#}}}

#{{{checkForCorruption
def checkForCorruption(restartFiles):
    #{{{docstring
    """
    Check for corruption by checking the field mean.

    WARNING: This is not very water-proof.
             Corruption could in theory still have occured, so use with care.

    Paramters
    ---------
    restartFiles : iterable
        Iterable containing the pahts to the restart files.
    """
    #}}}
    print("\nChecking for corrupted restart files by checking the field mean")
    for r in restartFiles:
        print("\nChecking {}".format(r))
        with DataFile(r) as restart:
            for var in restart.list():
                if restart.ndims(var) == 3:
                    mean = restart.read(var).mean()
                    if np.isclose(mean, 0):
                        message="{} has a zero mean. File could be corrupted.".\
                                format(var)
                        raise RuntimeError(message)
                    else:
                        print("{} PASSED with a mean of {}.".format(var, mean))
#}}}

#{{{repair
def repair(newFiles, dmpFiles, restartFiles):
    #{{{docstring
    """
    Will try to repair the corrupted dmpFiles.

    It will do so by:
        1. Copy available data from the restart
        2. Use the second last time point in the dmpFiles if the
        variable is not available in the restartFiles

    Parameters
    ----------
    newFiles : iterable
        Iterable containing the pahts to the new files to be created.
    dmpFiles : iterable
        Iterable containing the pahts to the dump files.
    restartFiles : iterable
        Iterable containing the pahts to the restart files.
    """
    #}}}

    for n, d, r in zip(newFiles, dmpFiles, restartFiles):
        print("\nRepairing {}".format(d))
        # Open the restart file in read mode and create the restart file
        with DataFile(n, write=True, create=True) as newF,\
             DataFile(d) as dmp,\
             DataFile(r) as restart:
            # Check that the problem exists
            time = dmp.read("t_array")
            if not np.isclose(time[-1], 0):
                message = "No problem found in {} as time[-1]={}.\nContinuing".\
                        format(d, time[-1])
                print(message)
                continue
            # Loop over the variables in the dmp file
            for var in dmp.list():
                # Read the data
                dData = dmp.read(var)
                restartList = restart.list()

                # Find 4D variables and time traces
                if dmp.ndims(var) == 4 or dmp.ndims(var)==1:
                    if var in restartList:
                        # Read from the restart and put it back into the dmp file

                        print("    Copying from the restart file "+var)
                        # Read from restart

                        rData   = np.expand_dims(restart.read(var), axis=0)
                        newData = np.concatenate((dData[:-1,...], rData),\
                                                 axis=0)
                    else:
                        if dmp.ndims(var) == 4 or dmp.ndims(var)==1:
                            print(("    Didn't find {} in the restart file, "
                                   "copying second last value").format(var))
                            newData = np.concatenate(\
                                (dData[:-1,...], dData[-1:,...]), axis=0)
                        else:
                            print(("    Didn't find {} in the restart file, "
                                   "copying from the dump").format(var))
                            newData = dData.copy()
                elif var == "iteration":
                    print("    Fixing 'iteration'")
                    newData = dData -1
                else:
                    print("    Nothing to be done for {}, copying from dmp".\
                          format(var))
                    newData = dData.copy()

                newF.write(var, newData)
        print("{} written".format(n))
#}}}

#{{{rmLastTime
def rmLastTime(newFiles, dmpFiles):
    #{{{docstring
    """
    Will remove the last time index in the dump files.

    Parameters
    ----------
    newFiles : iterable
        Iterable containing the pahts to the new files to be created.
    dmpFiles : iterable
        Iterable containing the pahts to the dump files.
    """
    #}}}
    for n, d in zip(newFiles, dmpFiles):
        print("\nRemoving last time point in {}".format(d))
        # Open the restart file in read mode and create the restart file
        with DataFile(n, write=True, create=True) as newF,\
             DataFile(d) as dmp:
            # Loop over the variables in the dmp file

            # Put a 4d variable in the front
            theList = dmp.list()
            ind = theList.index("lnN")
            theList[0], theList[ind] = theList[ind], theList[0]
            for var in theList:
                # Read the data
                dData = dmp.read(var)

                # Find 4D variables and time traces
                if dmp.ndims(var) == 4 or dmp.ndims(var)==1:
                    print("    Removing last time in "+var)
                    # Read from restart
                    newData = dData[:-1,...]
                elif var == "iteration":
                    print("    Fixing 'iteration'")
                    newData = dData -1
                else:
                    print("    Nothing to be done for {}".format(var))
                    newData = dData.copy()

                newF.write(var, newData, info=True)
        print("{} written".format(n))
#}}}

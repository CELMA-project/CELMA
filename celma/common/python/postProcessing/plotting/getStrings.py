#!/usr/bin/env python

"""
Collection of functions returning strings
"""

import datetime
import os

#{{{getTime
def getTime(depth = 'second'):
    """
    Gets the current time, and returns it as a string

    Input
    depth   - String giving the depth of the string
              'hour'   - hour will be the deepest level
              'minute' - minute will be the deepest level
              'second' - second will be the deepest level

    Output
    nowStr  - The string containing the current time
    """
    now = datetime.datetime.now()
    nowStr = "{}-{:02d}-{:02d}".format(now.year, now.month, now.day)

    if depth == 'hour' or depth == 'minute' or depth == 'second':
        nowStr += "-{:02d}".format(now.hour)
    if depth == 'minute' or depth == 'second':
        nowStr += "-{:02d}".format(now.minute)
    if depth == 'second':
        nowStr += "-{:02d}".format(now.second)

    return nowStr
#}}}

#{{{getSaveString
def getSaveString(fileName         ,\
                  simulationPath   ,\
                  timeFolder = None,\
                  prePaths   = None,\
                  postPaths  = None,\
                  ):
    """
    Returns a string which can be used in saving functions.

    Also ensures that the path exisits.

    Input
    fileName       - Name of the file
    simulationPath - Path of the simulation
    timeFolder     - Name of time folder.
                     If none is given, a folder will be made
    prePaths       - Iterable with folder names inserted in the front
    postPaths      - Iterable with folder names inserted in the back

    Output
    saveString - String which can be used to save functions
    timeFolder - String of the time folder
    """

    # Get the folderstructure before the time
    if prePaths == None:
        folders = []
    else:
        if type(prePaths) == str:
            folders = [prePaths]
        else:
            folders = list(prePaths)

    # Create the folder structure
    if timeFolder == None:
        timeFolder = getTime()
        folders.append(timeFolder)
    else:
        folders.append(timeFolder)
    folders.insert(0, simulationPath.split('/')[0])
    savePath = ''
    for folder in folders:
        savePath = os.path.join(savePath, folder)

    # Get the folderstructure after the time
    if postPaths == None:
        folders = []
    else:
        if type(postPaths) == str:
            folders = [postPaths]
        else:
            folders = list(postPaths)
    for folder in folders:
        savePath = os.path.join(savePath, folder)

    # Make dir if not exists
    if not os.path.exists(savePath):
        os.makedirs(savePath)

    # Append the fileName
    saveString = os.path.join(savePath, fileName)

    return saveString, timeFolder
#}}}

#!/usr/bin/env python

"""
Clollection of savePathFuncs
"""

#{{{scanWTagSaveFunc
def scanWTagSaveFunc(dmpFolder, theRunName = None):
    """
    Takes the tag and run name into account in the saveFolder

    Parameters
    ----------
    dmpFolder : [tuple|str]
        String with dump location
    theRunName : str
        The run name
    """
    if theRunName is None:
        message = "'theRunName' required as a keyword argument in scanWTagSaveFunc"
        raise RuntimeError(message)

    if hasattr(dmpFolder, "__iter__") and type(dmpFolder) != str:
        dmpFolder = dmpFolder[0]

    # Splits / and chooses the path with the scan parameter
    pathWTag = None
    allFolders = dmpFolder.split("/")
    for folder in allFolders:
        if "tag" in folder:
            pathWTag = folder
            break
    if pathWTag is None:
        message = "'tag' not found in {}".format(dmpFolder)
        raise RuntimeError(message)

    saveFolder =\
        theRunName + "-" +\
        pathWTag.replace("_tag_"+theRunName+"_0","")

    return saveFolder
#}}}

#{{{onlyScan
def onlyScan(dmpFolder, scanParameter=None):
    """
    Make folder with only the scan parameter and its value

    Parameters
    ----------
    dmpFolder : [tuple|str]
        String with dump location
    scanParameter : str
        The scan parameter
    """

    if scanParameter is None:
        message = "'scanParameter' required as a keyword argument in onlyScan"
        raise RuntimeError(message)

    if hasattr(dmpFolder, "__iter__") and type(dmpFolder) != str:
        dmpFolder = dmpFolder[0]

    # Splits / and chooses the path with the scan parameter
    pathWScanParameter = None
    allFolders = dmpFolder.split("/")
    for folder in allFolders:
        if scanParameter in folder:
            pathWScanParameter = folder
            break
    if pathWScanParameter is None:
        message = "Scan parameter {} not found in {}".\
            format(dmpFolder, scanParameter)
        raise RuntimeError(message)

    listScanPath = pathWScanParameter.split("_")
    # Find the scan parameter index
    ind = listScanPath.index(scanParameter)
    # The value is found in the preceding index
    val = listScanPath[ind+1]

    saveFolder = "_".join((scanParameter, val))

    return saveFolder
#}}}

#!/usr/bin/env python

"""
Contains the path merger which merges paths of same scan parameters
"""

#{{{pathMerger
def pathMerger(dmpFoldersDict, mergeThese):
    #{{{docstring
    """
    Merge the dmpFolders into a tuple of tuple containing paths of same
    scan parameters in ascending temporal order.

    Parameters
    ----------
    dmpFoldersDict : dict
        Dictionary pickled by the scanDriver
    mergeThese : tuple
        Tuple of keys to merge

    Returns
    -------
    mergedDict : dict
        Dictionary with keys "param0", "param1" etc., where the values
        are the paths of same scan parameters in ascending temporal order.
    """
    #}}}
    numberOfParameters = len(dmpFoldersDict["init"])
    mergedDicts = {"param{}".format(nr):[] for nr in range(numberOfParameters)}

    # NOTE: This may have DRY problems
    if "init" in mergeThese:
        for nr in range(numberOfParameters):
            mergedDicts["param{}".format(nr)].\
                    append(dmpFoldersDict["init"][nr])

    if "expand" in mergeThese:
        for nr in range(numberOfParameters):
            mergedDicts["param{}".format(nr)].\
                    append(dmpFoldersDict["expand"][nr])

    if "linear" in mergeThese:
        for nr in range(numberOfParameters):
            mergedDicts["param{}".format(nr)].\
                    append(dmpFoldersDict["linear"][nr])

    if "turbulence" in mergeThese:
        for nr in range(numberOfParameters):
            mergedDicts["param{}".format(nr)].\
                    append(dmpFoldersDict["turbulence"][nr])

    if "extraTurbulence" in mergeThese:
        for nr in range(numberOfParameters):
            for restartNr in range(len(dmpFoldersDict["extraTurbulence"][nr])):
                mergedDicts["param{}".format(nr)].\
                    append(dmpFoldersDict["extraTurbulence"][nr][restartNr])

    for key in mergedDicts.keys():
        mergedDicts[key] = tuple(mergedDicts[key])

    return mergedDicts
#}}}

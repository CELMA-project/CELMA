#!/usr/bin/env python

"""
Clollection of savePathFuncs
"""

#{{{scanWTagSaveFunc
def scanWTagSaveFunc(dmpFolder, theRunName = None):
    """
    Function which makses the saveFolder string (typically used for post
    processing) when a tag is given.

    Parameters
    ----------
    dmpFolder : [tuple|str]
        String with dump location
    theRunName : str
        The run name
    """

    if hasattr(dmpFolder, "__iter__") and type(dmpFolder) != str:
        dmpFolder = dmpFolder[0]

    saveFolder =\
        theRunName + "-" +\
        dmpFolder.split('/')[-1].replace('_tag_'+theRunName+'_0','')

    return saveFolder
#}}}

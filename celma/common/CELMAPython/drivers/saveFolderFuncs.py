#!/usr/bin/env python

"""
Clollection of saveFolderFuncs for driversCombined.py
"""

#{{{scanWTagSaveFunc
def scanWTagSaveFunc(path, theRunName = None, *args, **kwargs):
    """Function which makses the saveFolder string when a tag is given"""

    if hasattr(path, "__iter__"):
        path = path[0]

    saveFolder = theRunName + "-" + path.split('/')[-1].replace('_tag_'+theRunName+'_0','')

    return saveFolder
#}}}

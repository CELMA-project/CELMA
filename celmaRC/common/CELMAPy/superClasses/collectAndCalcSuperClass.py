#!/usr/bin/env python

"""
Contains the super class for the collect and calc
"""

from ..unitsConverter import UnitsConverter

#{{{CollectAndCalcSuperClass
class CollectAndCalcSuperClass(object):
    """
    The parent collect and calc class
    """

    #{{{Constructor
    def __init__(self                     ,\
                 dmp_folders              ,\
                 collectPaths      = None ,\
                 convertToPhysical = False,\
                 uc                = None ,\
                 ):
        #{{{docstring
        """
        This constructor

        * Sets the dmp_folders, collectPaths and convertToPhysical
        * Creates a UnitsConverter if none is given

        Parameters
        ----------
        dmp_folders: tuple
            This is the output dmp_folder from bout_runners.
            Typically, these are the folders in a given scan
        collectPaths : [None|tuple]
            Paths to collect from.
            If None dmp_folders will be set to collectPaths
        convertToPhysical : bool
            If the physical or normalized units should be plotted.
        uc : [None|UnitsConverter]
            The UnitsConverter will be set from the first element in
            collectPaths if not given.
        """
        #}}}

        # Set the member data
        self.dmp_folders = dmp_folders

        # Set the collection paths
        if collectPaths is None:
            collectPaths  = tuple(dmp_folders)
        else:
            collectPaths = collectPaths

        self._collectPaths = collectPaths

        # Make the unitsConverter
        if uc is None:
            self.uc = UnitsConverter(collectPaths[0], convertToPhysical)

        self.convertToPhysical = self.uc.convertToPhysical
#}}}

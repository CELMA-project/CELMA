#!/usr/bin/env python

"""
Contains the super class driver for points
"""

from .driverSuperClass import DriverSuperClass

#{{{DriverPointsSuperClass
class DriverPointsSuperClass(DriverSuperClass):
    """
    The parent driver for points
    """

    #{{{Constructor
    def __init__(self                    ,\
                 *args                   ,\
                 convertToPhysical = True,\
                 **kwargs):
        #{{{docstring
        """
        This constructor:
            * Calls the parent class

        Parameters
        ----------
        *args : positional arguments
            See parent class for details.
        convertToPhysical : bool
            Whether or not to convert to physical units.
        **kwargs : keyword arguments
            See parent class for details.
        """
        #}}}

        # Call the constructors of the parent class
        super().__init__(*args, **kwargs)

        # Set member data
        self.convertToPhysical = convertToPhysical
    #}}}
#}}}

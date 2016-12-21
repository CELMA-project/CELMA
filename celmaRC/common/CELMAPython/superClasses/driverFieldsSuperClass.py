#!/usr/bin/env python

"""
Contains the super class driver for fields 1D and fields 2D
"""

from .postProcessorDriver import PostProcessorDriver

#{{{DriverFieldsSuperClass
class DriverFieldsSuperClass(PostProcessorDriver):
    """
    The parent driver of 1D and 2D field plotting
    """

    #{{{Constructor
    def __init__(self                        ,\
                 *args                       ,\
                 collectPaths = None         ,\
                 xguards      = False        ,\
                 yguards      = False        ,\
                 xSlice       = slice(0,None),\
                 ySlice       = slice(0,None),\
                 zSlice       = slice(0,None),\
                 tSlice       = None         ,\
                 **kwargs):
        #{{{docstring
        """
        This constructor:
            * Calls the parent class
            * Sets the common memberdata

        Parameters
        ----------
        *args : str
            See parent class for details.
        collectPaths : tuple
            Paths to collect from.
            The corresponind 't_array' of the paths must be in ascending order.
        xguards : bool
            If xguards should be included when collecting.
        yguards : bool
            If yguards should be included when collecting.
        xSlice : slice
            How the data will be sliced in x.
        ySlice : slice
            How the data will be sliced in y.
        zSlice : slice
            How the data will be sliced in z.
        tSlice : slice
            How the data will be sliced in t.
        writer : str
            Writer to use if the plots are animated.
        **kwargs : keyword arguments
            See parent class for details.
        """
        #}}}

        # Call the constructor of the parent class
        super().__init__(*args, **kwargs)

        # Set the member data
        self._collectPaths = collectPaths
        self._xguards      = xguards
        self._yguards      = yguards
        self._xSlice       = xSlice
        self._ySlice       = ySlice
        self._zSlice       = zSlice
        self._tSlice       = tSlice
    #}}}
#}}}


# # FIXME: Move maxGradRho to field2D class, and make it a ghost there
#         # Get the current scan
#         if maxGradRhoFolder:
#             if self._scanParameters:
#                 self._maxGradRhoFolder = maxGradRhoFolder
#             else:
#                 self._maxGradRhoFolder =\
#                     convertToCurrentScanParameters(dmpFolder, maxGradRhoFolder, scanParameters)
#         else:
#             self._maxGradRhoFolder = None

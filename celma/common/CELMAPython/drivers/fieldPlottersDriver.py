#!/usr/bin/env python

"""
Contains the post processing driver class for plotting the fields
"""

from .postProcessorDriver import PostProcessorDriver

#{{{FieldPlottersDriver
class FieldPlottersDriver(PostProcessorDriver):
    """
    The parent driver of all field plotting functions

    Sets the memberdata
    """

    #{{{Constructor
    def __init__(self                             ,\
                 *args                            ,\
                 xguards           = False        ,\
                 yguards           = False        ,\
                 xSlice            = slice(0,None),\
                 ySlice            = slice(0,None),\
                 zSlice            = slice(0,None),\
                 tSlice            = None         ,\
                 **kwargs
                ):
        #{{{docstring
        """
        This constructor sets the memberdata.

        Parameters
        ----------
        *args : str
            See PostProcessorDriver for more details.
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
        **kwargs : keyword arguments
            Additional keyword arguments given as input to
            PostProcessorDriver.
        """
        #}}}

        # Call the constructor of the parent class
        super().__init__(*args, **kwargs)

        # Set the member data
        self._xguards = xguards
        self._yguards = yguards
        self._xSlice  = xSlice
        self._ySlice  = ySlice
        self._zSlice  = zSlice
        self._tSlice  = tSlice
#}}}

#!/usr/bin/env python

"""
Contains the Line class.
"""

#{{{class Line
class Line(object):
    """
    Class containing the data for the line which are going to be
    plotted.

    Provides the data which is connected to a line in a subplot.
    Specifically it contains information about:

    * name     - The name used for collecting the data from a simulation
    * label    - The label which is going to be used for in the legend
    * plotPos  - Index of the plot number (if any)
    * ax       - The axis object for the line
    * bottomAx - If the axis is the lowest
    * lineObj  - The lineobject of the plot
    * color    - The color of the line
    * field    - The data used to plot the line
    * dim      - The dimension of the data
    """

    #{{{Constructor
    def __init__(self          ,\
                 name          ,\
                 label         ,\
                 plotPos = None,\
                ):
        #{{{docstring
        """
        The constructor sets the member data

        Parameters
        ----------
        name : str
            Name which is going to be collected
        label : str
            Label which is going to be used in the plot
            NOTE: "$" will be added around this string
            NOTE: Remember to use raw strings
        plotPos : int
            If there is a preferred position of the plot.
            Given as an index number.
            NOTE: The user have to take care so that two plots
            does not share the same number.
        """
        #}}}

        # Set member data from input
        self.name    = name
        self.label   = r"$" + label + r"$"
        self.plotPos = plotPos

        # Initialize non-input members
        self.ax       = None
        self.bottomAx = False
        self.lineObj  = None
        self.color    = None
        self.field    = None
        self.dim      = None
    #}}}
#}}}

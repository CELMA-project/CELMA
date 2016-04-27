#!/usr/bin/env python

from .line import Line
import numpy as np
from matplotlib.gridspec import GridSpec
import matplotlib.cm as cm
import matplotlib.pyplot as plt

"""
Contains the organizer class
"""

#{{{class Organizer
class Organizer(object):
    """
    Class which organizes several lines in a plot.

    This class is responsible for

    * Organization of the lines
    * Make all axes
    * Setting proper names on figures and files
    * Making a combined line if useCombinedPlot is true
    """

    #{{{Constructor
    def __init__(self                   ,\
                 pltName                ,\
                 cols            = 2    ,\
                 useCombinedPlot = False,\
                 ):
        """
        The constructor initializes the list of lines

        pltName         - Name of the plot written in LaTeX format, but
                          without the $
        cols            - The total number of columns to be used in the plot
        useCombinedPlot - Toggles if a plot of combined lines are to be plotted
        """

        # Set member data from input
        self._cols           = cols
        self.useCombinedPlot = useCombinedPlot
        self.pltName         = pltName

        # Initialize non-input members
        self._pltSize             = (18,12)
        self.combLine             = None
        self.combLineLineObjs     = []
        self.lines                = []
        self.axes                 = []
    #}}}

    #{{{pltPrepare
    def pltPrepare(self):
        """
        Prepares the lines in a plot for plotting.
        Call this function before calling collect.

        1. Check if any plot pos have been given.
           If yes, lines will be rearranged.
        2. Set the color of each plot.
        3. Finds the bottom axes
        4. Creates the figure and axes.
        5. Returns the figure
        """

        # Organize the lines
        newLines = self.lines.copy()
        for line in self.lines:
            if line.plotPos:
                newInd = line.plotPos
                oldInd = self.lines.index(line)
                # Swap
                newLines[oldInd], newLines[newInd] =\
                        newLines[newInd], newLines[oldInd]

        # Reassign
        self.lines = newLines

        # Set the colors
        colorSpace = np.arange(len(self.lines))
        colors = cm.rainbow(np.linspace(0, 1, len(colorSpace)))

        for lineNr, line in enumerate(self.lines):
            line.color = colors[lineNr]

        # If a combined line is to be plotted
        if self.useCombinedPlot:
            # Make a line object
            self.combLine = Line(name  = 'combinedLine'               ,\
                                 label = r'\partial_t ' + self.pltName,\
                                 )
            # Make the lastline black, and append it to the lines
            self.combLine.color = 'k'
            self.lines.append(self.combLine)

        # Calculate the number of rows
        rows = int(np.ceil(len(self.lines)/self._cols))

        # Create the figure
        fig = plt.figure(figsize = self._pltSize)
        gs  = GridSpec(rows, self._cols)

        # Make the axes
        for lineNr, line in enumerate(self.lines):
            if lineNr == 0:
                # Need an initial line
                line.ax = plt.subplot(gs[lineNr])
                firstAx = line.ax
            else:
                line.ax = plt.subplot(gs[lineNr], sharex=firstAx)

        for col in range(1, self._cols+1):
            self.lines[-col].bottomAx = True

        # Pop the combined line in order not to collect it
        if self.useCombinedPlot:
            self.combLine = self.lines.pop()

        fig.canvas.set_window_title(self.pltName)

        return fig
    #}}}

    #{{{makeCombinedLine
    def makeCombinedLine(self):
        """
        Makes a combined line.
        To be called after all other lines are collected.
        """
        # Initialize the field
        self.combLine.field = np.zeros(self.lines[0].field.shape)

        for line in self.lines:
            self.combLine.field += line.field

        # Re-add the combLine to the list
        self.lines.append(self.combLine)
    #}}}
#}}}

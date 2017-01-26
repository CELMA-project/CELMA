#!/usr/bin/env python

""" Contains the SizeMaker class """

import matplotlib.pyplot as plt

#{{{SizeMaker
class SizeMaker(object):
    r"""
    Class which determines the size (in inches) and calculates the dpi.

    NOTE: For the text size to remain constant one must use
          {s\textwidth} in the plots in LaTeX, where 's' is the scale
          set in the function.
    """

    # \textwidth
    # \usepackage{layouts}
    # \printinunitsof{in}\prntlen{\textwidth}
    textWidth = 6.3
    # Golden ratio of heigth/width
    aspect = 0.7
    # Dots (default from matplotlib dpi = 100 for 6.8 inches)
    dots = 1000

    @staticmethod
    #{{{setTextWidth
    def setTextWidth(w):
        """
        Sets the static textWidth (in inches) with w
        """
        SizeMaker.textWidth = w
    #}}}

    @staticmethod
    #{{{setAspect
    def setAspect(a):
        """
        Sets the static aspect with a
        """
        SizeMaker.aspect = a
    #}}}

    @staticmethod
    #{{{setDots
    def setDots(d):
        """
        Sets the static dots (used to calculate dpi) with d
        """
        SizeMaker.dots = d
    #}}}

    @staticmethod
    #{{{standard
    def standard(s = 1, w = textWidth, a = aspect):
        #{{{docstring
        """
        Returns the standard plot size, recalculates the dpi.

        Parameters
        ----------
        s : float
            The scale of the width and heigth.
        w : float
            The width (in inches).
        a : float
            The aspec ratio heigth/width.

        Returns
        -------
        plotSize : tuple
            The plot size as (width, heigth) in inches.
        """
        #}}}

        w *= s
        plt.rc("figure", dpi = SizeMaker.dots/w)
        return (w, w*a)
    #}}}

    @staticmethod
    #{{{golden
    def golden(s = 1, w = textWidth, a = 1/1.618):
        #{{{docstring
        """
        Returns the golden ratio plot size, recalculates the dpi.

        Parameters
        ----------
        s : float
            The scale of the width and heigth.
        w : float
            The width (in inches).
        a : float
            The aspec ratio heigth/width.

        Returns
        -------
        plotSize : tuple
            The plot size as (width, heigth) in inches.
        """
        #}}}

        w *= s
        plt.rc("figure", dpi = SizeMaker.dots/w)
        return (w, w*a)
    #}}}

    @staticmethod
    #{{{square
    def square(s = 1, w = textWidth):
        #{{{docstring
        """
        Returns the square plot size, recalculates the dpi.

        Parameters
        ----------
        s : float
            The scale of the width and heigth.
        w : float
            The width (in inches).

        Returns
        -------
        plotSize : tuple
            The plot size as (width, heigth) in inches.
        """
        #}}}

        w *= s
        plt.rc("figure", dpi = SizeMaker.dots/w)
        return (w, w)
    #}}}

    @staticmethod
    #{{{array
    def array(cols, rows, wSpace=0.0, hSpace=0.0, w = textWidth, aSingle = aspect):
        #{{{docstring
        """
        Returns the plot size for an array plot, recalculates the dpi.

        Parameters
        ----------
        cols : int
            Number of columns in the array.
        rows : int
            Number of rows in the array.
        wSpace : float
            Spacing between the columns (in inches).
            Setting this will make the the total height more correct.
        hSpace : float
            Spacing between the rows (in inches).
            Setting this will make the the total height more correct.
        w : float
            Width of the entire figure.
        aSingle : float
            Aspect for a single plot in the array.

        Returns
        -------
        plotSize : tuple
            The plot size as (width, heigth) in inches.
        """
        #}}}

        wSingle = (w/cols) - wSpace
        hSingle = aSingle*wSingle
        h       = (hSingle + hSpace)*rows
        plt.rc("figure", dpi = SizeMaker.dots/w)
        return (w, h)
    #}}}
#}}}

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
    textwidth = 6.3
    # Golden ratio of heigth/width
    aspect = 0.7
    # Dots (default from matplotlib dpi = 100 for 6.8 inches)
    dots = 1000

    @staticmethod
    #{{{setTextWidth
    def setTextWidth(w):
        """
        Sets the static textwidth with w
        """
        SizeMaker.textwidth = w
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
    #{{{standard
    def standard(s = 1, w = textwidth, a = aspect):
        #{{{docstring
        """
        Returns the standard plot size

        Parameters
        ----------
        s : float
            The scale of the width and heigth.
        w : float
            The width.
        a : float
            The aspec ratio heigth/width.

        Returns
        -------
        plotSize : tuple
            The plot size as (width, heigth)
        """
        #}}}

        w *= s
        plt.rc("figure", dpi = SizeMaker.dots/w)
        return (w, w*a)
    #}}}

    @staticmethod
    #{{{golden
    def golden(s = 1, w = textwidth, a = 1/1.618):
        #{{{docstring
        """
        Returns the golden ratio plot size

        Parameters
        ----------
        s : float
            The scale of the width and heigth.
        w : float
            The width.
        a : float
            The aspec ratio heigth/width.

        Returns
        -------
        plotSize : tuple
            The plot size as (width, heigth)
        """
        #}}}

        w *= s
        plt.rc("figure", dpi = SizeMaker.dots/w)
        return (w, w*a)
    #}}}

    @staticmethod
    #{{{square
    def square(s = 1, w = textwidth):
        #{{{docstring
        """
        Returns the square plot size

        Parameters
        ----------
        s : float
            The scale of the width and heigth.
        w : float
            The width.

        Returns
        -------
        plotSize : tuple
            The plot size as (width, heigth)
        """
        #}}}

        w *= s
        plt.rc("figure", dpi = SizeMaker.dots/w)
        return (w, w)
    #}}}
#}}}

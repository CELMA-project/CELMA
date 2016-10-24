#!/usr/bin/env python

"""
Contains drivers for plotting 1D and 2D plots
"""

from .drivers1D import Drivers1D
from .drivers2D import Drivers2D
from multiprocessing import Process

#{{{Drivers1D2D
class Drivers1D2D(Drivers1D, Drivers2D):
    """
    Class which combines the 1D and 2D plots of the fields.
    """

    #{{{Constructor
    def __init__(self             ,\
                 *args            ,\
                 **kwargs):
        #{{{docstring
        """
        This constructor:

        1. Calls the parent constructor

        Parameters
        ----------
        *args : positional arguments
            See the constructor of Drivers1D and Drivers2D for details.
        **kwargs : keyword arguments
            See the constructor of Drivers1D and Drivers2D for details.
        """
        #}}}

        # Call the constructor of the parent class
        super(Drivers1D2D, self).__init__(*args, **kwargs)
    #}}}

    #{{{plot1DAnd2DDriver
    def plot1DAnd2DDriver(self):
        if self._useSubProcess:
            #{{{ Function call through subprocess
            proc = {}
            # Plot the 2D plots
            proc["allMainFields2DDriver"] =\
                Process(\
                        target = self.allMainFields2DDriver,\
                        args   = ()                        ,\
                        kwargs = {}
                       )

            # Plot the 1D plots
            proc["parPerpDriver"] =\
                Process(\
                        target = self.parPerpDriver,\
                        args   = ()                ,\
                        kwargs = {}
                       )

            # Start the processes
            for key in proc.keys():
                proc[key].start()
            # Wait for process to finish
            for key in proc.keys():
                proc[key].join()
            #}}}
        else:
            #{{{ Normal function call
            # Plot the 2D plots
            self.allMainFields2DDriver()

            # Plot the 1D plots
            self.parPerpDriver()
            #}}}

    #}}}

    #{{{plot1D2DAndFluctDriver
    def plot1D2DAndFluctDriver(self):
        if self._useSubProcess:
            #{{{ Function call through subprocess
            # Plot original fields
            self._setSubPolAvg(False)
            Process(\
                    target = self.plot1DAnd2DDriver,\
                    args   = ()                    ,\
                    kwargs = {}
                   ).start()

            # Plot fluctuation fields
            self._setSubPolAvg(True)
            Process(\
                    target = self.plot1DAnd2DDriver,\
                    args   = ()                    ,\
                    kwargs = {}
                   ).start()
            #}}}
        else:
            #{{{ Normal function call
            # Plot original fields
            self._setSubPolAvg(False)
            self.plot1DAnd2DDriver()

            # Plot fluctuation fields
            self._setSubPolAvg(True)
            self.plot1DAnd2DDriver()
            #}}}
    #}}}
#}}}

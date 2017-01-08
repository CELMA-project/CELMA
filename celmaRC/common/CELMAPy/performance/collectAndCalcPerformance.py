#!/usr/bin/env python

"""
Contains class for collecting and calculating the performance
"""

from ..collectAndCalcHelpers import collectTime
from ..logReader import collectiveGetLogNumbers
from ..superClasses import CollectAndCalcSuperClass
import numpy as np

#{{{CollectAndCalcPerformance
class CollectAndCalcPerformance(CollectAndCalcSuperClass):
    """
    Class for collecting and calcuating the performance
    """

    #{{{constructor
    def __init__(self   ,\
                 *args  ,\
                 **kwargs):
        #{{{docstring
        """
        This constructor will:
            * Call the parent constructor

        Parameters
        ----------
        *args : positional arguments
            See parent constructor for details.
        *kwargs : keyword arguments
            See parent constructor for details.
        """
        #}}}

        # Call the constructor of the parent class
        super().__init__(*args, **kwargs)
    #}}}

    #{{{executeCollectAndCalc
    def executeCollectAndCalc(self):
        #{{{docstring
        """
        Function which collects and calculates the performance.

        Returns
        -------
        performance : dict
            Although the keys may vary depending on the settings, the
            standard keys are:
                * Sim Time  - The normalized time in the simulation
                * RHS evals - Number of right hand side evaluations before a
                              time step
                * Wall Time - Physical time used on the step
                * Calc      - Percentage of time used on arithmetical
                              calculation
                * Inv       - Percentage of time used on laplace inversion
                * Comm      - Percentage of time used on communication
                * I/O       - Percentage of time used on inupt/output
                * SOLVER    - Percentage of time used in the solver
            The values are stored in tuples.
        """
        #}}}

        performance = collectiveGetLogNumbers(self._collectPaths)

        # Get the dt as we would like to plot RHS evals/timestep
        dt = None
        for path in self._collectPaths:
            time = collectTime((path,))
            # Minus one as the first step will have no dt
            curDt = np.empty(len(time)-1)
            for i in range(len(time)-1):
                curDt[i] = time[i+1] - time[i]
            if dt is None:
                dt = curDt
            else:
                dt=np.concatenate((dt, curDt), axis=0)

        # Calc RHSPrTime
        performance["RHSPrTime"] = performance.pop("RHSevals")/dt

        # Rename Calc
        performance["Arithmetic"] = performance.pop("Calc")

        # Convert the sim time, and rename it to time
        if self.uc.convertToPhysical:
            performance["time"]  =\
                self.uc.physicalConversion(performance["SimTime"], "t")
        else:
            performance["time"] = performance["SimTime"]

        performance.pop("SimTime")

        return performance
    #}}}
#}}}

#!/usr/bin/env python

"""
Contains function which collects and calculates the energies
"""

from ..superClasses import CollectAndCalcSuperClass
from ..collectAndCalcHelpers import collectiveCollect, collectTime

#{{{CollectAndCalcEnergy
class CollectAndCalcEnergy(CollectAndCalcSuperClass):
    """
    Class for collecting and calcuating the energies
    """

    #{{{constructor
    def __init__(self   ,\
                 *args  ,\
                 **kwargs):
        #{{{docstring
        """
        This constructor will:
            * Call the parent constructor
            * Initialized the timeSlice

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

        self.setTSlice()
    #}}}

    #{{{executeCollectAndCalc
    def executeCollectAndCalc(self):
        #{{{docstring
        """
        Function which collects and calculates the energies.

        Returns
        -------
        energies : dict
            Dictionary with the keys:
            * "perpKinEE" - Time trace of perpendicular kinetic electron energy
            * "parKinEE"  - Time trace of parallel kinetic electron energy
            * "perpKinEI" - Time trace of perpendicular kinetic ion energy
            * "parKinEI"  - Time trace of parallel kinetic ion energy
            * "polEE"     - Time trace of potential electron energy
            * "time"      - Time corresponding time
        """
        #}}}

        # Initialize output
        eKeys   = ("perpKinEE", "parKinEE")
        iKeys   = ("perpKinEI", "parKinEI")
        energies = {}

        # Set the slice
        if self._tSlice is not None:
            tStart = self._tSlice.start
            tEnd   = self._tSlice.stop
            tInd = (tStart, tEnd)
        else:
            tInd = None
        if self._tSlice is not None:
            # Slice the variables with the step
            # Make a new slice as the collect dealt with the start and
            # the stop of the slice
            newSlice = slice(None, None, self._tSlice.step)

        # Collect the energies
        for key in eKeys:
            var =\
                collectiveCollect(self._collectPaths, (key,), tInd = tInd)
            var = var[key][:,0,0,0]
            if self._tSlice is not None:
                var = var[newSlice]
            if self.uc.convertToPhysical:
                energies[key] = self.uc.physicalConversion(var, "eEnergy")
            else:
                energies[key] = self.uc.normalizedConversion(var, "eEnergy")

        for key in iKeys:
            var =\
                collectiveCollect(self._collectPaths, (key,), tInd = tInd)
            var = var[key][:,0,0,0]
            if self._tSlice is not None:
                var = var[newSlice]
            if self.uc.convertToPhysical:
                energies[key] = self.uc.physicalConversion(var, "iEnergy")
            else:
                energies[key] = var

        # Special treatment of the potential
        key = "particleNumber"
        var =\
            collectiveCollect(self._collectPaths, (key,), tInd = tInd)
        var = var[key][:,0,0,0]
        if self._tSlice is not None:
            var = var[newSlice]
        # NOTE: Te is a free variable, when normalized, it equals 1
        #       Hence the normalized potential energy equals the
        #       normalized particle number
        if self.uc.convertToPhysical:
            # NOTE: The eEnergy normalization and iEnergy normalization
            #       are the same apart from a factor 1/mu
            #       Neither the potential energy nT, nor the iEnergy has
            #       the factor, thus we use the iEnergy factor here
            energies["potEE"] = self.uc.physicalConversion(var, "iEnergy")
        else:
            energies["potEE"] = var

        # Collect the time
        time = collectTime(self._collectPaths, tInd=tInd)
        if self._tSlice is not None:
            time = time[newSlice]
        if self.uc.convertToPhysical:
            energies["time"] = self.uc.physicalConversion(time, "t")
        else:
            energies["time"] = time

        return energies
    #}}}

    #{{{setTSlice
    def setTSlice(self, tSlice=None):
        #{{{docstring
        """
        Sets the indices

        Parameters
        ----------
        tSlice : [None|slice]
            If given this is the slice of t to use when collecting.
        """
        #}}}

        # Set the tSlice
        self._tSlice = tSlice
    #}}}
#}}}

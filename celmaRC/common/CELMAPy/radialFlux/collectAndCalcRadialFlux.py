#!/usr/bin/env python

"""
Contains the radial flux calculation
"""

from ..calcVelocities import calcRadialExBPoloidal
from ..collectAndCalcHelpers import slicesToIndices

#{{{CollectAndCalcRadialFlux
class CollectAndCalcRadialFlux(object):
    """
    Class for collecting and calcuating the radial flux
    """

    #{{{constructor
    def __init__(self             ,\
                 collectPaths     ,\
                 slices           ,\
                 mode             ,\
                 dh               ,\
                 convertToPhysical,\
                 ):
        #{{{docstring
        """
        This constructor will:
            * Call the parent constructor

        Parameters
        ----------
        collectPaths : tuple
            Tuple from where to collect
        slices : tuple of tuples
            Tuple the indices to use.
            On the form ((xInd1,...), (yInd1,...), (zInd1,...), (tSlice1,...))
        mode : ["normal"|"fluct"]
            Whether to look at fluctuations or normal data
        dh : DimensionsHelper
            DimensionHelper object (used to find the Jacobian)
        convertToPhysical : bool
            Whether or not to convert to physical
        """
        #}}}

        # Set the member data
        self._collectPaths      = collectPaths
        self._mode              = mode
        self._dh                = dh
        self._convertToPhysical = convertToPhysical
        self._xInds, self._yInds, self._zInds, self._tSlice = slices
    #}}}

    #{{{getRadialExBTrace
    def getRadialExBTrace(self):
        #{{{docstring
        """
        Function which calculates the radial flux.

        NOTE: We are here translating a Field1D to a time trace.

        Returns
        -------
        radialExBTraces : dict
            Dictionary where the keys are on the form "rho,theta,z".
            The value is a dict containing of
            {"uExBRadial":radialExBtimeTrace, "time":time}.
            And additional key "zInd" will be given in addition to varName
            and "time" if mode is set to "fluct".
            The timeTrace is a 1d array.
        """
        #}}}

        # Initialize output
        raidalExBTraces = {}
        tCounter = 0
        for x, y, z in zip(self._xInds, self._yInds, self._zInds):
            # NOTE: The indices
            rho   = self._dh.rho     [x]
            theta = self._dh.thetaDeg[z]
            par   = self._dh.z       [y]

            # Add key and dict to timeTraces
            key = "{},{},{}".format(rho,theta,par)
            raidalExBTraces[key] = {}

            # Collect and slice
            if self._tSlice is not None:
                t = self._tSlice[tCounter]
            tCounter += 1

            slices = (x, y, t)
            radialExB, time = calcRadialExBPoloidal(\
                                self._collectPaths                         ,\
                                slices                                     ,\
                                mode = self._mode                          ,\
                                convertToPhysical = self._convertToPhysical,\
                               )

            # Save and cast to to 2d:
            raidalExBTraces[key]["uExBRadial"] = radialExB[:, 0, 0, z]
            raidalExBTraces[key]["time"] = time

        return raidalExBTraces
    #}}}

    @staticmethod
    #{{{calcRadialFlux
    def calcRadialFlux(timeTraces, radialExBTraces):
        #{{{docstring
        """
        Function which calculates the radial flux.

        Parameters
        ----------
        timeTraces : dict
            Dictionary where the keys are on the form "rho,theta,z".
            The value is a dict containing of
            {varName:timeTrace, "time":time}.
            And additional key "zInd" will be given in addition to varName
            and "time" if mode is set to "fluct".
            The timeTrace is a 1d array.
        raidalExBTraces : dict
            The same as the timeTraces, but where the varName is
            "uExBRadial".

        Returns
        -------
        radialFlux : dict
            The same as the timeTraces, but where the varName is
            varNameRadialFlux.
        """
        #}}}

        # Initialize the output
        radialFlux = {}

        # Obtain the varName
        ind  = tuple(timeTraces.keys())[0]
        keys = timeTraces[ind].keys()
        varName = tuple(var for var in keys if var != "time")[0]

        # Make the keys
        fluxKey = "{}RadialFlux".format(varName)

        # Obtain the radialFlux
        for key in timeTraces.keys():
            # Initialize the radialFlux
            radialFlux[key] = {}

            radialFlux[key]["time"]  = timeTraces[key]["time"]
            radialFlux[key][fluxKey] =\
                timeTraces[key][varName] * radialExBTraces[key]["uExBRadial"]

        # NOTE: If timeTraces was converted to physical units, then radialFlux
        #       is in physical units as well
        return radialFlux
    #}}}
#}}}

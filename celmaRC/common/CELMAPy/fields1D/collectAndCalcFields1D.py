#!/usr/bin/env python

"""
Contains class for collecting and calculating the 1D fields
"""

from ..superClasses import CollectAndCalcFieldsSuperClass
from ..collectAndCalcHelpers import (collectTime,\
                                     collectConstRho,\
                                     collectConstZ,\
                                     collectParallelProfile,\
                                     collectRadialProfile,\
                                     collectPoloidalProfile,\
                                     polAvg,\
                                     slicesToIndices,\
                                     timeAvg)

#{{{CollectAndCalcFields1D
class CollectAndCalcFields1D(CollectAndCalcFieldsSuperClass):
    """
    Class for collecting and calcuating 1D fields
    """

    #{{{constructor
    def __init__(self                   ,\
                 *args                  ,\
                 mode       = "parallel",\
                 processing = None      ,\
                 return2d   = True      ,\
                 **kwargs):
        #{{{docstring
        """
        This constructor will:
            * Call the parent constructor
            * Set member data

        Parameters
        ----------
        *args : positional arguments
            See parent constructor for details.
        mode : ["radial"|"parallel"]
            * "radial"   - Radial profiles will be used
            * "parallel" - Parallel profiles will be used
            * "poloidal" - Poloidal profiles will be used
        processing : [None|str]
            * None                 - Raw data will be used
            * "polAvg"             - <var>_theta will be used
            * "polAndTimeAvg"      - <<var>_theta>_t will be used
            * "polAvgFluct"        - var - <var>_theta will be used
            * "polAndTimeAvgFluct" - var - <<var>_theta>_t will be used
        return2d : bool
            Whether or not the variable should be returned as a
            2d-array (if True) or 4d-array (if False).
        *kwargs : keyword arguments
            See parent constructor for details.
        """
        #}}}

        # Guard
        implemented = ("parallel", "radial", "poloidal")
        if not(mode in implemented):
            message = "mode '{}' not implemented".format(mode)
            raise NotImplementedError(message)

        if processing is not None:
            implemented = (\
                           "polAvg",\
                           "polAndTimeAvg",\
                           "polAvgFluct",\
                           "polAndTimeAvgFluct",\
                          )
            if not(processing in implemented):
                message = "processing '{}' not implemented".format(processing)

        # Call the constructor of the parent class
        super().__init__(*args, **kwargs)

        self._mode       = mode
        self._processing = processing
        self._return2d   = return2d
    #}}}

    #{{{executeCollectAndCalc
    def executeCollectAndCalc(self):
        #{{{docstring
        """
        Function which collects and calculates 1D fields.

        Returns
        -------
        field1D : dict
            Dictionary with the keys:
                * var    - A 2d array (a 1d spatial array of each time)
                           of the collected variable if return2d is True.
                           If return2d is false, a 4d variable is
                           returned instead.
                * "X"    - The abscissa of the variable
                * "time" - The time trace
                * pos1   - The position of the first fixed index
                * pos2   - The position of the second fixed index
        """
        #}}}

        # Guard
        if len(self._notCalled) > 0:
            message = "The following functions were not called:\n{}".\
                        format("\n".join(self._notCalled))
            raise RuntimeError(message)

        # Initialize output
        field1D = {}

        # Obtain the mesh
        if self._mode == "parallel":
            X = self._dh.z
        elif self._mode == "radial":
            X = self._dh.rho
        elif self._mode == "poloidal":
            X = self._dh.thetaRad

        # Convert to indices
        xInd = slicesToIndices(self._collectPaths[0], self._xSlice, "x",\
                               xguards=self._xguards)
        yInd = slicesToIndices(self._collectPaths[0], self._ySlice, "y",\
                               yguards=self._yguards)
        zInd = slicesToIndices(self._collectPaths[0], self._zSlice, "z")
        tInd = slicesToIndices(self._collectPaths[0], self._tSlice, "t")

        collectGhost = True if (self._xguards or self._yguards) else False

        # Set keyword arguments
        collectKwargs = {\
            "collectGhost" : collectGhost,\
            "tInd"         : tInd        ,\
            "xInd"         : xInd        ,\
            "yInd"         : yInd        ,\
            "zInd"         : zInd        ,\
                }

        # Collect
        var, time = self._collectWrapper(collectKwargs)

        if self.convertToPhysical:
            var  = self.uc.physicalConversion(var , self._varName)
            time = self.uc.physicalConversion(time, "t")

        # Store the fields
        field1D["X"   ] = X
        field1D["time"] = time

        if "parallel" in self._mode:
            if self._return2d:
                field1D[self._varName] = var[:, 0, :, 0]
            else:
                field1D[self._varName] = var
            field1D["rhoPos"]   = self._dh.rho     [xInd[0]]
            field1D["thetaPos"] = self._dh.thetaDeg[zInd[0]]
        elif "radial" in self._mode:
            if self._return2d:
                field1D[self._varName] = var[:, :, 0, 0]
            else:
                field1D[self._varName] = var
            field1D["zPos"]     = self._dh.z       [yInd[0]]
            field1D["thetaPos"] = self._dh.thetaDeg[zInd[0]]
        elif "poloidal" in self._mode:
            if self._return2d:
                field1D[self._varName] = var[:, 0, 0, :]
            else:
                field1D[self._varName] = var
            field1D["rhoPos"]   = self._dh.rho[xInd[0]]
            field1D["zPos"]     = self._dh.z  [yInd[0]]

        return field1D
    #}}}

    #{{{_collectWrapper
    def _collectWrapper(self, collectKwargs):
        #{{{docstring
        """
        Collects the variable and the time.

        Parameters
        ----------
        collectKwargs : dict
           Keyword arguments to use in the collect

        Returns
        -------
        var : array-4d
            The collected array.
        time : array
            The time array.
        """
        #}}}

        # Prepare the collectKwargs
        if self._mode == "parallel":
            collectKwargs.pop("yInd")
        elif self._mode == "radial":
            collectKwargs.pop("xInd")
        elif self._mode == "poloidal":
            collectKwargs.pop("zInd")

        if self._processing:
            if "pol" in self._processing:
                collectKwargs.pop("zInd")

        if not(self._return2d):
             if "zInd" in collectKwargs.keys():
                collectKwargs.pop("zInd")

        # Decide if theta should be collected
        if not(self._mode == "poloidal"):
            if (self._processing is not None) and\
               (self._processing == "poloidal"):
                collectTheta = True
            else:
                collectTheta = False

        # Get collect function
        if self._mode == "parallel":
            if not(collectTheta) and self._return2d:
                collecter = collectParallelProfile
            else:
                collecter = collectConstRho
        elif self._mode == "radial":
            if not(collectTheta) and self._return2d:
                collecter = collectRadialProfile
            else:
                collecter = collectConstZ
        elif self._mode == "poloidal":
            collecter = collectPoloidalProfile

        # Collect
        var = collecter(self._collectPaths, self._varName, **collectKwargs)

        time = collectTime(self._collectPaths, collectKwargs["tInd"])

        # Process
        if self._processing is not None:
            if "pol" in self._processing:
                polAvgVar = polAvg(var)
                if "time" in self._processing.lower():
                    polAvgTimeAvgVar, timeAvgTime = timeAvg(polAvgVar, t=time)
            if self._processing == "polAvg":
                var = polAvg
            elif self._processing == "polAndTimeAvg":
                var  = polAvgTimeAvgVar
                time = timeAvgTime
            elif self._processing == "polAvgFluct":
                var = var - polAvgVar
            elif self._processing == "polAndTimeAvgFluct":
                var  = var - polAvgTimeAvgVar
                time = timeAvgTime

        # Slice
        if self._tSlice is not None:
            if type(self._tSlice) == slice:
                if self._tSlice.step is not None:
                    var  = var [::self._tSlice.step]
                    time = time[::self._tSlice.step]

        return var, time
    #}}}

    #{{{getDh
    def getDh(self):
        """
        Returns the dimension helper.

        Returns
        -------
        dh : DimensionHelper
            The dimension helper object
        """
        return self._dh
    #}}}
#}}}

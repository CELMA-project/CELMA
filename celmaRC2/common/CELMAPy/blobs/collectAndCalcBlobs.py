#!/usr/bin/env python

"""
Contains the blobs calculation
"""

from ..calcVelocities import calcRadialExBPoloidal


#{{{CollectAndCalcBlobs
class CollectAndCalcBlobs(object):
    """
    Class for collecting and calcuating the blobs
    """

    #{{{constructor
    def __init__(self             ,\
                 collectPaths     ,\
                 slices           ,\
                 convertToPhysical,\
                 ):
        #{{{docstring
        """
        This constructor will:
            * Set the member data

        Parameters
        ----------
        collectPaths : tuple
            Tuple from where to collect
        slices : tuple of tuples
            Tuple the indices to use.
            On the form (xInd, yInd, zInd, tSlice)
        convertToPhysical : bool
            Whether or not to convert to physical
        """
        #}}}

        # Set the member data
        self._collectPaths      = collectPaths
        self._convertToPhysical = convertToPhysical
        self._xInd, self._yInd, self._zInd, self._tSlice = slices
    #}}}

    #{{{executeCollectAndCalc
    def executeCollectAndCalc(self):
# FIXME: Clean this
# FIXME: Consider to have getters for the output
        #{{{docstring
        """
        This function will:
            * Collect the density flux at the probe place.
            * Find the times where the fluxs meets the condition.
            * Collect the density at the time where the condition is met.
            * Sort the blobs from the holes (negative perturbation
              transported inwards).
            * Make the averages.

        NOTE: The treshold condition is made in the flux as:
                1. One can have high perturbation as a consequence of
                   high azimuthal flux.
                2. The plasma can be alongated and spin in the azimuthal
                   direction.

        Organization of ouput dict
        --------------------------
        The output dicts contains the keys
            * par2D     - The parallel 2D slice
            * perp2D    - The perpednicular 2D slice
            * timeTrace - The time trace in the point
        All the keys will have one of the following suffixes
            * Blobs     - A tuple of the blobs bins dictionary
            * BlobsAvg  - A dictionary where the bins have been averaged
            * Holes     - A tuple of the holes bins dictionary
            * HolesAvg  - A dictionary where the bins have been averaged
        The sub dictionaries will contain the following keys
            * "n"
            * "time"
            * ""

        Returns
        -------
        blobBinsDict : dict
            Dict of the blob-bins dicts.
            See "Organization of ouput dict" for more info.
        holeBinsDict : dict
            Dict of the holes-bins dicts.
            See "Organization of ouput dict" for more info.
        blobsAvgDict : dict
            Dict of the average blobs dicts.
            See "Organization of ouput dict" for more info.
        holesAvgDict : dict
            Dict of the average holes dicts.
            See "Organization of ouput dict" for more info.
        """
        #}}}

        # Collect flux
        radialFlux, self._uc = self._collectRadialFlux()
        self._dh = DimensionHelper(self._collectPaths[0], self._uc)

        # Initialize
        key = radialFlux.keys()[0]
        flux = radialFlux[key]["nRadialFlux"]
        condition = flux.std()*3
        rhoPos, thetaPos, zPos = key.split(",")
        time = radialFlux[key]["time"]
        dt = time[1] - time[0]

        # Get indices and window sizes
        indices = self._getIndicesMeetingCondition(flux, condition)
        waitingTimes, pulseWidths =\
                self._getWaitingTimesAndPulseWidth(indices, dt)
        windowSize = self._getWindowSize(indices, self._pctPadding)
        slices =\
            self._transformContiguousIndicesToSlices(contiguousIndices,\
                                                     windowSize)

        # Collect the bins
        par2DBins, perp2DBins = self._collect2DBins(slices)
        timeTraceBins         = self._getTimeTraceBins(perp2DBins)

        # Sort into blobs and holes
        midIndex = int((slices[0].stop - slices[0].start)/2)
        blobsIndices, holesIndices =\
            self._idendtifyBlobsAndHoles(self, timeTraceBins, midIndex)
        par2DBlobs, par2DHoles =\
            self._extractBlobsAndHoles(blobsIndices, holesIndices, par2DBins)
        perp2DBlobs, perp2DHoles =\
            self._extractBlobsAndHoles(blobsIndices, holesIndices, perp2DBins)
        timeTraceBlobs, timeTraceHoles =\
            self._extractBlobsAndHoles(blobsIndices, holesIndices, timeTraceBins)

        # Make averages
        # +1 for symmetry
        self._timeSlice = np.array(range(-windowSize,windowSize+1))*dt
        blobs = (par2DBlobs, perp2DBlobs, timeTraceBlobs)
        holes = (par2DHoles, perp2DHoles, timeTraceHoles)

        par2DBlobAvg, perp2DBlobAvg, timeTraceBlobAvg =\
            self._calcAverages(blobs)

        par2DBHolesAvg, perp2DHolesAvg, timeTraceHolesAvg =\
            self._calcAverages(blobs)

        # Make return dicts
        blobBinsDict = {\
            "par2DBlobs"     : par2DBlobs    ,\
            "perp2DBlobs"    : perp2DBlobs   ,\
            "timeTraceBlobs" : timeTraceBlobs,\
                   }

        holeBinsDict = {\
            "par2DHoles"     : par2DHoles    ,\
            "perp2DHoles"    : perp2DHoles   ,\
            "timeTraceHoles" : timeTraceHoles,\
                       }

        blobsAvgDict = {\
            "par2DBlobAvg"     : par2DBlobAvg    ,\
            "perp2DBlobAvg"    : perp2DBlobAvg   ,\
            "timeTraceBlobAvg" : timeTraceBlobAvg,\
                      }

        holesAvgDict = {\
            "par2DBHolesAvg"    : par2DBHolesAvg   ,\
            "perp2DHolesAvg"    : perp2DHolesAvg   ,\
            "timeTraceHolesAvg" : timeTraceHolesAvg,\
                   }

        return blobBinsDict, holeBinsDict, blobsAvgDict, holesAvgDict
    #}}}

    #{{{_collectRadialFlux
    def _collectRadialFlux(self):
        #{{{docstring
        """
        Collects the radial flux in the probe point

        Returns
        -------
        radialFlux : dict
            Dictionary where the keys are on the form "rho,theta,z".
            The value is a dict containing of
            {varName:radialFlux, "time":time, "timeIndices":timeIndices}
        uc : UnitsConverter
            The UnitsConverter object.
        """
        #}}}
        varName = "n"
        mode    = "fluct"

        indicesArgs   = (self._xInd, self._yInd, self._zInd)
        indicesKwargs = {"tSlice"          : self._tSlice,\
                         "nPoints"         : 1           ,\
                         }

        radialFlux, uc = getRadialFlux(self._collectPaths     ,\
                                       varName                ,\
                                       self._convertToPhysical,\
                                       mode                   ,\
                                       indicesArgs            ,\
                                       indicesKwargs          ,\
                                      )

        return radialFlux, uc
    #}}}

    #{{{_getIndicesMeetingCondition
    def _getIndicesMeetingCondition(self, var, condition):
        #{{{docstring
        """
        Returns indices where the condition is meet.

        The indices are corrected for the tSlice.start.

        Parameters
        ----------
        var : array-1d
            Array of to find where the condition is meet.
        condition : float
            The condition to check for.

        Returns
        -------
        contiguousIndices : tuple
            Tuple of arrays, where each element in the tuple is a series
            of indices, contiguous, in ascending order where the
            condition is meet.
        """
        #}}}

        # Get indices where the values are above the conditions
        indicesMeetingCond = np.where(var >= condition)[0]

        # Get contiguous indices
        contiguousIndices = []
        while len(indicesMeetingCond) !=0:
            # -1 in order not to get out of index on last point
            for i in range(len(indicesMeetingCond)-1):
                if indicesMeetingCond[i]+1 != indicesMeetingCond[i+1]:
                    break
            contiguous = indicesMeetingCond[0:i+1]
            indicesMeetingCond = indicesMeetingCond[i+1:]
            contiguousIndices.append(contiguous)

        # As i ran to len(indicesMeetingCond)-1, we make a special treatment on
        # the last index
        contiguousIndices[-2] =\
            np.append(contiguousIndices[-2], contiguousIndices[-1][0])
        contiguousIndices = contiguousIndices[:-1]

        # Correct for tSlices.start
        for i in range(len(contiguousIndices)):
            contiguousIndices[i] += self._tSlice.start

        contiguousIndices = tuple(contiguousIndices)

        return contiguousIndices
    #}}}

    #{{{_getWaitingTimesAndPulseWidth
    def _getWaitingTimesAndPulseWidth(self, contiguousIndices, dt):
        #{{{docstring
        """
        Returns the waiting times and the pulse widths

        Parameters
        ----------
        contiguousIndices : tuple
            Tuple of arrays, where each element in the tuple is a series
            of indices, contiguous, in ascending order where the
            condition is meet.
        dt : float
            The time incrementation between two indices.

        Returns
        -------
        waitingTimes : tuple
            Tuple of waiting times in time units.
        waitingTimes : tuple
            Tuple of pulse widths in time units.
        """
        #}}}

        contiguousTimes = tuple(indices*dt for indices in contiguousindices)
        waitingTimes = []
        pulseWidths  = []
        for times in contiguousTimes:
            pulseWidths.append(times[-1]-times[0])
            mid = int(len(times)/2)
            if len(pulseWidths) == 1:
                # Initialize waiting times
                waitingStart = times[mid]
                continue

            waitingEnd = times[mid]
            waitingTimes.append(waitingEnd - waitingStart)
            waitingStart = waitingEnd.copy()

        waitingTimes = tuple(waitingTimes)
        pulseWidths  = tuple(pulseWidths)

        return waitingTimes, pulseWidths
    #}}}

    #{{{_getWindowSize
    def _getWindowSize(self, contiguousIndices, pctPadding):
        #{{{docstring
        """
        Get the window size which will define the size of one bin.

        Parameters
        ----------
        contiguousIndices : tuple
            Tuple of arrays, where each element in the tuple is a series
            of indices, contiguous, in ascending order where the
            condition is meet.
        pctPadding : float
            Padding (in percent) which will be added to the max length
            of the contiguousIndices

        Returns
        -------
        windowSize : int
            Number of indices to defining the size of one bin.
        """
        #}}}

        maxLen = 0
        for indices in contiguousIndices:
            maxLen = len(indices) if len(indices) > maxLen else maxLen
        windowSize = int(maxLen*(1+(pctPadding/100)))

        return windowSize
    #}}}

    #{{{_transformContiguousIndicesToSlices
    def _transformContiguousIndicesToSlices(self, contiguousIndices, windowSize):
        #{{{docstring
        """

        Parameters
        ----------
        contiguousIndices : tuple
            Tuple of arrays, where each element in the tuple is a series
            of indices, contiguous, in ascending order where the
            condition is meet.
        windowSize : int

        Returns
        -------
        slices : tuple
            Tuple of slices, where the individual slices are the slices
            which will be used to collect a bin.
        windowSize : int
            Number of indices to defining the size of one bin.
        """
        #}}}

        # Transform to slices
        slices = []
        for indices in contiguousIndices:
            # Find the mid of the indices
            mid = int(len(indices)/2)
            # +1 for symmetry
            curSlice = slice(indices[mid]-windowSize, indices[mid]+windowSize+1)
            # Guard for the beginning
            if curSlice.start >= 0:
                # Guard for the end
                if curSlice.stop <= len(dens):
                    slices.append(curSlice)
        slices = tuple(slices)

        return slices
    #}}}

    #{{{_collect2DBins
    def _collect2DBins(self, slices):
        #{{{docstring
        """
        Collects the bins which will be used in the average.

        Parameters
        ----------
        slices : tuple
            Tuple of slices, where the individual slices are the slices
            which will be used to collect a bin.

        Returns
        -------
        par2DBins : tuple
            A tuple of the parallel 2D bins stored as dicts with the
            keys:
                * "n"    - A 3d array (a 2d spatial array of each time)
                           of the collected variable.
                * "nPPi" - The field at pi away from the varName field
                * "X"    - The cartesian x mesh to the field
                * "Y"    - The cartesian Y mesh to the field
                * "time" - The time trace
                * pos    - The position of the fixed index
        perp2DBins : tuple
            A tuple of the perpendicular 2D bins stored as dicts with
            the keys:
                * "n"    - A 3d array (a 2d spatial array of each time)
                           of the collected variable.
                * "X"    - The cartesian x mesh to the field
                * "Y"    - The cartesian Y mesh to the field
                * "time" - The time trace
                * pos    - The position of the fixed index
        """
        #}}}

        fluct             = True
        xSlice            = None
        ySlice            = self._ySlice
        zSlice            = self._zSlice
        convertToPhysical = self._convertToPhysical
        varName           = "n"

        par2DBins  = []
        perp2DBins = []

        for tSlice in slices:
            # Pependicular collection
            ccf2D = CollectAndCalcFields2D(\
                        collectPaths              ,\
                        fluct             = fluct ,\
                        mode              = "perp",\
                        convertToPhysical = convertToPhysical)

            ccf2D.setSlice(xSlice, yInd, zSlice, tSlice)
            ccf2D.setVarName(varName)
            perp2DBins.append(ccf2D.executeCollectAndCalc())

            # Parallel collection
            ccf2D = CollectAndCalcFields2D(\
                        collectPaths             ,\
                        fluct             = fluct,\
                        mode              = "par",\
                        convertToPhysical = convertToPhysical)
            ccf2D.setSlice(xSlice, ySlice, zInd, tSlice)
            ccf2D.setVarName(varName)
            par2DBins.append(ccf2D.executeCollectAndCalc())

        return par2DBins, perp2DBins
    #}}}

    #{{{_getTimeTraceBins
    def _getTimeTraceBins(self, perp2DBins):
        #{{{docstring
        """
        Get the time trace bins from perp2DBins

        Parameters
        ----------
        perp2DBins : tuple
            A tuple of the perpendicular 2D bins stored as dicts with
            the keys:
                * "n"    - A 3d array (a 2d spatial array of each time)
                           of the collected variable.
                * "X"    - The cartesian x mesh to the field
                * "Y"    - The cartesian Y mesh to the field
                * "time" - The time trace
                * pos    - The position of the fixed index

        Returns
        -------
        timeTraceBins : tuple
            A tuple of the time traces, stored as dicts with the keys:
                * "n"    - A 3d array (a 2d spatial array of each time)
                           of the collected variable.
                * "time" - The time trace
                * "pos"  - Tuple of the (rho, theta, z) fixed positions.
        """
        #}}}

        timeTraceBins = []
        rhoPos   = self._dh.rho     [self._xInd]
        thetaPos = self._dh.thetaRad[self._yInd]
        for perp2DBin in perp2DBins:
            curDict         = {}
            curDict["n"]    = perp2DBins["n"][:, self._xInd, self._yInd, self._zInd]
            curDict["time"] = perp2DBins["time"]
            curDict["pos"]  = (rhoPos, thetaPos, perp2DBins["z"])
            timeTraceBins.append(curDict)

        return timeTraceBins
    #}}}

    #{{{_idendtifyBlobsAndHoles
    def _idendtifyBlobsAndHoles(self, timeTraceBins, midIndex):
        #{{{docstring
        """
        Identifies whether the bin contains a blob or a hole.

        An outward transportation of positive fluctuations is a positive
        flux.
        A negative transportation of a negative fluctuation is also a
        positive flux.

        Parameters
        ----------
        timeTraceBins : tuple
            A tuple of the time traces, stored as dicts with the keys:
                * "n"    - A 3d array (a 2d spatial array of each time)
                           of the collected variable.
                * "time" - The time trace
                * "pos"  - Tuple of the (rho, theta, z) fixed positions.
        midIndex : int
            Index of the mid of the time trace.

        Returns
        -------
        blobsIndices : tuple
            Tuple containing the indices identified as blobs
        holesIndices : tuple
            Tuple containing the indices identified as holes
        """
        #}}}

        blobsIndices  = []
        holesIndices = []
        for nr, timeTrace in enumerate(timeTraceBins):
            curDens = timeTrace["n"]
            if curDens[midIndex] >= 0:
                blobsIndices.append(nr)
            else:
                holesIndices.append(nr)

        blobsIndices = tuple(blobsIndices)
        holesIndices = tuple(holesIndices)

        return blobsIndices, holesIndices
    #}}}

    #{{{_extractBlobsAndHoles
    def _extractBlobsAndHoles(self, blobsIndices, holesIndices, bins):
        #{{{docstring
        """
        Extracts blobs and holes from a tuple of bins.

        Parameters
        ----------
        blobsIndices : tuple
            Tuple of blob indices.
        holesIndices : tuple
            Tuple of holes indices.
        bins : tuple
            Bins to be separated.

        Returns
        -------
        blobs : tuple
            Tuple of the bins which are blobs.
        holes : tuple
            Tuple of the bins which are holes.
        """
        #}}}

        blobs = []
        holes = []

        for i in blobsIndices:
            blobs.append(bins[i])
        for i in holesIndices:
            holes.append(bins[i])

        return blobs, holes
    #}}}

    #{{{_calcAverages
    def _calcAverages(self, tupleOfDicts):
        #{{{docstring
        """
        Returns a tuple of dict where "n" and "time" in the dicts has
        been replaced with averaged values.

        Parameters
        ----------
        tupleOfDicts : tuple
            Tuple of dicts, where each dict contains the keys "n" and
            "time".

        Returns
        -------
        averages : tuple
            Tuple of the dicts where "n" and "time" has been averaged.
        """
        #}}}
        averages = []

        for curDict in tupleOfDicts:
            avgDict         = curDict.copy()
            avgDict["n"]    = np.array(curDict["n"]).mean(axis=0)
            avgDict["time"] = self._timeSlice
            averages.append(avgDict)

        return averages
    #}}}
#}}}

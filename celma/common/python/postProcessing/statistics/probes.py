#!/usr/bin/env python

"""
Contains classes which probes the data
"""

import numpy as np
from scipy.stats import kurtosis, skew
from scipy.signal import periodogram
from boutdata import collect
from .polAvg import polAvg
from .collectiveCollect import collectiveCollect
from .derivatives import DDZ, findLargestRadialGrad

#{{{class Probes
class Probes(object):
    """
    Parent class which stores the Jacobian J, the variable and the
    time traces.

    Contains analysis functions.
    """

    #{{{Constructor
    def __init__(self, var, varName, time,\
                 steadyStatePath=None, radialProbeIndices=None,\
                 collectPath=None):
        """
        Constructor for the Probes class

        Parameters
        ----------
        var : array
            Variable that will be stored as a member data
        varName : string
            Name of the variable.
        time : array
            The time array
        steadyStatePath : str
            What path to use when collecting J. If radialProbeIndices is
            None, this will also be the path for finding the largest
            gradient. Default is None.
        radialProbeIndices : [None|array]
            What radial indices to probe. If set to None, it will be
            selected from the largest gradient in the steady state
            variable.
        collectPath : str
            Path to collect J from. Not effective if steadyStatePath is set.
        """

        # Guard
        if steadyStatePath is None and radialProbeIndices is None:
            message = ("steadyStatePath and radialProbeIndices cannot "
                       "be None simultaneously")
            raise ValueError(message)
        if steadyStatePath is None and collectPath is None:
            message = ("steadyStatePath and collectPath cannot "
                       "be None simultaneously")
            raise ValueError(message)

        self._var  = var
        self._time = time

        # Find the fluctuations in var
        self._varFluct = var - polAvg(var)

        self._varName = varName
        if varName == "n":
            varName = "lnN"

        if collectPath is None:
            collectPath = steadyStatePath

        # Set the Jacobian
        # Contains the ghost points as we are using this in DDZ
        self._J = collect("J", path=collectPath,\
                          xguards=True, yguards=True, info=False)

        if radialProbeIndices == None:
            # Note that the ghost cells are collected, as we are taking
            # derivatives of the field
            self._varSteadyState =\
                collect(varName, path=collectPath,\
                        xguards=True, yguards=True, info=False)

            if self._varName == "n":
                self._varSteadyState = np.exp(self._varSteadyState)

        # Set uninitialized variables to None
        self.results              = None
        self._timeTraceOfVar      = None
        self._timeTraceOfVarFluct = None
        self._xInds               = None
        self._yInd                = None
        self._yInds               = None
        self._zInds               = None
    #}}}

    #{{{initializeInputOutput
    def initializeInputOutput(self, xInds, yInds, zInds):
        """
        Get the time trace and the fluctuation time trace.
        Initializes the results dictionary.
        Sets the member indices

        Parameters
        ----------
        xInds : array
            x indices for the time trace
        yInds : [int|array]
            y index or indices for the time trace
        zInds : array
            z indices for the time trace
        """

        self._timeTraceOfVar      = {}
        self._timeTraceOfVarFluct = {}

        # Set member data
        if self._yInd is not None:
            self._yInds = yInds
        else:
            # If self._yInds has been set in PerpPlaneProbes
            self._yInds = [0]

        self._xInds = xInds
        self._zInds = zInds

        for xInd in self._xInds:
            for zInd in self._zInds:
                self._timeTraceOfVar\
                        ["{}{}{}".format(xInd, self._yInd, zInd)] =\
                            self._var[:, xInd, self._yInd, zInd]
                self._timeTraceOfVarFluct\
                        ["{}{}{}".format(xInd, self._yInd, zInd)] =\
                            self._varFluct[:, xInd, self._yInd, zInd]

                # Initialize the result as a dictionary
                self.results["{}{}{}".format(xInd, self._yInd, zInd)] = {}
    #}}}

    #{{{calcStatsMoments
    def calcStatsMoments(self):
        """
        Calculates the first statistical moments.

        Output
        ------
        The output will be stored in the keys (specified below) under
        self._result[indexString], where indexString is the string of
        the index under consideration.

        mean : array
            The mean of the fluctuations for each time (should be 0).
        var : array
            The variance of the fluctuations for each time.
        kurtosis : array
            The kurtosis of the fluctuations for each time (is 3.0 for a
            normal distribution).
        skew : array
            The skewness of the the fluctuations for each time (is 0 for
            a normal distribution)
        """

        for xInd in self._xInds:
            for yInd in self._yInds:
                for zInd in self._zInds:
                    key = "{}{}{}".format(xInd, yInd, zInd)
                    self.results[key]["mean"] =\
                        self._timeTraceOfVarFluct[key].mean()
                    self.results[key]["var"] =\
                        self._timeTraceOfVarFluct[key].var()
                    self.results[key]["kurtosis"] =\
                        kurtosis(self._timeTraceOfVarFluct[key])
                    self.results[key]["skew"] =\
                        skew(self._timeTraceOfVarFluct[key])
    #}}}

    #{{{calcPDFs
    def calcPDFs(self):
        """
        Calculates the probability distribution function (PDF) of the
        fluctuations.

        Probability distribution function
        ---------------------------------
        Probability that the measurement falls within an infinite small interval.

        Output
        ------
        The output will be stored in the keys (specified below) under
        self._result[indexString], where indexString is the string of
        the index under consideration.

        pdfX : float
            The bins of the PDF.
        pdfY : float
            The counts of the PDF
        """

        # Guard
        if self.results is None:
            message = ("Results not initiazled. Please do so by calling "
                       "initializeInputOutput")
            raise RuntimeError(message)

        # n is a measure of count per bin
        # Histogram counts the occurences of values within a specific interval
        # Density normalizes so that the integral (of the continuous variable)
        # equals one, note that the sum of histograms is not necessarily 1)
        # http://docs.scipy.org/doc/numpy/reference/generated/numpy.histogram.html
        # http://stackoverflow.com/questions/36150257/probability-distribution-function-python/36248810
        for xInd in self._xInds:
            for yInd in self._yInds:
                for zInd in self._zInds:
                    key = "{}{}{}".format(xInd, yInd, zInd)

                    self.results[key]["pdfX"], bins =\
                              np.histogram(self._timeTraceOfVarFluct,\
                                           bins="auto",\
                                           density=True)
                    # Initialize y
                    self.results[key]["pdfY"] =\
                            np.zeros(self.results[key]["pdfX"].size)

                    for k in range(self.results[key]["pdfX"].size):
                        # Only the bin edges are saved. Interpolate to the center of the bin
                        self.results[key]["pdfY"] = 0.5*(bins[k]+bins[k+1])
    #}}}

    #{{{calcPSDs
    def calcPSDs(self):
        """
        Estimates the power spectral density (PSD) of the fluctuations
        in a non-parametric way (no assumption for the model of the
        random process). Using a periodogram estimate.

        Power spectral density
        ----------------------
            * Tells us what frequencies are present in the signal.
            * The average power of the signal is the integrated of the PSD
            * The bandwidth of the process (turbulence) is defined as
              the frequency width where the signal is within 3dB of its
              peak value.
            * PSD is a deterministic description of the spectral
              characteristic of the signal. Cannot use fourier transform on
              random variables as this is not necessarily definied etc.
            * PSD is the fourier transformed of the auto-correlation
              function, and cross spectral density is the fourier
              transformed of the cross correlation function
            * The periodogram estimate is the same as autocorrelation with a
              triangular window

        Output
        ------
        The output will be stored in the keys (specified below) under
        self._result[indexString], where indexString is the string of
        the index under consideration.

        psdX : float
            The frequency of the PSD
        psdY : float
            The power spectral density measured in variableUnits**2/Hz
        """

        # Guard
        if self.results is None:
            message = ("Results not initiazled. Please do so by calling "
                       "initializeInputOutput")
            raise RuntimeError(message)

        fs = self._time[1] - self._time[0]
        for xInd in self._xInds:
            for yInd in self._yInds:
                for zInd in self._zInds:
                    key = "{}{}{}".format(xInd, yInd, zInd)

                    # window = None => window = "boxcar"
                    # scaling = density gives the correct units
                    self.results[key]["psdX"], self.results[key]["psdY"] =\
                        periodogram(self._timeTraceOfVarFluct[key],\
                                    fs=fs, window=None, scaling="density")
    #}}}

    # TODO: FIXME: Resolve the issue of fluctuations
    #{{{calcFluxThroughVolumeElement
    def calcFluxThroughVolumeElement(self, u, uName):
        """
        Gives the radial flux through a volume element.

        Note that this is different from the total flux through a surface,
        which is the surface integral of the same quantity.
        However, this is usually the flux through a voluem element which are
        measured with probes in experiments.

        The flux of the variable and the flux of the fluctuating quantity is
        returned.

        FIXME: Resolve two next paragraphs
        A note on poloidal averaging:
        ab = (a_avg + a_fluct)(b_avg + b_fluct)
           =  a_avg b_avg + a_avg b_fluct + a_fluct b_ avg + a_fluct b_fluct
        Usually we are looking at time averaged quantities of this
        (ab)_avg = (a_avg b_avg)_avg + (a_avg b_fluct)_avg + (a_fluct b_ avg)_avg + (a_fluct b_fluct)_avg
        (a_avg b_avg)_avg = a_avg b_avg
        (a_fluct)_avg = 0 if a is truly random, and a sufficient average is taken
        c(a_fluct)_avg = 0
        b_avg = constant => (a_avg b_fluct)_avg = 0

        We are doing a poloidal average
        Assume no preferred direction => poloidal average is comparable to a
        time average long enough for the averages to be the same (turbulence
        have "walked" around poloidally)


        Parameters
        ----------
        var : array
            The variable to find the flux of (the xInd and yInd must be
            specified)
        u : array
            The velocity of the flux in the direction of the flux (xInd and yInd must
            be specified)
        uName : str
            Name of the velocity

        Output
        ------
        The output will be stored in the keys (specified below) under
        self._result[indexString], where indexString is the string of
        the index under consideration.

        varFlux* : array
            The flux at the (xInd, yInd, zInd) position for each time.
        fluctVarFlux* : array
            The flux arising from the fluctuation at the (xInd, yInd, zInd)
            position for each time.

        * The keys below will be appended with uName.
        """

        # Guard
        if self.results is None:
            message = ("Results not initiazled. Please do so by calling "
                       "initializeInputOutput")
            raise RuntimeError(message)

        # Find the fluctuating velocity
        uFluct = u - polAvg(u)


        for xInd in self._xInds:
            for yInd in self._yInds:
                for zInd in self._zInds:
                    key = "{}{}{}".format(xInd, yInd, zInd)
                    self.results[key]["varFlux" + uName.capitalize()] =\
                        self._timeTraceOfVar[key] * u[:, xInd, yInd, zInd]
                    self.results[key]["varFluxFluct" + uName.capitalize()] =\
                        self._timeTraceOfVarFluct[key] * uFluct[:, xInd, yInd, zInd]
    #}}}

    #{{{calcFFTs
    def calcFFTs(self):
        """
        Calculates the FFT of the poloidal profile belonging to the
        point under consideration.

        Output
        ------
        The output will be stored in the keys (specified below) under
        self._result[indexString], where indexString is the string of
        the index under consideration.

        zFFT : array
            The fourier transformed of the z-direction for each time. Notice that the
            result will be the same for every z-index for a fixed x- and
            y-index.
        """

        # Guard
        if self.results is None:
            message = ("Results not initiazled. Please do so by calling "
                       "initializeInputOutput")
            raise RuntimeError(message)

        for xInd in self._xInds:
            for yInd in self._yInds:

                varFluctFFT =\
                    np.abs(np.fft.fft(self._varFluct[xInd][yInd], axis=-1))
                for zInd in self._zInds:
                    key = "{}{}{}".format(xInd, yInd, zInd)
                    self.results[key]["zFFT"] = varFluctFFT
    #}}}
#}}}

#{{{class PerpPlaneProbes
class PerpPlaneProbes(Probes):
    """
    Child class of probes. Calls probes constructor, finds
    radialProbeIndice based on highest gradient (if not set).
    Also calculates the ExB velocity.

    This class makes it possible to collect only a z-plane instead of
    the whole 3D variable.
    """

    #{{{Constructor
    def __init__(self, varName, paths, yInd,\
                 nProbes=5, physicalUnits=False,\
                 steadyStatePath=None, radialProbeIndices=None):
        """
        Constructor for the PerpPlaneProbes class

        Note that ghost points are collected, as we are finding
        derivatives of the fields.

        Parameters
        ----------
        varName : string
            Name of the variable.
        paths : string
            What path to use when collecting the variable. Must be in
            ascending temporal order as the variable will be
            concatenated.
        yInd : int
            yInd to collect from
        nProbes : int
            Number of probes. Default is 5. Will be overridden by
            radialProbeIndices if set.
        physicalUnits : bool
            Will normalize the time array and collect the normalization
            constants if True. Note that normalization of the other
            variables must be recalculated manually in order to obtain
            the quantities in physical units. Default is False.
        steadyStatePath : string
            What path to use when collecting J. If radialProbeIndices is
            None, this will also be the path for finding the largest
            gradient. Default is None.
        radialProbeIndices : [None|array]
            What radial indices to probe. If set to None, it will be
            selected from the largest gradient in the steady state
            variable.
        """

        # Guard
        if steadyStatePath is None and radialProbeIndices is None:
            message = ("steadyStatePath and radialProbeIndices cannot "
                       "be None simultaneously")
            raise ValueError(message)

        if varName == "n":
            varName = "lnN"

        # Collect the variable
        varTimeDict = collectiveCollect(paths, [varName, "t_array"],\
                                        collectGhost = False,\
                                        yInd = [yInd, yInd],\
                                        )

        var  = varTimeDict[varName]
        time = varTimeDict["t_array"]

        if varName == "n":
            var = np.exp(var)

        # Collect and recalculate the time
        if physicalUnits:
            self.n0   = collect("n0" ,  path = paths[-1])
            self.Te0  = collect("Te0",  path = paths[-1])
            self.Ti0  = collect("Ti0",  path = paths[-1])
            self.B0   = collect("B0" ,  path = paths[-1])
            self.omCI = collect("omCI", path = paths[-1])
            self.rhoS = collect("rhoS", path = paths[-1])

            time /= self.omCI

        # Call the parent class
        super().__init__(var,\
                         varName,\
                         time,\
                         steadyStatePath,\
                         radialProbeIndices)

        self._yInd = yInd

        if radialProbeIndices is None:
            self._MXG = collect("MXG", path = paths[-1])
            dx = collect("dx", path = paths[-1], info = False)
            # Find the max gradient of the variable
            _, maxGradInd =\
                findLargestRadialGrad(\
                  self._varSteadyState[0:1, :, self._yInd:self._yInd+1, 0:1],\
                  dx,\
                  self._MXG)
            self.radialProbesIndices =\
                self.getRadialProbeIndices(maxGradInd, nProbes)
        else:
            self.radialProbesIndices = radialProbeIndices

        # Get the radial ExB drift for this plane
        phi = collectiveCollect(paths, ["phi"],\
                                collectGhost = True,\
                                yInd = [self._yInd, self._yInd],\
                                )["phi"]

        # The ExB velocity for a Clebsch system can be found in section B.5.1
        # in the BOUT++ coordinates manual. However, the cylindrical
        # coordinate system is not a Clebsch system, but the metric overlaps.
        # In order to get the cylindrical coordinates ExB, we must multiply
        # the ExB velocity in BOUT++ with B (i.e. divide by rho). Thus, the
        # radial ExB velocity is the cylindrical theta derivative of phi
        self._radialExB = DDZ(phi, self._J)
    #}}}

    #{{{getRadialProbeIndices
    def getRadialProbeIndices(self, indexIn, nProbes = 5):
        """
        Get rho indices for nProbes located in an symmetric, equidistant way
        around the input indexIn

        Parameters
        ----------
        indexIn : int
            The index to put the probes around
        nProbes : int
            Number of probes (including indexIn). Note that even
            numbers here will be converted to nearest odd number below.

        Returns
        -------
        indices : list
            A list of the rho index to put the probes on (including indexIn)
        """

        # Find out if we are above half
        innerLen = self._var.shape[1] - 2*self._MXG
        pctOfInd = indexIn/innerLen

        # We here find the span of available indices to put the probes at
        if pctOfInd < 0.5:
            indexSpan = indexIn
        else:
            indexSpan = innerLen - indexIn

        # Floor in order not to get indices out of bounds
        halfNProbes = int(np.floor((nProbes - 1)/2))

        # Pad with one probe which we do not use
        halfNProbes += 1

        # Index sepration from indexIn
        # Floor in order not to get indices out of bounds
        indexSep = int(np.floor(indexSpan/halfNProbes))

        indices = [indexIn]

        for i in range(1, halfNProbes):
            # Insert before
            indices.insert(0, indexIn - i*indexSep)
            # Insert after
            indices.append(indexIn + i*indexSep)

        return indices
    #}}}
#}}}

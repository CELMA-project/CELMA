#!/usr/bin/env python

"""
Contains classes which probes the data
"""

from .polAvg import polAvg
from ..plotHelpers import (PlotHelper,\
                           collectiveCollect,\
                           DDZ,\
                           findLargestRadialGrad)
import numpy as np
from scipy.stats import kurtosis, skew
from scipy.signal import periodogram
from boutdata import collect

#{{{class Probes
class Probes(object):
    #{{{docstring
    """
    Parent class which stores the Jacobian J, the variable and the
    time traces.

    Contains analysis functions.
    """
    #}}}

    #{{{Constructor
    def __init__(self                       ,\
                 var                        ,\
                 varName                    ,\
                 time                       ,\
                 tIndSaturatedTurb   = None ,\
                 steadyStatePath     = None ,\
                 radialProbesIndices = None ,\
                 collectPath         = None ,\
                 convertToPhysical   = False,\
                 ):
        #{{{docstring
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
        tIndSaturatedTurb : [int|None]
            Index at where the turbulence saturates. This should be set
            to be after the overshoot in the energy.
        steadyStatePath : str
            What path to use when collecting J. If radialProbesIndices is
            None, this will also be the path for finding the largest
            gradient.
            Default is None.
        radialProbesIndices : [None|array]
            What radial indices to probe. If set to None, it will be
            selected from the largest gradient in the steady state
            variable.
        collectPath : str
            Path to collect J, coordinates and the normalization constants
            from. Not effective if steadyStatePath is set.
        convertToPhysical : str
            Whether normalized or physical units should be used.
        """
        #}}}

        # Guard
        if steadyStatePath is None and radialProbesIndices is None:
            message = ("steadyStatePath and "
                       "radialProbesIndices cannot be None simultaneously")
            raise ValueError(message)
        if steadyStatePath is None and collectPath is None:
            message = ("steadyStatePath and "
                       "collectPath cannot be None simultaneously")
            raise ValueError(message)

        # Set the collect path
        if collectPath is None:
            collectPath = steadyStatePath

        # Make the PlotHelper object
        # Public as used in the driver
        self.helper = PlotHelper(collectPath                           ,\
                                  t                 = time             ,\
                                  xguards           = False            ,\
                                  yguards           = False            ,\
                                  convertToPhysical = convertToPhysical,\
                                 )

        # Get the units (eventually convert to physical units)
        self._var, self.varNormalization, self.varUnits =\
            self.helper.physicalUnitsConverter(var, varName)

        self.time      = time
        self.fluctTime = time[tIndSaturatedTurb:]

        # Find the fluctuations in var
        self._varAvg   = polAvg(var[tIndSaturatedTurb:, :, :, :])
        self._varFluct = var[tIndSaturatedTurb:, :, :, :] - self._varAvg

        if varName == "n":
            collectVarName = "lnN"

        # Set the Jacobian
        # Contains the ghost points as we are using this in DDZ
        self._J = collect("J", path=collectPath,\
                          xguards=True, yguards=True, info=False)

        if radialProbesIndices == None:
            # Note that the ghost cells are collected, as we are taking
            # derivatives of the field
            self._varSteadyState =\
                collect(collectVarName, path=collectPath,\
                        xguards=True, yguards=True, info=False)

            if varName == "n":
                self._varSteadyState = np.exp(self._varSteadyState)

        self.varName           = varName
        self.tIndSaturatedTurb = tIndSaturatedTurb

        # Set uninitialized variables to None
        self.results             = None
        self.timeTraceOfVar      = None
        self.timeTraceOfVarAvg   = None
        self.timeTraceOfVarFluct = None
        self._xInds              = None
        self.yInd                = None
        self._yInds              = None
        self._zInds              = None
    #}}}

    #{{{initializeInputOutput
    def initializeInputOutput(self, xInds, yInds, zInds):
        #{{{docstring
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
        #}}}

        self.results             = {}
        self.timeTraceOfVar      = {}
        self.timeTraceOfVarAvg   = {}
        self.timeTraceOfVarFluct = {}

        self._xInds = xInds
        if self.yInd is not None:
            # If the perpPlaneClass has selected the variable
            actualYInds = [self.yInd]
            yInds       = [0]
        else:
            actualYInds = [yInds]

        self._yInds       = yInds
        self._actualYInds = actualYInds
        self._zInds       = zInds

        # Copy the coordinates so that they will have the same
        # dimensions as timetraces
        rho   = self.helper.rho  .copy()
        theta = self.helper.theta.copy()
        z     = self.helper.z    .copy()

        self.rho   = {}
        self.theta = {}
        self.z     = {}

        for xInd in self._xInds:
            for yInd, actualYInd in zip(self._yInds, self._actualYInds):
                for zInd in self._zInds:
                    # Initialize the time traces
                    self.timeTraceOfVar\
                            ["{},{},{}".format(xInd, actualYInd, zInd)] =\
                                self._var[:, xInd, yInd, zInd]

                    self.timeTraceOfVarAvg\
                            ["{},{},{}".format(xInd, actualYInd, zInd)] =\
                                self._varAvg[:, xInd, yInd, zInd]

                    self.timeTraceOfVarFluct\
                            ["{},{},{}".format(xInd, actualYInd, zInd)] =\
                                self._varFluct[:, xInd, yInd, zInd]

                    # Initialize the result as a dictionary
                    self.results["{},{},{}".format(xInd, actualYInd, zInd)] =\
                        {}

                    # Set the coordinates
                    self.rho["{},{},{}".format(xInd, actualYInd, zInd)] =\
                        rho[xInd]
                    self.theta["{},{},{}".format(xInd, actualYInd, zInd)] =\
                        theta[zInd]
                    self.z["{},{},{}".format(xInd, actualYInd, zInd)] =\
                        z[actualYInd]

        # Organize the keys
        probesKeys = list(self.results.keys())

        # Sort the probesKeys list (this cannot be done alphabethically
        # Instead we cast the key into a number, and sort it from that
        probesKeysForSorting = [int(el.replace(",","")) for el in probesKeys]
        self.probesKeys =\
            [list(el) for el in zip(*sorted(\
                zip(probesKeysForSorting, probesKeys), key=lambda pair:
                pair[0]))][1]
    #}}}

    #{{{calcStatsMoments
    def calcStatsMoments(self):
        #{{{docstring
        """
        Calculates the first statistical moments.

        Output
        ------
        The output will be stored in the keys (specified below) under
        self._result[indexString], where indexString is the string of
        the index under consideration.

        fluctMean : array
            The mean of the fluctuations for each time (should be 0).
        fluctVar : array
            The variance of the fluctuations for each time.
        fluctKurt : array
            The kurtosis of the fluctuations for each time (is 3.0 for a
            normal distribution).
        fluctSkew : array
            The skewness of the the fluctuations for each time (is 0 for
            a normal distribution)
        """
        #}}}

        for xInd in self._xInds:
            for yInd, actualYInd in zip(self._yInds, self._actualYInds):
                for zInd in self._zInds:
                    key = "{},{},{}".format(xInd, actualYInd, zInd)
                    self.results[key]["fluctMean"] =\
                        self.timeTraceOfVarFluct[key].mean()
                    self.results[key]["fluctVar"] =\
                        self.timeTraceOfVarFluct[key].var()
                    self.results[key]["fluctKurt"] =\
                        kurtosis(self.timeTraceOfVarFluct[key])
                    self.results[key]["fluctSkew"] =\
                        skew(self.timeTraceOfVarFluct[key])
    #}}}

    #{{{calcPDFs
    def calcPDFs(self):
        #{{{docstring
        """
        Calculates the probability distribution function (PDF) of the
        fluctuations.

        Probability distribution function
        ---------------------------------
        Probability that the measurement falls within an infinite small
        interval.

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
        #}}}

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
            for yInd, actualYInd in zip(self._yInds, self._actualYInds):
                for zInd in self._zInds:
                    key = "{},{},{}".format(xInd, actualYInd, zInd)

                    try:
                        self.results[key]["pdfY"], bins =\
                                  np.histogram(self.timeTraceOfVarFluct[key],\
                                               bins="auto",\
                                               density=True)
                        # Initialize y
                        self.results[key]["pdfX"] =\
                                np.zeros(self.results[key]["pdfY"].size)

                        for k in range(self.results[key]["pdfY"].size):
                            # Only the bin edges are saved. Interpolate to the
                            # center of the bin
                            self.results[key]["pdfX"][k] = 0.5*(bins[k]+bins[k+1])
                    except MemoryError:
                        message = ("{0}{1}WARNING: MemoryError in histogram. "
                                   "Setting manually to 1{1}{0}".format("\n"*2, "!"*3))
                        self.results[key]["pdfX"] = self.results[key]["pdfY"] = [0.99,1]
                        print(message)
    #}}}

    #{{{calcPSDs
    def calcPSDs(self):
        #{{{docstring
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
        #}}}

        # Guard
        if self.results is None:
            message = ("Results not initiazled. Please do so by calling "
                       "initializeInputOutput")
            raise RuntimeError(message)

        try:
            fs = self.fluctTime[1] - self.fluctTime[0]
        except IndexError as ie:
            if "out of bounds" in ie.args[0]:
                message = ("{0}{1}WARNING Specified tIndSaturatedTurb was out of "
                           "range when calculating PSD.{1}{0}")
                print(message.format("\n", "!"*3))
            setToNone = True

        if setToNone:
            for xInd in self._xInds:
                for yInd, actualYInd in zip(self._yInds, self._actualYInds):
                    for zInd in self._zInds:
                        key = "{},{},{}".format(xInd, actualYInd, zInd)

                        # window = None => window = "boxcar"
                        # scaling = density gives the correct units
                        self.results[key]["psdX"] = None
        else:
            for xInd in self._xInds:
                for yInd, actualYInd in zip(self._yInds, self._actualYInds):
                    for zInd in self._zInds:
                        key = "{},{},{}".format(xInd, actualYInd, zInd)

                        # window = None => window = "boxcar"
                        # scaling = density gives the correct units
                        self.results[key]["psdX"], self.results[key]["psdY"] =\
                            periodogram(self.timeTraceOfVarFluct[key],\
                                        fs=fs, window=None, scaling="density")
    #}}}

    #{{{calcAvgFluxThroughVolumeElement
    def calcAvgFluxThroughVolumeElement(self, u, uName):
        #{{{docstring
        """
        Gives the average flux through a volume element.

        Note that this is different from the total flux through a surface,
        which is the surface integral of the same quantity.
        However, this is usually the flux through a voluem element which are
        measured with probes in experiments.

        A note on poloidal averaging:
        -----------------------------
        <ab> = < (<a> + a_fluct)(<b> + b_fluct) >
             = < <a> <b>+ a_fluct<b> + <a> b_fluct + a_fluct b_fluct>
             = < <a> <b> > + < a_fluct<b> > + < <a> b_fluct> +<a_fluct b_fluct>
             = < <a> <b> > + <b>< a_fluct > + <a>< b_fluct> + <a_fluct b_fluct>
             = < <a> <b> > + <b>0 + <a>0 + <a_fluct b_fluct>
             = < <a> <b> > + <a_fluct b_fluct>

        Parameters
        ----------
        u : array
            The velocity of the flux in the direction of the flux
            (xInd and yInd must be specified)
        uName : str
            Name of the velocity

        Output
        ------
        The output will be stored in the keys (specified below) under
        self._result[indexString], where indexString is the string of
        the index under consideration.
        Notice that the result will be the same for every z-index for a
        fixed x- and y-index.

        varFlux* : array
            The average flux at the (xInd, yInd, zInd) position for each
            time.
        varFluxAvg* : array
            The average flux arising from the averaged fields at the (xInd,
            yInd, zInd) position for each time.
        varFluxFluct* : array
            The average flux arising from the fluctuations at the (xInd,
            yInd, zInd) position for each time.

        * The keys below will be appended with uName.
        """
        #}}}

        # Guard
        if self.results is None:
            message = ("Results not initiazled. Please do so by calling "
                       "initializeInputOutput")
            raise RuntimeError(message)

        # Convert to physical, and get units
        u, self._uNormalization, self._uUnits =\
            self.helper.physicalUnitsConverter(u, "u")

        # Find the fluctuating velocity
        uAvg         = polAvg(u)
        uFluct       = u - uAvg

        for xInd in self._xInds:
            for yInd, actualYInd in zip(self._yInds, self._actualYInds):
                avgFlux      = polAvg(\
                                self._var     [:,xInd:xInd+1,yInd:yInd+1,:]*\
                                u             [:,xInd:xInd+1,yInd:yInd+1,:])
                avgFluxAvg   = polAvg(\
                                self._varAvg  [:,xInd:xInd+1,yInd:yInd+1,:]*\
                    uAvg[self.tIndSaturatedTurb:,xInd:xInd+1,yInd:yInd+1,:])
                avgFluxFluct = polAvg(\
                                self._varFluct[:,xInd:xInd+1,yInd:yInd+1,:]*\
                    uFluct[self.tIndSaturatedTurb:,xInd:xInd+1,yInd:yInd+1,:])
                for zInd in self._zInds:

                    key = "{},{},{}".format(xInd, actualYInd, zInd)

                    self.results[key]["varAvgFlux" + uName.capitalize()] =\
                        avgFlux[:, 0, 0, zInd]

                    self.results[key]["varAvgFluxAvg" + uName.capitalize()] =\
                        avgFluxAvg[:, 0, 0, zInd]

                    self.results[key]["varAvgFluxFluct"+ uName.capitalize()] =\
                        avgFluxFluct[:, 0, 0, zInd]
    #}}}

    #{{{calcFFTs
    def calcFFTs(self):
        #{{{docstring
        """
        Calculates the FFT of the poloidal profile belonging to the
        point under consideration.

        Output
        ------
        The output will be stored in the keys (specified below) under
        self._result[indexString], where indexString is the string of
        the index under consideration.

        zFFT : array
            The fourier transformed of the z-direction for each time. Notice
            that the result will be the same for every z-index for a fixed x-
            and y-index.
        zFFTNonSaturatedIndex : int
            Index for the end of the non-saturated phase. The end is here
            defined as the first point where the max of all the modes is
            15% of the max value.
        zFFTLinearIndex : int
            Index for the end of the linear phase. The end is here
            defined as the index where the first mode which reaches
            1e-7 (excluding the clip).
        """
        #}}}

        # Guard
        if self.results is None:
            message = ("Results not initiazled. Please do so by calling "
                       "initializeInputOutput")
            raise RuntimeError(message)

        for xInd in self._xInds:
            for yInd, actualYInd in zip(self._yInds, self._actualYInds):
                varFFT =\
                    np.fft.fft(\
                        self._var[:, xInd:xInd+1, yInd:yInd+1,:],\
                        axis=-1)
                for zInd in self._zInds:
                    key = "{},{},{}".format(xInd, actualYInd, zInd)
                    # Save the results and reshape the data
                    self.results[key]["zFFT"] = varFFT[:,0,0,:]

        # Find the non saturated phase end and linear phase end for each key
        fracOfMax           = 0.15
        firstIndexEndLinear = self.results[key]["zFFT"].shape[0]
        # Do not take into account the three first time steps (which is
        # after the initial perturbation)
        clip = 3
        for key in self.results.keys():
            curMax = 0
            # Skip the offset mode in range
            for mode in range(1, self.results[key]["zFFT"].shape[-1]):
                #{{{ NOTE: We are dealing with a real signal:
                #          As the fourier transform breaks the signal up
                #          in cisoids there will be one part of the
                #          signal in the positive rotating ciscoid and
                #          one in the negative (negative frequencies)
                #          for a given mode number. We need to take into
                #          account both in order to calculate the
                #          amplitude. As the signal is real only one of
                #          the phase sifts are needed. Notice that for a
                #          real signal the imaginary part occurs as a
                #          complex conjugate pair
                # http://dsp.stackexchange.com/questions/431/what-is-the-physical-significance-of-negative-frequencies?noredirect=1&lq=1
                # http://dsp.stackexchange.com/questions/4825/why-is-the-fft-mirrored
                #}}}
                # Magnitude of the signal
                # https://en.wikipedia.org/wiki/Discrete_Fourier_transform#Definition
                magnitude = (np.abs(np.abs(self.results[key]["zFFT"][clip:, mode]))
                            +np.abs(np.abs(self.results[key]["zFFT"][clip:,-mode]))
                            )/self.results[key]["zFFT"].shape[-1]
                # Non saturated phase
                maxOfThisMode = np.max(magnitude)
                if maxOfThisMode > curMax:
                    curMax      = maxOfThisMode
                    modeWithMax = mode
                # Linear phase
                curIndicesEndLinear = np.where(magnitude >= 1e-7)
                if len(curIndicesEndLinear) > 0:
                    if curIndicesEndLinear < firstIndexEndLinear:
                        firstIndexEndLinear = curIndicesEndLinear 

            self.results[key]["zFFTLinearIndex"] = firstIndexEndLinear

            try:
                # Find the first occurence where the mode is above or
                # equal to 15%
                self.results[key]["zFFTNonSaturatedIndex"] =\
                    int(np.where(np.abs(\
                        self.results[key]["zFFT"][clip:,modeWithMax]) >\
                        curMax*fracOfMax\
                    )[0])
            except TypeError as er:
                if "only length-1 arrays" in er.args[0]:
                    # Need to subscript once more 
                    self.results[key]["zFFTNonSaturatedIndex"] =\
                        int(np.where(np.abs(\
                            self.results[key]["zFFT"][clip:,modeWithMax]) >\
                            curMax*fracOfMax\
                        )[0][0])
                else:
                    raise er
    #}}}
#}}}

#{{{class PerpPlaneProbes
class PerpPlaneProbes(Probes):
    #{{{docstring
    """
    Child class of probes. Calls probes constructor, finds
    radialProbeIndice based on highest gradient (if not set).
    Also calculates the ExB velocity.

    This class makes it possible to collect only a z-plane instead of
    the whole 3D variable.
    """
    #}}}

    #{{{Constructor
    def __init__(self                          ,\
                 varName                       ,\
                 paths                         ,\
                 yInd                          ,\
                 nProbes                = 5    ,\
                 radialProbesIndices    = None ,\
                 **kwargs):
        #{{{docstring
        """
        Constructor for the PerpPlaneProbes class

        Note that ghost points are collected, as we are finding
        derivatives of the fields.

        Parameters
        ----------
        varName : string
            Name of the variable.
        paths : iterable of strings
            What path to use when collecting the variable. Must be in
            ascending temporal order as the variable will be
            concatenated.
        yInd : int
            yInd to collect from
        nProbes : int
            Number of probes. Default is 5. Will be overridden by
            radialProbesIndices if set.
        radialProbesIndices : [None|array]
            What radial indices to probe. If set to None, it will be
            selected from the largest gradient in the steady state
            variable.
        **kwargs : keyword arguments
            Keyword arguments used in the parent constructor. See the
            Probes constructor for details.
        """
        #}}}

        if varName == "n":
            collectVarName = "lnN"

        # Collect the variable
        varTimeDict = collectiveCollect(paths, [collectVarName, "t_array"],\
                                        collectGhost = False,\
                                        yInd = [yInd, yInd],\
                                        )

        var  = varTimeDict[collectVarName]
        time = varTimeDict["t_array"]

        if varName == "n":
            var = np.exp(var)

        # Call the parent class
        super().__init__(var                                      ,\
                         varName                                  ,\
                         time                                     ,\
                         radialProbesIndices = radialProbesIndices,\
                         **kwargs)

        self.yInd = yInd

        if radialProbesIndices is None:
            # Collect dx and MXG
            # xguards will be collected as derivatives will be taken
            dx        = collect("dx",  path=paths[0], xguards=True, info=False)
            self._MXG = collect("MXG", path=paths[0], info=False)
            # Find the max gradient of the variable (subtracts the guard cells)
            _, maxGradInd =\
                findLargestRadialGrad(\
                  self._varSteadyState[0:1, :, self.yInd:self.yInd+1, 0:1],\
                  dx,\
                  self._MXG)
            self.radialProbesIndices =\
                self.getRadialProbeIndices(maxGradInd, nProbes)
        else:
            self.radialProbesIndices = radialProbesIndices

        # Get the radial ExB drift for this plane
        # NOTE: We do not convert this to physical units but the ExB
        #       drift
        phi = collectiveCollect(paths, ["phi"],\
                                collectGhost = True,\
                                yInd = [self.yInd, self.yInd],\
                                )["phi"]

        # The ExB velocity for a Clebsch system can be found in section B.5.1
        # in the BOUT++ coordinates manual. However, the cylindrical
        # coordinate system is not a Clebsch system, but the metric overlaps.
        # In order to get the cylindrical coordinates ExB, we must multiply
        # the ExB velocity in BOUT++ with B (i.e. divide by rho). Thus, the
        # radial ExB velocity is the cylindrical theta derivative of phi
        self.radialExB = DDZ(phi, self._J)
    #}}}

    #{{{getRadialProbeIndices
    def getRadialProbeIndices(self, indexIn, nProbes = 5):
        #{{{docstring
        """
        Get rho indices for nProbes located in an symmetric, equidistant way
        around the input indexIn.

        NOTE: Does not work if indexIn = 0.

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
        #}}}

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

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
    #{{{docstring
    """
    Parent class which stores the Jacobian J, the variable and the
    time traces.

    Contains analysis functions.
    """
    #}}}

    #{{{Constructor
    def __init__(self, var, varName, time, tIndSaturatedTurb=None,\
                 steadyStatePath=None, radialProbeIndices=None,\
                 collectPath=None):
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
            Index at where the turbulence saturates. This can be found
            from plotting self.results[index]['zFFT'] after calculation
        steadyStatePath : str
            What path to use when collecting J. If radialProbeIndices is
            None, this will also be the path for finding the largest
            gradient. Default is None.
        radialProbeIndices : [None|array]
            What radial indices to probe. If set to None, it will be
            selected from the largest gradient in the steady state
            variable.
        collectPath : str
            Path to collect J and the coordinates from. Not effective
            if steadyStatePath is set.
        """
        #}}}

        # Guard
        if steadyStatePath is None and radialProbeIndices is None:
            message = ("steadyStatePath and radialProbeIndices cannot "
                       "be None simultaneously")
            raise ValueError(message)
        if steadyStatePath is None and collectPath is None:
            message = ("steadyStatePath and collectPath cannot "
                       "be None simultaneously")
            raise ValueError(message)

        self._var      = var
        self.time      = time
        self.fluctTime = time[tIndSaturatedTurb:]

        # Find the fluctuations in var
        self._varAvg   = polAvg(var)
        self._varFluct = var - self._varAvg

        # Clip the fluctuation part
        self._varFluct = self._varFluct[tIndSaturatedTurb:,:,:,:]

        if varName == "n":
            collectVarName = "lnN"

        if collectPath is None:
            collectPath = steadyStatePath

        # Set the Jacobian
        # Contains the ghost points as we are using this in DDZ
        self._J = collect("J", path=collectPath,\
                          xguards=True, yguards=True, info=False)

        # Sets the normalized coordinates
        # Get the coordinates
        #{{{rho
        self._dx = collect('dx', path = collectPath,\
                           xguards = True, yguards = True, info = False)
        self._MXG = collect('MXG', path = collectPath,\
                            xguards = True, yguards = True, info = False)

        nPoints  = self._dx.shape[0]
        dx       = self._dx[0,0]

        innerPoints = nPoints - 2*self._MXG
        self.rho    = dx * np.array(np.arange(0.5, innerPoints))

        # Insert the first and last grid point due to the guards
        self.rho = np.insert(self.rho, 0, - 0.5*dx)
        self.rho = np.append(self.rho, self.rho[-1] + dx)
        #}}}

        #{{{z
        dy = collect('dy', path = collectPath,\
                     xguards = True, yguards = True, info = False)
        MYG = collect('MYG', path = collectPath,\
                      xguards = True, yguards = True, info = False)

        nPoints  = dy.shape[1]
        self._dy = dy[0,0]

        innerPoints = nPoints - 2*MYG
        self.z = self._dy * np.array(np.arange(0.5, innerPoints))

        # Insert the first and last grid point
        self.z = np.insert(self.z, 0, - 0.5*self._dy)
        self.z = np.append(self.z, self.z[-1] + self._dy)
        #}}}

        #{{{theta
        dz = collect('dz', path = collectPath,\
                     xguards = True, yguards = True, info = False)
        MZ = collect('MZ', path = collectPath,\
                     xguards = True, yguards = True, info = False)

        # Subtract the unused plane
        innerPoints = MZ - 1

        self.theta = dz * np.array(np.arange(0.0, innerPoints))
        #}}}

        if radialProbeIndices == None:
            # Note that the ghost cells are collected, as we are taking
            # derivatives of the field
            self._varSteadyState =\
                collect(collectVarName, path=collectPath,\
                        xguards=True, yguards=True, info=False)

            if varName == "n":
                self._varSteadyState = np.exp(self._varSteadyState)

        self.varName            = varName
        self._tIndSaturatedTurb = tIndSaturatedTurb

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
        rho   = self.rho  .copy()
        theta = self.theta.copy()
        z     = self.z    .copy()

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
                    self.results["{},{},{}".format(xInd, actualYInd, zInd)] = {}

                    # Set the coordinates
                    self.rho["{},{},{}".format(xInd, actualYInd, zInd)] =\
                        rho[xInd]
                    self.theta["{},{},{}".format(xInd, actualYInd, zInd)] =\
                        theta[zInd]
                    self.z["{},{},{}".format(xInd, actualYInd, zInd)] =\
                        z[actualYInd]
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
        #}}}

        for xInd in self._xInds:
            for yInd, actualYInd in zip(self._yInds, self._actualYInds):
                for zInd in self._zInds:
                    key = "{},{},{}".format(xInd, actualYInd, zInd)
                    self.results[key]["mean"] =\
                        self.timeTraceOfVarFluct[key].mean()
                    self.results[key]["var"] =\
                        self.timeTraceOfVarFluct[key].var()
                    self.results[key]["kurtosis"] =\
                        kurtosis(self.timeTraceOfVarFluct[key])
                    self.results[key]["skew"] =\
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

                    self.results[key]["pdfY"], bins =\
                              np.histogram(self.timeTraceOfVarFluct[key],\
                                           bins="auto",\
                                           density=True)
                    # Initialize y
                    self.results[key]["pdfX"] =\
                            np.zeros(self.results[key]["pdfY"].size)

                    for k in range(self.results[key]["pdfY"].size):
                        # Only the bin edges are saved. Interpolate to the center of the bin
                        self.results[key]["pdfX"][k] = 0.5*(bins[k]+bins[k+1])
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

        fs = self.time[1] - self.time[0]
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

# FIXME: Bug here
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
             = < <a> <b> > + < a_fluct<b> > + < <a> b_fluct> + <a_fluct b_fluct>
             = < <a> <b> > + <b>< a_fluct > + <a>< b_fluct> + <a_fluct b_fluct>
             = < <a> <b> > + <b>0 + <a>0 + <a_fluct b_fluct>
             = < <a> <b> > + <a_fluct b_fluct>

        Parameters
        ----------
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
                                uAvg          [:,xInd:xInd+1,yInd:yInd+1,:])
                avgFluxFluct = polAvg(\
                                self._varFluct[:,xInd:xInd+1,yInd:yInd+1,:]*\
                                uFluct        [:,xInd:xInd+1,yInd:yInd+1,:])
                for zInd in self._zInds:

                    key = "{},{},{}".format(xInd, actualYInd, zInd)

                    self.results[key]["varAvgFlux" + uName.capitalize()] =\
                        avgFlux[:, 0, 0, zInd]

                    self.results[key]["varAvgFluxAvg" + uName.capitalize()] =\
                        avgFluxAvg[:, 0, 0, zInd]

                    self.results[key]["varAvgFluxFluct" + uName.capitalize()] =\
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
            The fourier transformed of the z-direction for each time. Notice that the
            result will be the same for every z-index for a fixed x- and
            y-index.
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
                    np.abs(np.fft.fft(\
                        self._var[:, xInd:xInd+1, yInd:yInd+1,:],\
                        axis=-1))
                for zInd in self._zInds:
                    key = "{},{},{}".format(xInd, actualYInd, zInd)
                    self.results[key]["zFFT"] = varFFT
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
    def __init__(self, varName, paths, yInd,\
                 nProbes=5, physicalUnits=False, tIndSaturatedTurb=None,\
                 steadyStatePath=None, radialProbeIndices=None):
        #{{{docstring
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
        tIndSaturatedTurb : [int|None]
            Index at where the turbulence saturates. This can be found
            from plotting self.results[index]['zFFT'] after calculation
        steadyStatePath : string
            What path to use when collecting J. If radialProbeIndices is
            None, this will also be the path for finding the largest
            gradient. Default is None.
        radialProbeIndices : [None|array]
            What radial indices to probe. If set to None, it will be
            selected from the largest gradient in the steady state
            variable.
        """
        #}}}

        # Guard
        if steadyStatePath is None and radialProbeIndices is None:
            message = ("steadyStatePath and radialProbeIndices cannot "
                       "be None simultaneously")
            raise ValueError(message)

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

        # Collect and recalculate the time
        if physicalUnits:
            try:
                self.n0   = collect("n0" ,  path = paths[-1])
                self.Te0  = collect("Te0",  path = paths[-1])
                self.Ti0  = collect("Ti0",  path = paths[-1])
                self.B0   = collect("B0" ,  path = paths[-1])
                self.omCI = collect("omCI", path = paths[-1])
                self.rhoS = collect("rhoS", path = paths[-1])

                time /= self.omCI
            except ValueError:
                # An OSError is thrown if the file is not found
                message = ("{0}{1}WARNING: Normalized quantities not found. "
                           "The time remains normalized".format("\n"*3,"!"*3))
                print(message)

        # Call the parent class
        super().__init__(var,\
                         varName,\
                         time,\
                         tIndSaturatedTurb,\
                         steadyStatePath,\
                         radialProbeIndices)

        self.yInd = yInd

        if radialProbeIndices is None:
            # Find the max gradient of the variable
            _, maxGradInd =\
                findLargestRadialGrad(\
                  self._varSteadyState[0:1, :, self.yInd:self.yInd+1, 0:1],\
                  self._dx,\
                  self._MXG)
            self.radialProbesIndices =\
                self.getRadialProbeIndices(maxGradInd, nProbes)
        else:
            self.radialProbesIndices = radialProbeIndices

        # Get the radial ExB drift for this plane
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

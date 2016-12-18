#!/usr/bin/env python

""" Contains the UnitsConverter class """

from ..collectAndCalcHelpers import safeCollect
import scipy.constants as cst

#{{{UnitsConverter
class UnitsConverter(object):
    """
    Class which deals with units and conversion
    """

    #{{{__init__
    def  __init__(self                    ,\
                  path                    ,\
                  convertToPhysical = True):
        #{{{docstring
        """
        The constructor for UnitsConverter, which:

        * Sets the member data
        * Collects the coordinates (excluding t)
        * Makes the conversion dictionary

        Parameters
        ----------
        path : str
            The path to collect from.
        t : array
            The time array.
        useSpatial : bool
            Whether or not to collect the spatial domain.
        xguards : bool
            If xguards should be included when collecting.
        yguards : bool
            If yguards should be included when collecting.
        convertToPhysical : bool
            Whether or not to convert to physical units.
        """
        #}}}

        # Set the member data
        self._path             = path
        self.convertToPhysical = convertToPhysical

        # Get the normalizer dict
        self._normDict = self._getNormalizerDict()

        # Makes the conversion dict
        self._makeConversionDict()
    #}}}

    #{{{_getNormalizerDict
    def _getNormalizerDict(self):
        #{{{docstring
        """
        Get the conversion dict.

        Sets convertToPhysical to False if the normalization parameters
        are not found.

        Returns
        -------
        normalizerDict : dictionary
        """
        #}}}

        normalizerDict = {}
        normalizers = ("omCI", "rhoS", "n0", "Te0")

        # Caclulate mi from mu
        normalizerDict["mi"] =\
                safeCollect("mu", path=self._path, info=False)*cst.m_e

        if self.convertToPhysical:
            # If no errors are encountered, the function returns normalizerDict
            try:
                for normalizer in normalizers:
                    normalizerDict[normalizer] =\
                            safeCollect(normalizer, path=self._path, info=False)

                # The collected Te0 is given in eV, we convert this to J
                normalizerDict["Te0"].setflags(write=True)
                normalizerDict["Te0"] *= cst.e
                normalizerDict["Te0"].setflags(write=False)

                return normalizerDict

            except ValueError as ve:
                if "not found" in ve.args[0]:
                    # An OSError is thrown if the file is not found
                    message = ("{0}{1}WARNING: {2} not found. "\
                               "Variables remains normalized"\
                               ).format("\n"*3, "!"*3, normalizer)
                    print(message)

                    # Reset convertToPhysical
                    self.convertToPhysical = False
                else:
                    raise ve

        # If convertToPhysical and no errors occured, the function have
        # already returned
        for normalizer in normalizers:
            normalizerDict[normalizer] = 1

        return normalizerDict
    #}}}

    #{{{_makeConversionDict
    def _makeConversionDict(self):
        #{{{docstring
        """
        Makes the physical converter dict
        The dictionary contains of the keys (and does not contain the
        '$' for LaTeX)
            * units         - String with the units of the physical quantity
            * normalization - String with the normalization
            * factor        - Multiplying factor to convert from
                              normalized quantity to physical units
                              quantity
            * normFactor    - Multiplying factor to convert a quantity
                              which is saved as a physical unit in the
                              simulation to a normalized quantity.
        """
        #}}}

        self.conversionDict = {\
        # NOTE: logarithm is without dimension
        "lnN"       :{"units"        :r"{}",\
                      "normalization":r"{}",\
                      "factor"       :1,\
                     },\
        # NOTE: n0 is input parameter, but n is from an evolving field
        "n"         :{"units"        :r"\mathrm{m}^{-3}",\
                      "normalization":r"/n_0",\
                      "factor"       :self._normDict["n0"],\
                     },\
        # NOTE: nn is an input parameter, and is thus already given in physcial
        #       units
        "nn"        :{"units"        :r"\mathrm{m}^{-3}",\
                      "normalization":r"/n_0",\
                      "factor"       :1,\
                      "normFactor"   :1/self._normDict["n0"],\
                     },\
        "vort"      :{"units"        :r"\mathrm{s}^{-1}",\
                      "normalization":r"/\omega_{{ci}}",\
                      "factor"       :self._normDict["omCI"],\
                     },\
        "vortD"     :{"units"        :r"\mathrm{m}^{-3}\mathrm{s}^{-1}",\
                      "normalization":r"/\omega_{{ci}}n_0",\
                      "factor"       :self._normDict["omCI"]*\
                                      self._normDict["n0"],\
                     },\
        # NOTE: B0 is an input parameter, and is thus already given in physcial
        #       units
        "B0"        :{"units"        :r"\mathrm{T}",\
                      "normalization":"",\
                      "factor"       :1,\
                     },\
        "phi"       :{"units"        :r"\mathrm{J}\mathrm{C}^{-1}",\
                      "normalization":r" q/T_{{e,0}}",\
                      "factor"       :self._normDict["Te0"]/cst.e,\
                     },\
        "jPar"      :{"units"        :r"\mathrm{C}\mathrm{s}^{-1}",\
                      "normalization":r"/n_0c_sq",\
                      "factor"       :cst.e*\
                                      self._normDict["rhoS"]*\
                                      self._normDict["omCI"]*\
                                      self._normDict["n0"],\
                     },\
        #NOTE: momDensPar is divided by mi, so we need to multiply by mi again
        #      here
        "momDensPar":{"units"        :r"\mathrm{kg\;m}^{-2}\mathrm{\;s}^{-1}",\
                      "normalization":r"/m_in_0c_s",\
                      "factor"       :self._normDict["mi"]*\
                                      self._normDict["rhoS"]*\
                                      self._normDict["omCI"]*\
                                      self._normDict["n0"],\
                     },\
        "uIPar"     :{"units"        :r"\mathrm{ms}^{-1}",\
                      "normalization":r"/c_s",\
                      "factor"       :self._normDict["rhoS"]*\
                                      self._normDict["omCI"],\
                     },\
        "uEPar"     :{"units"        :r"\mathrm{ms}^{-1}",\
                      "normalization":r"/c_s",\
                      "factor"       :self._normDict["rhoS"]*\
                                      self._normDict["omCI"],\
                     },\
        "u"         :{"units"        :r"\mathrm{ms}^{-1}",\
                      "normalization":r"/c_s",\
                      "factor"       :self._normDict["rhoS"]*\
                                      self._normDict["omCI"],\
                     },\
        "S"         :{"units"        :r"\mathrm{m}^{-3}\mathrm{s}^{-1}",\
                      "normalization":r"/\omega_{{ci}}n_0",\
                      "factor"       :self._normDict["omCI"]*\
                                      self._normDict["n0"],\
                     },\
        "t"         :{"units"        :r"\mathrm{s}",\
                      "normalization":r"\omega_{{ci}}",\
                      "factor"       :1/self._normDict["omCI"],\
                     },\
        # NOTE: The growth rates are in physical units if the time is in
        #       physical units.
        #       We are here assuming that the time is in physical units
        "growthRate":{"units"        :r"1/\mathrm{s}",\
                      "normalization":r"\omega_{{ci}}",\
                      "factor"       :1,\
                     },\
        "rho"       :{"units"        :r"\mathrm{m}",\
                      "normalization":r"/\rho_s",\
                      "factor"       :self._normDict["rhoS"],\
                     },\
        "z"         :{"units"        :r"\mathrm{m}",\
                      "normalization":r"/\rho_s",\
                      "factor"       :self._normDict["rhoS"],\
                     },\
        # NOTE: len is an input parameter (Ly), and is thus already given in
        #       physcial units
        "length"    :{"units"        :r"\mathrm{m}",\
                      "normalization":r"/\rho_s",\
                      "factor"       :1,\
                      "normFactor"   :1/self._normDict["rhoS"],\
                     },\
        # NOTE: The masses are not included in the integral from the simulations
        "eEnergy"   :{"units"        :r"\mathrm{kg\; m}^2\mathrm{\; s}^{-2}",\
                      "normalization":r"/n_0T_e\rho_s^3",\
                      "factor"       :(cst.m_e/self._normDict["mi"])*\
                                      self._normDict["n0"]*\
                                      self._normDict["Te0"]*\
                                      (self._normDict["rhoS"])**3,\
                      "normFactor"   :cst.m_e/self._normDict["mi"],\
                      },\
        # NOTE: The masses are not included in the integral from the simulations
        # NOTE: mi/mi = 1
        "iEnergy"   :{"units"        :r"\mathrm{kg\; m}^2\mathrm{\; s}^{-2}",\
                      "normalization":r"/n_0T_e\rho_s^3",\
                      "factor"       :self._normDict["n0"]*\
                                      self._normDict["Te0"]*\
                                      (self._normDict["rhoS"])**3,\
                     },\
        }
    #}}}

    #{{{physicalConversion
    def physicalConversion(self, var, key):
        #{{{docstring
        """
        Convert a variable from normalized to physical units.
        Will do nothing if convertToPhysical == False

        **NOTE**: This temporarily gives write access to var

        Parameters
        ----------
        var : array
            The variable.
        key : str
            Key to use in self.conversionDict

        Returns
        -------
        var : array
            The variable after eventual processing.
        """
        #}}}
        if self.convertToPhysical:
            # Give temporarily write access
            if hasattr(var, "setflags"):
                if var.flags.owndata:
                    var.setflags(write = True)
                else:
                    # Have to make a copy
                    var = var.copy()
                    var.setflags(write = True)

            var *= self.conversionDict[key]["factor"]

            # Turn off write access
            if hasattr(var, "setflags"):
                var.setflags(write = False)

        return var
    #}}}

    #{{{normalizedConversion
    def normalizedConversion(self, var, key):
        #{{{docstring
        """
        Convert a variable from physical units to normalized.

        **NOTE**: This temporarily gives write access to var

        Parameters
        ----------
        var : array
            The variable.
        key : str
            Key to use in self.conversionDict

        Returns
        -------
        var : array
            The variable after eventual processing.
        """
        #}}}
        if self.convertToPhysical:
            # Give temporarily write access
            if hasattr(var, "setflags"):
                if var.flags.owndata:
                    var.setflags(write = True)
                else:
                    # Have to make a copy
                    var = var.copy()
                    var.setflags(write = True)

            var *= self.conversionDict[key]["normFactor"]

            # Turn off write access
            if hasattr(var, "setflags"):
                var.setflags(write = False)

        return var
    #}}}

    #{{{getUnits
    def getUnits(self, key):
        #{{{docstring
        """
        Returns the units string.

        If physicalUnitsConverter is False an empty string is returned.
        The string does not contain the '$' for LaTeX.

        Parameters
        ----------
        key : str
            Key to use in self.conversionDict

        Returns
        -------
        units : str
            The units for the input key
        """
        #}}}

        if self.convertToPhysical:
            return self.conversionDict[key]["units"]
        else:
            return ""
    #}}}

    #{{{getNormalization
    def getNormalization(self, key):
        #{{{docstring
        """
        Returns the normalization string.

        If physicalUnitsConverter is True an empty string is returned.
        The string doesnot contain the '$' for LaTeX.

        Parameters
        ----------
        key : str
            Key to use in self.conversionDict

        Returns
        -------
        normalization : str
            The normalization used for the input key
        """
        #}}}

        if self.convertToPhysical:
            return ""
        else:
            return self.conversionDict[key]["normalization"]
    #}}}
#}}}

#!/usr/bin/env python

"""
Contains a factory like function which calls the post processing
function.
"""

from .drivers1D     import Drivers1D
from .drivers2D     import Drivers2D
from .drivers1D2D   import Drivers1D2D
from .driversProbes import DriversProbes

#{{{postBoutRunner
def postBoutRunner(path, driverName = None, **kwargs):
    #{{{docstring
    """
    Function which call the driver given by driverName.

    Creates the appropriate object, and calls the appropriate driver.

    Parameters
    ----------
    path : str
        The path to collect from (required in order to be a bout_runners
        post processor).
    driverName : str
        Name of the driver to be called. All the input parameters needed
        in the function is obtained from the memberdata in the created
        object, which is created based on **kwargs.
    **kwargs : keyword arguments
        Given as input to the constructor of the object to be created.
    """
    #}}}

    #{{{Drivers1D2D
    if driverName == "plot1DAnd2DDriver"\
        or driverName == "plot1D2DAndFluctDriver":

        # Make the driver object
        driver = Drivers1D2D(path, **kwargs)

        if driverName == "plot1DAnd2DDriver":
            driver.plot1DAnd2DDriver()
        elif driverName == "plot1D2DAndFluctDriver":
            driver.plot1D2DAndFluctDriver()
    #}}}
    #{{{Drivers2D
    elif driverName == "allMainFields2DDriver"\
        or driverName == "single2DDriver":

        # Make the driver object
        driver = Drivers2D(path, **kwargs)

        # Call the driver
        if driverName == "allMainFields2DDriver":
            driver.allMainFields2DDriver()
        elif driverName == "single2DDriver":
            driver.single2DDriver()
    #}}}
    #{{{Drivers1D
    elif driverName == "parPerpDriver"\
        or driverName == "parPerpDriver"\
        or driverName == "perpDriver"\
        or driverName == "parDriver"\
        or driverName == "single1DDriver":

        # Make the driver object
        driver = Drivers1D(path, **kwargs)

        # Call the driver
        if driverName == "parPerpDriver":
            driver.parPerpDriver()
        elif driverName == "perpDriver":
            driver.perpDriver()
        elif driverName == "parDriver":
            driver.parDriver()
        elif driverName == "single1DDriver":
            driver.single1DDriver()
    #}}}
    #{{{DriversProbes
    elif driverName == "plotProbes":
        # Make the driver object
        driver = DriversProbes(path, **kwargs)
        # Call the driver
        driver.plotProbes()
    #}}}
    else:
        message = "The driverName {} is not implemented".format(driverName)
        raise NotImplementedError(message)
#}}}

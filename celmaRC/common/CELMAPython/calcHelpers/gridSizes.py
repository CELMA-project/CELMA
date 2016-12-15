"""
Contains functions dealing with sizes of the grid
"""

from boututils.datafile import DataFile
import numpy as np
import os

#{{{getSizes
def getSizes(path, coordinate, varName="lnN", includeGhost=False):
    #{{{docstring
    """
    Fastest way to obtain coordinate sizes excluding ghost points

    Parameters
    ----------
    path : str
        Path to read from
    coordinate : str
        Coordinate to return size of
    varName : str
        Field to get the size of
    includeGhost : bool
        If the ghost points should be included

    Returns
    -------
    coordinateSize : int
        Size of the desired coordinate
    """
    #}}}
    with DataFile(os.path.join(path, "BOUT.dmp.0.nc")) as f:
        if coordinate == "x":
            # nx
            coordinateSize =\
                    (f.size(varName)[1] - 2*int(f.read("MXG")))*f.read("NXPE")
            if includeGhost:
                coordinateSize += 2*int(f.read("MXG"))
        elif coordinate == "y":
            # ny
            coordinateSize =\
                    (f.size(varName)[2] - 2*int(f.read("MYG")))*f.read("NYPE")
            if includeGhost:
                coordinateSize += 2*int(f.read("MYG"))
        elif coordinate == "z":
            # nz
            coordinateSize = (f.size(varName)[3]) - 1
        else:
            raise ValueError("Unknown coordinate {}".format(coordinate))

    return coordinateSize
#}}}

#{{{getUniformSpacing
def getUniformSpacing(path, coordinate, xguards=False, yguards=False):
    #{{{docstring
    """
    Fastest way to obtain the grid spacing assuming equidistant grid

    Parameters
    ----------
    path : str
        Path to read from
    coordinate : str
        Coordinate to return size of
    xguards : bool
        If the ghost points in x should be included
    yguards : bool
        If the ghost points in y should be included

    Returns
    -------
    spacing : array-like
        The grid spacing
    """
    #}}}
    with DataFile(os.path.join(path, "BOUT.dmp.0.nc")) as f:
        if coordinate == "x" or coordinate == "y":
            if coordinate == "x":
                # dx
                spacing = f.read("dx")
            elif coordinate == "y":
                # dy
                spacing = f.read("dy")

            shape = spacing.shape
            xSize = shape[0] - 2*int(f.read("MXG"))
            ySize = shape[1] - 2*int(f.read("MYG"))
            if not(xguards) and not (yguards):
                spacingEmpty = np.empty((xSize*int(f.read("NXPE")),\
                                         ySize*int(f.read("NYPE"))))
            elif xguards and not(yguards):
                spacingEmpty = np.empty((xSize*int(f.read("NXPE"))+\
                                             2*int(f.read("MXG")),\
                                         ySize*int(f.read("NYPE"))))
            elif not(xguards) and yguards:
                spacingEmpty = np.empty((xSize*int(f.read("NXPE")),\
                                         ySize*int(f.read("NYPE"))+\
                                             2*int(f.read("MYG"))))
            elif xguards and yguards:
                spacingEmpty = np.empty((xSize*int(f.read("NXPE"))+\
                                             2*int(f.read("MXG")),\
                                         ySize*int(f.read("NYPE"))+\
                                             2*int(f.read("MYG"))))
            spacingEmpty.fill(spacing[0,0])
            spacing = spacingEmpty
        elif coordinate == "z":
            # dz
            spacing = f.read("dz")
        else:
            raise ValueError("Unknown coordinate {}".format(coordinate))

    return spacing
#}}}

#{{{getEvenlySpacedIndices
def getEvenlySpacedIndices(path, coordinate, indexIn, nPoints = 5):
    #{{{docstring
    """
    Get indices for nPoints located in an symmetric, equidistant way
    in the "coordinate" direction around the index "indexIn".

    NOTE: Does not work if indexIn = 0.

    Parameters
    ----------
    path : str
        Path to find total coordinate length from.
    coordinate : ["x"|"y"|"z"]
        The coordinate to use.
    indexIn : int
        The central index, where 0 is the first inner point.
    nPoints : int
        Number of probes (including indexIn).
        Note that even numbers here will be converted to nearest odd
        number below.

    Returns
    -------
    indices : tuple
        A tuple of the rho index to put the probes on (including indexIn)
    """
    #}}}

    # Find out if we are above half
    innerLen = getSizes(path, coordinate)
    pctOfInd = indexIn/innerLen

    # We here find the span of available indices to put the probes at
    if pctOfInd < 0.5:
        indexSpan = indexIn
    else:
        indexSpan = innerLen - indexIn

    # Floor in order not to get indices out of bounds
    halfNPoints = int(np.floor((nPoints - 1)/2))

    # Pad with one probe which we do not use
    halfNPoints += 1

    # Index sepration from indexIn
    # Floor in order not to get indices out of bounds
    indexSep = int(np.floor(indexSpan/halfNPoints))

    indices = [indexIn]

    for i in range(1, halfNPoints):
        # Insert before
        indices.insert(0, indexIn - i*indexSep)
        # Insert after
        indices.append(indexIn + i*indexSep)

    # Cast to tuple to make immutable
    indices = tuple(indices)

    return indices
#}}}

#{{{getMXG
def getMXG(path):
    #{{{docstring
    """
    Fastest way to obtain MXG

    Parameters
    ----------
    path : str
        Path to read from

    Returns
    -------
    MXG : int
        Number of ghost points in x
    """
    #}}}
    with DataFile(os.path.join(path, "BOUT.dmp.0.nc")) as f:
            return f.read("MXG")
#}}}

#{{{getMYG
def getMYG(path):
    #{{{docstring
    """
    Fastest way to obtain MYG

    Parameters
    ----------
    path : str
        Path to read from

    Returns
    -------
    MYG : int
        Number of ghost points in y
    """
    #}}}
    with DataFile(os.path.join(path, "BOUT.dmp.0.nc")) as f:
            return f.read("MYG")
#}}}

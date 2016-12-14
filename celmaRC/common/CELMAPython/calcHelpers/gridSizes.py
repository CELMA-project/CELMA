"""
Contains functions dealing with sizes of the grid
"""

#{{{getSizes
def getSizes(path, coordinate, varName="lnN"):
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

    Returns
    -------
    coordinateSize : int
        Size of the desired coordinate
    """
    #}}}
    with DataFile("BOUT.dmp.0.nc") as f:
        if coordinate == "x":
            # nx
            coordinateSize =\
                    (f.size(varName)[1] - 2*int(f.read("MXG")))*f.read("NXPE")
        elif coordinate == "y":
            # ny
            coordinateSize =\
                    (f.size(varName)[2] - 2*int(f.read("MYG")))*f.read("NYPE")
        elif coordinate == "z":
            # nz
            coordinateSize = (f.size(varName)[3]) - 1
        else:
            raise ValueError("Unknown coordinate {}".format(coordinate))

    return coordinateSize
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

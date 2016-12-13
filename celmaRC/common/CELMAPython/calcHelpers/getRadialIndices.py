"""
Contains function to obtain equidistant radial indices
"""

#{{{getRadialIndices
def getRadialIndices(nRho, MXG, indexIn, nPoints = 5):
    #{{{docstring
    """
    Get rho indices for nPoints located in an symmetric, equidistant way
    around the input indexIn.

    NOTE: Does not work if indexIn = 0.

    Parameters
    ----------
    nRho   : int
        Number of points in rho
    MXG: int
        Number of ghost points in nRho
    indexIn : int
        The index to put the probes around
    nPoints : int
        Number of probes (including indexIn). Note that even
        numbers here will be converted to nearest odd number below.

    Returns
    -------
    indices : tuple
        A tuple of the rho index to put the probes on (including indexIn)
    """
    #}}}

    # Find out if we are above half
    innerLen = nRho - 2*MXG
    pctOfInd = indexIn/innerLen

    # We here find the span of available indices to put the probes at
    if pctOfInd < 0.5:
        indexSpan = indexIn
    else:
        indexSpan = innerLen - indexIn

    # Floor in order not to get indices out of bounds
    halfNProbes = int(np.floor((nPoints - 1)/2))

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

    # Cast to tuple to make immutable
    indices = tuple(indices)

    return indices
#}}}

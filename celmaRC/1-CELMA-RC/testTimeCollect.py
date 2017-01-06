#!/usr/bin/env python

"""Checkt that the two methods of collecting the time gives the same result"""

import os, sys
# If we add to sys.path, then it must be an absolute path
commonDir = os.path.abspath("./../common")
# Sys path is a list of system paths
sys.path.append(commonDir)

from CELMAPy.collectAndCalcHelpers import collectiveCollect, collectTime
import numpy as np
import pdb;

#{{{getTimes
def getTimes(paths, tInd):
    #{{{docstring
    """
    Returns the two times

    Parameters
    ----------
    paths : tuple
        Paths to collect from
    tInd : [list|None]
        List of the time indices

    Returns
    -------
    tcc : array
        The time collected with collectiveCollect
    tct : array
        The time collected with collectTime
    """
    #}}}
    tcc =\
        collectiveCollect(paths              ,\
                          ("t_array",)       ,\
                          tInd         = tInd,\
                          )
    tcc = tcc["t_array"].flatten()
    tct = collectTime(paths, tInd = tInd)
    return tcc, tct
#}}}

#{{{verify
def verify(a, b):
    #{{{docstring
    """
    Verifies that arrays are the same

    Parameters
    ----------
    a : array
        One array
    b : array
        Second array
    """
    #}}}
    if np.allclose(a, b):
        print("Success!\n\n")
    else:
        print("Test failed, launching pdb for investigation!\n\n")
        pdb.set_trace()
        print("\n\nWill continue")
#}}}

#{{{timeCollectTest
def timeCollectTest():
    """
    Runs the test and checks if the time is elementwise the same
    """

    collectPaths =\
        ("CSDXMagFieldScanAr/nout_1000_timestep_1/geom_Lx_7.8633_geom_Ly_275.2144_input_B0_0.1_switch_forceAddNoise_False_switch_includeNoise_False_switch_saveTerms_False_switch_useHyperViscAzVortD_True_tag_CSDXMagFieldScanAr-2-linearPhase1_0",\
         "CSDXMagFieldScanAr/nout_5000_timestep_1/geom_Lx_7.8633_geom_Ly_275.2144_input_B0_0.1_switch_saveTerms_False_switch_useHyperViscAzVortD_True_tag_CSDXMagFieldScanAr-3-turbulentPhase1_0/restart_1",\
         "CSDXMagFieldScanAr/nout_5000_timestep_1/geom_Lx_7.8633_geom_Ly_275.2144_input_B0_0.1_switch_saveTerms_False_switch_useHyperViscAzVortD_True_tag_CSDXMagFieldScanAr-3-turbulentPhase1_0")
    tLenFirstFolder  = 322
    tLenSecondFolder = 371
    tLenThirdFolder  = 491

    print("\n\nTesting time collects for one folder (no slicing)")
    tcc, tct = getTimes((collectPaths[0],), None)
    verify(tcc, tct)

    print("\n\nTesting time collects for one folder (with slicing)")
    tcc, tct = getTimes((collectPaths[0],), tInd = [3, 5])
    verify(tcc, tct)

    print("\n\nTesting time collects for one folder (using None first)")
    tcc, tct = getTimes((collectPaths[0],), tInd = [None, 5])
    verify(tcc, tct)

    print("\n\nTesting time collects for one folder (using None last)")
    try:
        tcc, tct = getTimes((collectPaths[0],), tInd = [3, None])
        verify(tcc, tct)
    except TypeError as te:
        if "unsupported operand" in te.args[0]:
            print("Collect rasied {}\nExpected, continuing\n\n".format(te.args[0]))
        else:
            raise te

    print("\n\nTesting time collects for two folders (no slicing)")
    tcc, tct = getTimes(collectPaths[:2], tInd = None)
    verify(tcc, tct)

    print("\n\nTesting time collects for two folders (cutting the first folder)")
# FIXME: BUG FOUND HERE
    tcc, tct = getTimes(collectPaths[:2], tInd = [tLenFirstFolder+1, None])
    verify(tcc, tct)

    print("\n\nTesting time collects for two folders (cutting the last folder)")
    tcc, tct = getTimes(collectPaths[:2], tInd = [None, tLenFirstFolder-1])
    verify(tcc, tct)

    print("\n\nTesting time collects for two folders (slice between)")
    tcc, tct = getTimes(collectPaths[:2], tInd = [1, tLenFirstFolder-1])
    verify(tcc, tct)

    print("\n\nTesting time collects for 3 folders (no slicing)")
    tcc, tct = getTimes(collectPaths, tInd = None)
    verify(tcc, tct)

    print("\n\nTesting time collects for three folders (cutting the first folder)")
    tcc, tct = getTimes(collectPaths, tInd = [tLenFirstFolder+1, None])
    verify(tcc, tct)

    print("\n\nTesting time collects for three folders (cutting the two first folders)")
    tcc, tct = getTimes(collectPaths, tInd = [tLenFirstFolder+tLenSecondFolder+1, None])
    verify(tcc, tct)

    print("\n\nTesting time collects for three folders (cutting the last folder)")
    tcc, tct = getTimes(collectPaths, tInd = [None, tLenFirstFolder + tLenSecondFolder-1])
    verify(tcc, tct)

    print("\n\nTesting time collects for three folders (cutting the two last folder)")
    tcc, tct = getTimes(collectPaths, tInd = [None, tLenFirstFolder-1])
    verify(tcc, tct)

    print("\n\nTesting time collects for three folders (slicing the full array)")
    tcc, tct = getTimes(collectPaths, tInd = [1, tLenFirstFolder + tLenSecondFolder + tLenThirdFolder-1])
    verify(tcc, tct)

    print("\n\nTesting time collects for three folders (slicing the two first arrays)")
    tcc, tct = getTimes(collectPaths, tInd = [1, tLenFirstFolder + tLenSecondFolder-1])
    verify(tcc, tct)

    print("\n\nTesting time collects for three folders (slicing the two last arrays)")
    tcc, tct = getTimes(collectPaths, tInd = [tLenFirstFolder + 1, tLenFirstFolder + tLenSecondFolder-1])
    verify(tcc, tct)
#}}}

if __name__ == "__main__":
    timeCollectTest()

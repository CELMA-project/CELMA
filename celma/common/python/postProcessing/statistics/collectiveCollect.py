#!/usr/bin/env python

"""
Contains function which collects a variable over several output timesteps
"""


# Collctive collect
paths = [\
"nout_100_timestep_1/switch_forceAddNoise_True_switch_includeNoise_True_switch_saveTerms_False_switch_useHyperViscAzVortD_True_tag_2-a-0-linearPhase1_0/",\
"nout_100_timestep_1/switch_saveTerms_False_switch_useHyperViscAzVortD_True_tag_2-a-1-linearPhase2_0/",\
"nout_100_timestep_1/switch_saveTerms_False_switch_useHyperViscAzVortD_True_tag_2-a-2-linearPhase3_0/",\
"nout_100_timestep_1/switch_saveTerms_False_switch_useHyperViscAzVortD_True_tag_2-a-3-linearPhase4_0/",\
"nout_100_timestep_1/switch_useHyperViscAzVortD_True_tag_2-a-4-linearPhase5_0/",\
"nout_300_timestep_1/switch_useHyperViscAzVortD_True_tag_3-a-0-turbulentPhase_0/",\
"nout_300_timestep_1/switch_saveTerms_False_switch_useHyperViscAzVortD_True_tag_3-a-1-turbulentPhase2_0/",\
"nout_300_timestep_1/switch_saveTerms_False_switch_useHyperViscAzVortD_True_tag_3-a-2-turbulentPhase3_0/",\
"nout_300_timestep_1/switch_saveTerms_False_switch_useHyperViscAzVortD_True_tag_3-a-3-turbulentPhase4_0/",\
"nout_300_timestep_1/switch_saveTerms_False_switch_useHyperViscAzVortD_True_tag_3-a-4-turbulentPhase5_0/",\
"nout_1000_timestep_1/switch_saveTerms_False_switch_useHyperViscAzVortD_True_tag_3-a-5-turbulentPhase6_0/",\
"nout_1000_timestep_1/switch_saveTerms_False_switch_useHyperViscAzVortD_True_tag_3-a-6-turbulentPhase7_0/",\
"nout_1000_timestep_1/switch_saveTerms_False_switch_useHyperViscAzVortD_True_tag_3-a-7-turbulentPhase8_0/",\
        ]

#{{{collectiveCollect
def collectiveCollect(paths, varStrings,\
                      collectGhost=False,\
                      tInd=None, yInd=None, xInd=None, zInd=None):
    """
    Collects variables from several paths

    Parameters
    ----------
    paths : iterable of strings
        The paths to collect from. Must be in ascending order of the
        simulation time, as the variables are being concatenated
    varStrings : iterable of strings
        The variables to be collected
    collectGhost : bool
        If the ghost is to be collected
    tind : [None|2d array]
        Time index range to collect. The first index is the start, and
        the second is the end of the range (inclusive)
    xind : [None|2d array]
        x index range to collect. The first index is the start, and the
        second is the end of the range (inclusive)
    yind : [None|2d array]
        y index range to collect. The first index is the start, and the
        second is the end of the range (inclusive)
    zind : [None|2d array]
        z index range to collect. The first index is the start, and the
        second is the end of the range (inclusive)

    Return
    ------
    data : dict
        A dictionary of the concatenated variables
    """

    data = {var: np.array([]) for var in varStrings}

    for path in paths:
        for var in varStrings:
            # Make a local var which is reused for every interation,
            # then concatenate the dictionary
            localVar =\
                collect(var,path=path,\
                        tind = tInd,\
                        xind = xInd,\
                        yind = yInd,\
                        zind = zInd,\
                        xguards=guards,\
                        yguards=guards,\
                        info=False)
            data[var] = np.concatenate((data[var], localVar))

    return data
#}}}

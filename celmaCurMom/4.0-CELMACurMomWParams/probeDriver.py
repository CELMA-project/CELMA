#!/usr/bin/env python

"""Driver which post-process using probes."""

import os, sys
import pickle
# If we add to sys.path, then it must be an absolute path
commonDir = os.path.abspath('./../common')
# Sys path is a list of system paths
sys.path.append(commonDir)
from CELMAPython.statsAndSignals import PerpPlaneProbes

# paths = [\
# "a-data/nout_100_timestep_1/switch_forceAddNoise_True_switch_includeNoise_True_switch_saveTerms_False_switch_useHyperViscAzVortD_True_tag_2-a-0-linearPhase1_0/",\
# "a-data/nout_100_timestep_1/switch_saveTerms_False_switch_useHyperViscAzVortD_True_tag_2-a-1-linearPhase2_0/",\
# "a-data/nout_100_timestep_1/switch_saveTerms_False_switch_useHyperViscAzVortD_True_tag_2-a-2-linearPhase3_0/",\
# "a-data/nout_100_timestep_1/switch_saveTerms_False_switch_useHyperViscAzVortD_True_tag_2-a-3-linearPhase4_0/",\
# "a-data/nout_100_timestep_1/switch_useHyperViscAzVortD_True_tag_2-a-4-linearPhase5_0/",\
# "a-data/nout_300_timestep_1/switch_useHyperViscAzVortD_True_tag_3-a-0-turbulentPhase_0/",\
# "a-data/nout_300_timestep_1/switch_saveTerms_False_switch_useHyperViscAzVortD_True_tag_3-a-1-turbulentPhase2_0/",\
# "a-data/nout_300_timestep_1/switch_saveTerms_False_switch_useHyperViscAzVortD_True_tag_3-a-2-turbulentPhase3_0/",\
# "a-data/nout_300_timestep_1/switch_saveTerms_False_switch_useHyperViscAzVortD_True_tag_3-a-3-turbulentPhase4_0/",\
# "a-data/nout_300_timestep_1/switch_saveTerms_False_switch_useHyperViscAzVortD_True_tag_3-a-4-turbulentPhase5_0/",\
# "a-data/nout_1000_timestep_1/switch_saveTerms_False_switch_useHyperViscAzVortD_True_tag_3-a-5-turbulentPhase6_0/",\
# "a-data/nout_1000_timestep_1/switch_saveTerms_False_switch_useHyperViscAzVortD_True_tag_3-a-6-turbulentPhase7_0/",\
# "a-data/nout_1000_timestep_1/switch_saveTerms_False_switch_useHyperViscAzVortD_True_tag_3-a-7-turbulentPhase8_0/",\
# "a-data/nout_1000_timestep_1/switch_saveTerms_False_switch_useHyperViscAzVortD_True_tag_3-a-8-turbulentPhase9_0/",\
#         ]

# Get all the paths
# paths = [\
#  "newKiwiFlat/nout_500_timestep_1/switch_forceAddNoise_True_switch_includeNoise_True_switch_saveTerms_False_switch_useHyperViscAzVortD_True_tag_a1-KiwiFlat-2-linearPhase1_0/",\
#  "newKiwiFlat/nout_5000_timestep_1/switch_saveTerms_False_switch_useHyperViscAzVortD_True_tag_a1-KiwiFlat-3-turbulentPhase1_0/",\
#         ]

paths = [\
 "newKiwiFlat/nout_500_timestep_1/switch_forceAddNoise_True_switch_includeNoise_True_switch_saveTerms_False_switch_useHyperViscAzVortD_True_tag_a1-KiwiFlat-2-linearPhase1_0/",\
 "newKiwiFlat/nout_5000_timestep_1/switch_saveTerms_False_switch_useHyperViscAzVortD_True_tag_a1-KiwiFlat-3-turbulentPhase1_0/",\
        ]

# Create the probe
ppp = PerpPlaneProbes('n',\
                      paths=paths,\
                      yInd=16,\
                      nProbes=5,\
                      physicalUnits=False,\
                      steadyStatePath="newKiwiFlat/nout_20_timestep_50/nz_256/ownFilters_type_none_tag_a1-KiwiFlat-1-expand_0/",\
                      radialProbeIndices=None)

# Create the probe
ppp.initializeInputOutput(ppp.radialProbesIndices, [ppp.yInd], [0])
ppp.calcStatsMoments()
ppp.calcPDFs()
ppp.calcPSDs()
ppp.calcAvgFluxThroughVolumeElement(ppp.radialExB, "ExB")
ppp.calcFFTs()

# Pickle the result
with open('n.pickle', 'wb') as f:
    # Pickle the 'data' dictionary using the highest protocol available.
    pickle.dump(ppp, f, pickle.HIGHEST_PROTOCOL)

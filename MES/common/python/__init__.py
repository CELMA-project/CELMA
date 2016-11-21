"""Makes a shortcut to CELMAPython.plotHelpers rcParams"""

import os, sys
# If we add to sys.path, then it must be an absolute path
commonDir = os.path.abspath("./../../../celma/common")
# Sys path is a list of system paths
sys.path.append(commonDir)

from CELMAPython.plotHelpers import __init__

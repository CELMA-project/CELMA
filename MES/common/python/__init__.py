"""Makes a shortcut to CELMAPython.plotHelpers rcParams"""

import os, sys

# Find the root of the CELMA folder
path = os.path.dirname(os.path.abspath(__file__))

# Infinite loop which runs until break
while True:
    base = os.path.basename(path)
    if base != "MES":
        path = os.path.dirname(path)
    else:
        # Step one back and join celma, common
        path      = os.path.dirname(path)
        commonDir = os.path.join(path, "celma", "common")
        break

# Sys path is a list of system paths
sys.path.append(commonDir)

from CELMAPython.plotHelpers import __init__

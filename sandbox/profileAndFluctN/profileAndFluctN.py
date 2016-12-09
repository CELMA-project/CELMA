#%%
import os, sys
sys.path.append("/home/mmag/Desktop/CELMA-dev/celma/common/")
from CELMAPython.statsAndSignals.averages import timeAvg

import pickle
import numpy as np
with open("nt", "rb") as f:
   n, t = pickle.load(f)

tAvgN, tAvgT = timeAvg(n,t)

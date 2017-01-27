#!/usr/bin/env python

"""
Script which edits the dmpFolders.

Nifty when not all jobs have finished
"""

import pickle

with open("dmpFoldersDict.pickle", "rb") as f:
    a = pickle.load(f)

a["extraTurbulence"] = list(a["extraTurbulence"])
for nr in range(len(a["extraTurbulence"])):
    a["extraTurbulence"][nr] = list(a["extraTurbulence"][nr])
    a["extraTurbulence"][nr].pop(-1)
    a["extraTurbulence"][nr] = tuple(a["extraTurbulence"][nr])
a["extraTurbulence"] = tuple(a["extraTurbulence"])

a["turbulence"] = list(a["turbulence"])
for nr in range(len(a["turbulence"])):
    a["turbulence"][nr] += "/restart_3"
    print(a["turbulence"][nr])
a["turbulence"] = tuple(a["turbulence"])

with open("dmpFoldersDict.pickle", "wb") as f:
    pickle.dump(a,f)

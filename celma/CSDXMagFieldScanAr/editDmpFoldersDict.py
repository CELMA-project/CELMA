#!/usr/bin/env python

"""
Script which edits the dmpFolders.

Nifty when not all jobs have finished
"""

import pickle

with open("dmpFoldersDict.pickle", "rb") as f:
    dmpFolders = pickle.load(f)

dmpFolders["extraTurbulence"] = list(dmpFolders["extraTurbulence"])
for nr in range(len(dmpFolders["extraTurbulence"])):
    dmpFolders["extraTurbulence"][nr] = list(dmpFolders["extraTurbulence"][nr])
    # Remove the current running folder
    dmpFolders["extraTurbulence"][nr].pop(-1)
    dmpFolders["extraTurbulence"][nr] = tuple(dmpFolders["extraTurbulence"][nr])
dmpFolders["extraTurbulence"] = tuple(dmpFolders["extraTurbulence"])

dmpFolders["turbulence"] = list(dmpFolders["turbulence"])
for nr in range(len(dmpFolders["turbulence"])):
    # Set the last finished folder to the last finished folder
    dmpFolders["turbulence"][nr] += "/restart_10"
dmpFolders["turbulence"] = tuple(dmpFolders["turbulence"])

with open("dmpFoldersDict.pickle", "wb") as f:
    pickle.dump(dmpFolders,f)

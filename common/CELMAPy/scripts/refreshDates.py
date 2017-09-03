#!/usr/bin/env python

"""Refresh dates of files to prevent automatic deletion  by cluster"""

import pathlib

allFiles = list(pathlib.Path(".").glob("**/*"))
for file in allFiles:
    try:
        with open(file, "r") as f:
            lines = f.readlines()
        with open(file, "w") as f:
            for l in lines:
                f.write(l)
    except UnicodeDecodeError:
        with open(file, "rb") as f:
            lines = f.readlines()
        with open(file, "wb") as f:
            for l in lines:
                f.write(l)

#!/usr/bin/env python

"""Refresh dates of files to prevent automatic deletion  by cluster"""

import pathlib

def refreshDates():
    """Opens and saves all files recursively"""

    allFiles = list(pathlib.Path(".").glob("**/*"))
    for path in allFiles:
        if not path.is_file():
            continue
        path = str(path)
        try:
            with open(path, "r") as f:
                lines = f.readlines()
            with open(path, "w") as f:
                for l in lines:
                    f.write(l)
        except UnicodeDecodeError:
            with open(path, "rb") as f:
                lines = f.readlines()
            with open(path, "wb") as f:
                for l in lines:
                    f.write(l)

if __name__ == "__main__":
    refreshDates()

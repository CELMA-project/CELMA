"""Refresh dates of files to prevent automatic deletion by cluster"""

import pathlib
from os import utime

def refreshDates():
    """Opens and saves all files recursively"""

    allFiles = list(pathlib.Path(".").glob("**/*"))
    for path in allFiles:
        if not path.is_file():
            continue
        path = str(path)
        try:
            utime(path)
        except PermissionError:
            print("Permission denied for {}".format(path))

if __name__ == "__main__":
    refreshDates()

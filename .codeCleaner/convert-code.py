#!/usr/bin/env python

"""
Contains a function to convert the code
"""

from glob import glob
from subprocess import call
from os.path import expanduser

def convert(path):
    """
    Converts the code using the binary specified in the input path
    """

    extension = ("*.cxx", "*.hxx")
    for ex in extension:
        files = glob("./../**/"+ex, recursive = True)
        for f in files:
            cmd = "python {} -r {}".format(path, f)
            print(cmd)
            call(cmd.split())

if __name__ == "__main__":
    home = expanduser("~")
    path = "{}/BOUT-dev/bin/bout_3to4.py".format(home)
    convert(path)

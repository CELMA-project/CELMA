#!/usr/bin/env python

"""
Script which formats all .cxx and .hxx according to the .clang-format
style (based upon the one found in BOUT-dev).
"""

from glob import glob
from subprocess import call
import shutil

# Check if clang-format is installed
path = shutil.which("clang-format")
if path is None:
    raise RuntimeError("This script requires 'clang-format'.")

extension = ("*.cxx", "*.hxx")
for ex in extension:
    files = glob("./../**/"+ex, recursive = True)
    for f in files:
        cmd = "clang-format -style=file {}".format(f)
        print(cmd)
        with open("tmp", "w") as tmp:
            call(cmd.split(), stdout=tmp)
            shutil.move("tmp", f)

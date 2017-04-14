#!/usr/bin/env python

"""
This script removes all outputs from jupyter notebooks

NOTE:
Must be run from the top level of working three

To clean the repository:

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! WARNING:                 !
!    This rewrites history !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

git filter-branch --tree-filter 'python <absPathToFilter>/.filters/clear_all_ipynb.py' HEAD --all
"""

import os
import shutil
from glob import glob
from subprocess import call

# Change to the file's path
path = os.path.realpath(__file__)
os.chdir(os.path.dirname(path))
# Find all notebooks
notebooks = glob("./**/*.ipynb", recursive=True)

for nb in notebooks:
    nb = os.path.realpath(nb)
    new = nb+".new"
    call("{}/rm_outputs_jupyter_nb < {} > {}".format(path, nb, new), shell=True)
    os.remove(nb)
    shutil.move(new, nb)

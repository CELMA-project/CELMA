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
from glob import glob
from rm_outputs_jupyter_nb import rm_outputs

# Change to the file's super path
path = os.path.realpath(__file__)
os.chdir(os.path.dirname(path))
os.chdir("..")
# Find all notebooks
notebooks = glob("./**/*.ipynb", recursive=True)

for nb in notebooks:
    rm_outputs(nb)

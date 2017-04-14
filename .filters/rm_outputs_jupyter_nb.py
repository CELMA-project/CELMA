#!/usr/bin/env python

"""
Suppress output and prompt numbers in git version control.

Modified from
https://gist.github.com/pbugnion/ea2797393033b54674af

This script will tell git to ignore prompt numbers and cell output.

See also this blogpost: http://pascalbugnion.net/blog/ipython-notebooks-and-git.html.
"""

import sys
import json

def rm_outputs(filename):
    """
    Removes the outputs by overwriting the file.
    """

    with open(filename, "r") as nb:
        json_in = json.load(nb)

    ipy_version = int(json_in["nbformat"])-1 # nbformat is 1 more than actual version.

    def strip_output_from_cell(cell):
        if "outputs" in cell:
            cell["outputs"] = []
        if "prompt_number" in cell:
            del cell["prompt_number"]

    if ipy_version == 2:
        for sheet in json_in["worksheets"]:
            for cell in sheet["cells"]:
                strip_output_from_cell(cell)
    else:
        for cell in json_in["cells"]:
            strip_output_from_cell(cell)

    with open(filename, "w") as nb:
        json.dump(json_in, nb, sort_keys=True, indent=1, separators=(",",": "))

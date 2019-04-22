#!/usr/bin/env python

"""
Strip outputs from an IPython Notebook

Opens a notebook, strips its output, and writes the outputless version to the original file.
Useful mainly as a git filter or pre-commit hook for users who don"t want to track output in VCS.
This does mostly the same thing as the `Clear All Output` command in the notebook UI.
LICENSE: Public Domain

Source: https://gist.github.com/minrk/6176788
"""

import io
import sys
import re
import fileinput

try:
    # Jupyter >= 4
    from nbformat import read, write, NO_CONVERT
except ImportError:
    # IPython 3
    try:
        from IPython.nbformat import read, write, NO_CONVERT
    except ImportError:
        # IPython < 3
        from IPython.nbformat import current

        def read(f, as_version):
            return current.read(f, "json")

        def write(nb, f):
            return current.write(nb, f, "json")


def _cells(nb):
    """Yield all cells in an nbformat-insensitive manner"""
    if nb.nbformat < 4:
        for ws in nb.worksheets:
            for cell in ws.cells:
                yield cell
    else:
        for cell in nb.cells:
            yield cell


def strip_output(nb):
    """Strips the outputs from a notebook object"""
    nb.metadata.pop("signature", None)
    for cell in _cells(nb):
        if "outputs" in cell:
            cell["outputs"] = []
        if "prompt_number" in cell:
            cell["prompt_number"] = None
    return nb


def repair_corrupted(filename):
    """Repairs corrupted notebooks with unmatched double quotes"""
    print("Repairing corrupted Notebook...", end="")

    pattern = re.compile(r"^[^\"]*\"[^\"]*$")

    with fileinput.FileInput(filename, inplace=True, backup=None) as f:
        for line in f:
            match = re.search(pattern, line)
            if match is None:
                # Last character is a newline
                print(line, end="")
            else:
                # Add the missing double quote (last character is a newline)
                print(match.group(0)[:-1] + "\"")
    print("done")


if __name__ == "__main__":
    # In case the filename contains spaces, we need all sys.argv
    # (excluding the filename)
    filename = " ".join(sys.argv[1:])
    try:
        with io.open(filename, 'r', encoding='utf8') as f:
            nb = read(f, as_version=NO_CONVERT)
    except Exception as e:
        if "Notebook does not appear to be JSON" in e.args[0]:
            repair_corrupted(filename)
            with io.open(filename, 'r', encoding='utf8') as f:
                nb = read(f, as_version=NO_CONVERT)
        else:
            raise e
    nb = strip_output(nb)
    with io.open(filename, "w", encoding="utf8") as f:
        write(nb, f)

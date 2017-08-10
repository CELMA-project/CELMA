#!/usr/bin/env python

"""
Saves all *.log.* and *.restart.* files of a directory to a zip.

The folder structure will be preserved.
The zip will be mailed if a mail is set.
"""

import os, shutil, pathlib, subprocess

# INPUT
# =============================================================================
# Directory to find log and restart files from
directory = "CSDXMagFieldScanAr"
# If you would like a mail on finished job enter your mail here
# Example: "john@doe.com"
mail = None
# =============================================================================

logRestartDir = "logAndRestartFilesFrom" + directory

# Find all the log files
logFiles = list(pathlib.Path(directory).glob("**/*.log.*"))

# Find all the restart files
restartFiles = list(pathlib.Path(directory).glob("**/*.restart.*"))

# Combine
allFiles = (*logFiles, *restartFiles)

for f in allFiles:
    # Add logRestartDir to directory without its root
    dst = pathlib.PurePath(logRestartDir, f.relative_to(*f.parts[:2]))
    # Clean folders marked with BAK
    if "BAK" not in str(dst):
        # Make the paths
        newPath = pathlib.Path(*dst.parts[:-1])
        if not os.path.exists(str(newPath)):
            newPath.mkdir(parents = True)

        shutil.copy2(str(f), str(dst))

# Compress and remove the folder
shutil.make_archive(logRestartDir, "zip", logRestartDir)
shutil.rmtree(logRestartDir)


if mail:
    # Sends mail through the terminal
    theZip = logRestartDir + ".zip"
    cmd = (\
        'echo "See attachment" | mutt -a "{}" -s "Log and restart files" -- {}'\
          ).format(theZip, mail)

    process = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()
    print("{} sent to {}".format(theZip, mail))

    # Clean-up
    shutil.rm(theZip)

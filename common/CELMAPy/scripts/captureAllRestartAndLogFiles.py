#!/usr/bin/env python

"""
Contains the captureAllRestartAndLogFiles function.
"""

import os, shutil, pathlib, subprocess

def captureAllRestartAndLogFiles(directory, mail=None):
    """
    Saves all *.log.* and *.restart.* files of a directory to a zip.

    The folder structure will be preserved.
    The zip will be mailed if a mail is set.

    Parameters
    ----------
    directory : str
        Directory to recursively find log and restart files from.
    mail : [None | str]
        If set to a string, the zip will be mailed if mutt mailing is
        supported.
    """

    if directory == ".":
        logRestartDir = "logAndRestartFilesFromRoot"
    else:
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
            pathlib.Path(*dst.parts[:-1]).mkdir(parents = True, exist_ok=True)
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

if __name__ == "__main__":
    captureAllRestartAndLogFiles(".")

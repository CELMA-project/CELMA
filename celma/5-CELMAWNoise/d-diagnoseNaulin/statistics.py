#!/usr/bin/env python

"""Driver which finds max, min, avg, std-dev of the Naulin Solver"""

import os, glob
import numpy as np

superFolders = [ "nout_1_timestep_1/ny_24", "nout_1_timestep_10/ny_24" ]

for d in superFolders:
    # Find all subfolders
    subDirs = [os.path.join(d,el) for el in os.listdir(d) if os.path.isdir(os.path.join(d,el))]
    for sd in subDirs:
        print("Checking in {}".format(sd))
        # Find the log file of stdout (0-*.log)
        searchFor = "0-*.log"
        hit = glob.glob(os.path.join(sd, searchFor))[0]

        # Read the file line for line, and store it into a list
        with open(hit) as f:
            lines = f.readlines()

        # Filter the lines where the iteration count is given
        solverLine = [el for el in lines if "NaulinSolver success" in el]

        # List to store iteration count in
        it = []
        for line in solverLine:
            # A typical line looks like
            # NaulinSolver success: Iterations=3 abserr=1.9742e-10 relerr=2.79877e-06
            it.append(int(line.split("=")[1].split(" ")[0]))

        # Make a np array
        it = np.array(it)

        print("Min - {}, max - {}, mean - {}, std - {}".\
                format(it.min(), it.max(), it.mean(), it.std()))

#!/usr/bin/env python

"""Function which generates equilibrium profiles for the rewritten equations."""

import glob, os
from boututils.datafile import DataFile
import numpy as np

def createRestart(path="data", output="./", informat="nc", outformat=None):
    """
    Rewrites the restart files in order to be compatible with the new
    set of equations

    Returns True on success
    """

    if outformat == None:
        outformat = informat

    if path == output:
        raise ValuError("Can't overwrite restart files when reforming")

    file_list = glob.glob(os.path.join(path, "BOUT.restart.*."+informat))
    file_list.sort()
    nfiles = len(file_list)

    if nfiles == 0:
        print("ERROR: No data found")
        return False

    print("Number of files found: " + str(nfiles))

    for f in file_list:
        new_f = os.path.join(output, f.split('/')[-1])
        print("Changing {} => {}".format(f, new_f))

        # Open the restart file in read mode and create the new file
        with DataFile(f) as old,\
             DataFile(new_f, write=True, create=True) as new:

            # 3D-vars in old field
            oldFields = {}

            # Loop over the variables in the old file
            for var in old.list():
                # Read the data
                data = old.read(var)

                # Find 3D variables
                if old.ndims(var) == 3:
                    oldFields[var] = data
                else:
                    print("    Copying "+var)
                    new.write(var, data)

            print("    Copying lnN")
            new.write("lnN", oldFields['lnN'])
            print("    Copying vortD")
            new.write("vortD", oldFields['vortD'])
            print("Generating new variables")
            jPar       = np.exp(oldFields['lnN'])*\
                         (oldFields['uIPar'] - oldFields['uEPar'])
            momDensPar = np.exp(oldFields['lnN'])*oldFields['uIPar']
            print("    Writing jPar")
            new.write("jPar"      , jPar      )
            print("    Writing momDensPar")
            new.write("momDensPar", momDensPar)

    print("\n"*2+"Restart files successfully created")
    return True

if __name__ == "__main__":
    path = "./../../celma/9-CELMAUpdateArtViscCleanUp/b-longerNeumannO4/nout_20_timestep_101/nz_1/cst_artViscParVortD_0.0_cst_artViscPerpVortD_0.0_ownFilters_type_none_switch_useHyperViscAzVortD_False_tag_longNeumannO4Initializer_0/"
    output = "equilibriumRestart"
    createRestart(path, output)

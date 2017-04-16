#!/usr/bin/env bash

# Installs the master branch of BOUT-dev


# Options
# ==============================================================================
VERBOSE="false"          # Only set to "false" in order for travis to work
INSTALL_CONDA="true"     # Needed for post-processing
INSTALL_CMAKE="false"    # Needed for sundials if CMAKE is below 2.8.11
INSTALL_FFMPEG="true"    # Needed for post-processing if x264 is not present
INCL_SUNDIALS="true"     # The preferred time solver
INCL_PETSC_SLEPC="false" # Only needed for fancy features
OPTIMIZING="true"        # Good for speed
DEBUG="false"            # Good for debugging, bad for speed
# ==============================================================================


# Preparations
# ==============================================================================
# exit on error
set -e

# Get the path of the calling script
# http://stackoverflow.com/questions/59895/getting-the-source-directory-of-a-bash-script-from-within
CURDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Ensure paths are available when building
export PATH="$HOME/local/bin:$PATH"
export LD_LIBRARY_PATH=$HOME/local/lib:$LD_LIBRARY_PATH

# Extra packages and flags for BOUT-dev
EXTRA_PACKAGES=""
EXTRA_FLAGS=""

# Set appropriate flags
if [ "$OPTIMIZING" = "true" ]; then
    EXTRA_FLAGS="${EXTRA_FLAGS} --enable-checks=no --enable-optimize=3"
elif [ "$DEBUG" = "true" ]; then
    EXTRA_FLAGS="${EXTRA_FLAGS} --enable-debug --enable-checks=3 --enable-track"
fi

# Travis 4MB workaround
# http://stackoverflow.com/questions/26082444/how-to-work-around-travis-cis-4mb-output-limit/26082445#26082445
# ..............................................................................
export PING_SLEEP=30s
export CURDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
export BUILD_OUTPUT=$CURDIR/build.out

touch $BUILD_OUTPUT

dump_output() {
   echo Tailing the last 500 lines of output:
   tail -500 $BUILD_OUTPUT
}
error_handler() {
  kill $PING_LOOP_PID
  echo ERROR: An error was encountered with the build.
  dump_output
  exit 1
}

# If an error occurs, run our error handler to output a tail of the build
trap 'error_handler' ERR

# Set up a repeating loop to send some output to Travis.
bash -c "while true; do echo \$(date) - building ...; sleep $PING_SLEEP; done" &
PING_LOOP_PID=$!
# ..............................................................................

# Set the verbosity
if [ "$VERBOSE" = "true" ]; then
    # Could also have used " > /dev/stdout"
    STDOUT=""
else
    STDOUT=" &>> $BUILD_OUTPUT"
fi
# ==============================================================================


# Install packages needed for BOUT-dev
# ==============================================================================
if [ "$INSTALL_CONDA" = "true" ]; then
    echo "bash $CURDIR/condaInstall.sh $STDOUT"
    bash $CURDIR/condaInstall.sh $STDOUT
fi

if [ "$INSTALL_CMAKE" = "true" ]; then
    echo "bash $CURDIR/cmakeInstall.sh $STDOUT"
    bash $CURDIR/cmakeInstall.sh $STDOUT
fi

if [ "$INSTALL_FFMPEG" = "true" ]; then
    echo "bash $CURDIR/ffmpegInstall.sh $STDOUT"
    bash $CURDIR/ffmpegInstall.sh $STDOUT
fi

# Install mpi
echo "bash $CURDIR/mpiInstall.sh $STDOUT"
bash $CURDIR/mpiInstall.sh $STDOUT

# Install fftw
echo "bash $CURDIR/fftwInstall.sh $STDOUT"
bash $CURDIR/fftwInstall.sh $STDOUT

# Install hdf5
echo "bash $CURDIR/hdf5Install.sh $STDOUT"
bash $CURDIR/hdf5Install.sh $STDOUT

# Install netcdf
echo "bash $CURDIR/netcdfInstall.sh $STDOUT"
bash $CURDIR/netcdfInstall.sh $STDOUT

if [ "$INCL_SUNDIALS" = "true" ]; then
    # Install sudials
    echo "bash $CURDIR/sundialsInstall.sh $STDOUT"
    bash $CURDIR/sundialsInstall.sh $STDOUT
    EXTRA_PACKAGES="${EXTRA_PACKAGES} --with-sundials"
fi

if [ "$INCL_PETSC_SLEPC" = true ]; then
    # Install PETSc
    echo "bash $CURDIR/PETScInstall.sh $STDOUT"
    bash $CURDIR/PETScInstall.sh $STDOUT

    # Install SLEPc
    echo "bash $CURDIR/SLEPcInstall.sh $STDOUT"
    bash $CURDIR/SLEPcInstall.sh $STDOUT
    EXTRA_PACKAGES="${EXTRA_PACKAGES} --with-petsc --with-slepc"
fi
# ==============================================================================


# Builiding BOUT-dev master
# ==============================================================================
echo -e "\n\n\nInstalling BOUT-dev\n\n\n"
cd $HOME
git clone https://github.com/boutproject/BOUT-dev.git
cd BOUT-dev
# NOTE: Explicilty state netcdf and hdf5 in order not to mix with anaconda
./configure ${EXTRA_FLAGS} ${EXTRA_PACKAGES} --with-netcdf=$HOME/local --with-hdf5=$HOME/local
make
echo -e "\n\n\nDone installing BOUT-dev\n\n\n"
echo -e "\n\n\nChecking installation\n\n\n"
cd examples/bout_runners_example
make
python 6a-run_with_MMS_post_processing_specify_numbers.py
# ==============================================================================


# Write what needs to be exported
# ==============================================================================
echo -e "\n\n\nInstallation complete.\n"
echo -e "Make sure the following lines are present in your .bashrc:"
echo -e "export PYTHONPATH=\"\$HOME/BOUT-dev/tools/pylib/:\$PYTHONPATH\""
echo -e "export PATH=\"\$HOME/local/bin:\$PATH\""
echo -e "export LD_LIBRARY_PATH=\"\$HOME/local/lib:\$LD_LIBRARY_PATH\""

if [ "$INCL_PETSC_SLEPC" = true ]; then
    echo -e "export PETSC_DIR=\"\$HOME/petsc-\${PETSC_VERSION}\""
    echo -e "export PETSC_ARCH=\"arch-linux2-cxx-debug\""
    echo -e "export SLEPC_DIR=\"\$HOME/slepc-\${SLEPC_VERSION}\""
fi

if [ "$INSTALL_CONDA" = true ]; then
    echo -e "export PATH=\"\$HOME/anaconda3/bin:\$PATH\""
fi
# ==============================================================================


# Termination of travis workaround
# ==============================================================================
# The build finished without returning an error so dump a tail of the output
if [ ! "$VERBOSE" = "true" ]; then
    dump_output
fi

# nicely terminate the ping output loop
kill $PING_LOOP_PID
# ==============================================================================

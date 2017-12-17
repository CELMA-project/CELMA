#!/usr/bin/env bash

# NOTE: No workaround for the new hdf5-1.10.1 has been found.
#       Apparent download location:
#       https://www.hdfgroup.org/downloads/hdf5/source-code/#
#       Or login via
#       https://www.hdfgroup.org/downloads/hdf5/

# Installs hdf5
HDF5_MAJOR="1"
HDF5_MINOR="8"
HDF5_PATCH="20"

# exit on error
set -e

HDF5_MAJOR_MINOR=${HDF5_MAJOR}.${HDF5_MINOR}
HDF5_VERSION=${HDF5_MAJOR_MINOR}.${HDF5_PATCH}

# Install hdf5
echo -e "\n\n\nInstalling hdf5\n\n\n"
cd $HOME
mkdir -p local
cd local
mkdir -p examples
cd ..
mkdir -p install
cd install
wget https://support.hdfgroup.org/ftp/HDF5/current18/bin/hdf5-${HDF5_VERSION}-linux-centos7-x86_64-gcc485-shared.tar.gz -O hdf5-${HDF5_VERSION}.tar.gz
tar -xzvf hdf5-${HDF5_VERSION}.tar.gz
cd hdf5-${HDF5_VERSION}
# : is the no-op command
make clean || :
./configure --prefix=$HOME/local --enable-cxx=yes
make
make install
echo -e "\n\n\nDone installing hdf5\n\n\n"


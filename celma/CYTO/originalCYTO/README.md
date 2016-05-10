# CYTO
This is a non-bloated version of 'https://svn.fysik.dtu.dk/vona/check/TYR',
originally checked out at r42

## Install
Disclaimer: Additional dependencies could occur, but is not found since
`BOUT-dev` is already installed.

```
sudo apt-get install libjpeg-dev libfftw3-dev
mkdir install
cd install
wget https://www.hdfgroup.org/ftp/HDF/releases/HDF4.2.7/src/hdf-4.2.7-patch1.tar.gz
tar -xzvf hdf-4.2.7-patch1.tar.gz
cd hdf-4.2.7-patch1
./configure --prefix=$PWD/../..
make
make install
cd ..
wget https://www.open-mpi.org/software/ompi/v1.10/downloads/openmpi-1.10.2.tar.gz
tar -xzvf openmpi-1.10.2.tar.gz
cd openmpi-1.10.2
./configure --prefix=$PWD/../..
make
# Grab coffee, as this can take a while
make install
cd ../../cyto
```

## Run

```
rm -f ini.00* HDF_DS_p0_l* cyto.0* cyto.1* cyto.2* cyto.3* database/cyto && make clean && make DEBUG=1 && mv cyto.ubuntu cyto
./cyto -D
```

**Note**: Other ways of running it seems to raise a NaN number

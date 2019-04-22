# install

Contains scripts for locally installing the `master`-branch of `BOUT-dev`
together with its dependencies.
`BOUT-dev` will be installed in `$HOME`, whereas the dependencies will be
installed in `local`.

For a costume install of `BOUT++`, see the [`BOUT++`
user-manual](http://bout-dev.readthedocs.io/en/latest/)

* [boutPPInstall.sh](boutPPInstall.sh) - Builinding and installing the
  `master`-branch of `BOUT-dev` together with its dependencies.
  The installation can be modified by altering the options flags in the script.
  In particular, it is recommended to set `VERBOSE="true"`
* [cmakeInstall.sh](cmakeInstall.sh) - Builinding and installing `CMAKE`
* [condaInstall.sh](condaInstall.sh) - Builinding and installing the `python`
  package manager
* [ffmpegInstall.sh](ffmpegInstall.sh) - Builinding and installing `yasm`, `x264` and `ffmpeg`
* [fftwInstall.sh](fftwInstall.sh) - Builinding and installing `fftw`
* [gccLocalInstall.sh](gccLocalInstall.sh) - Builinding and installing `gcc`.
  Not used in `boutPPinstall.sh`.
  Only included as some clusters have archaic compilers.
  Note that these clusters usually have up-to-date compilers available through
  `modules`.
  Based on http://luiarthur.github.io/gccinstall
* [hdf5Install.sh](hdf5Install.sh) - Builinding and installing `hdf5`
* [mpiInstall.sh](mpiInstall.sh) - Builinding and installing `mpi`
* [netcdfInstall.sh](netcdfInstall.sh) - Builinding and installing `netcdf`
* [PETScInstall.sh](PETScInstall.sh) - Builinding and installing `PETSc`
* [SLEPcInstall.sh](SLEPcInstall.sh) - Builinding and installing `SLEPc`
* [sundialsInstall.sh](sundialsInstall.sh) - Builinding and installing `sundials`

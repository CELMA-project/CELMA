# Use ubuntu as the base layer
FROM ubuntu:18.04

# Install apt dependencies
RUN apt-get -yqq update
RUN apt-get -yqq upgrade
# Python dependencies (python2.7 is required for installing PETSc)
RUN apt-get -yqq install python2.7 python3 python3-pip
# Needed to make python2 a runnable command
RUN ln /usr/bin/python2.7 /usr/bin/python2
# Needed for bout-install dependencies
RUN apt-get -yqq install git
RUN apt-get -yqq install build-essential cmake gfortran m4 zlib1g-dev
# Needed for mayavi
RUN apt-get -yqq install locales libsm6 libxext6 libxrender-dev libxt6
RUN apt-get -yqq install libgl1-mesa-glx libfontconfig
RUN apt-get -yqq install libxkbcommon-x11-0
ENV QT_QPA_PLATFORM offscreen
RUN locale-gen en_US.UTF-8
ENV LANG en_US.UTF-8
ENV LANGUAGE en_US:en
ENV LC_ALL en_US.UTF-8

# Install pip dependecies
RUN pip3 install bout-install boututils boutdata

# Install BOUT-dev
RUN bout_install -a

# Dependencies which requires packages installed by bout_install
RUN pip3 install h5py netcdf4
# Nice to have packages
RUN pip3 install pytest ipython python-dotenv

# Set environment variables
ENV PATH="/root/local/bin:${PATH}"
ENV LD_LIBRARY_PATH="/root/local/lib:${LD_LIBRARY_PATH}"

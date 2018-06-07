FROM ubuntu:latest

# Install Ubuntu 18.04.01 LTS (Bionic Beaver)
# Install python 2.7.12 and 3.6.5, and pip 2 and 3

RUN apt-get update && apt-get install -y \
    python2 \
    python3 \
    python2-pip \
    python3-pip \
    libhdf5-serial-dev \
    hdf5-tools 


# Configure HDF5 include and lib directories
ENV C_INCLUDE_PATH /usr/include/hdf5/serial
ENV LD_LIBRARY_PATH /usr/lib/x86_64-linux-gnu/hdf5/serial
ENV HDF5_DIR /usr/lib/x86_64-linux-gnu/hdf5/serial


# Update pip
# Install mercurial for hgtools (needed to install for testing)
RUN python2 -m pip install --upgrade \ 
    pip \
    setuptools \
    mercurial \

RUN python3 -m pip install --upgrade \
    pip3 \
    setuptools \
    mercurial \



FROM ubuntu:latest

# Install Ubuntu 18.04.01 LTS (Bionic Beaver)
# Install python 2.7.12 and 3.6.5 and pip 2 and 3
# Install HDF5
# Install curl (for retrieving GMTK)
RUN apt-get update && apt-get install -y \
    python \
    python3 \
    python-pip \
    libhdf5-serial-dev \
    hdf5-tools \
    curl \
&& rm -rf /var/lib/apt/lists/*
# RUN apt-get install -y libbz2-dev liblzo2-dev libblosc-dev libatlas-dev liblapack-dev

# Configure HDF5 include and lib directories
ENV C_INCLUDE_PATH /usr/include/hdf5/serial
ENV LD_LIBRARY_PATH /usr/lib/x86_64-linux-gnu/hdf5/serial
ENV HDF5_DIR /usr/lib/x86_64-linux-gnu/hdf5/serial


# Install GMTK
# Downloads location: http://melodi.ee.washington.edu/cgi-bin/gmtk-download.pl
# Ubuntu 18.04 defaults to a more recent c++ version; need to specify c++11 with CXXFLAGS
# TODO: Check to see if the image would be smaller if using ADD of the URL,
# skipping install of curl and then using RUN for the remaining gmtk commands
RUN curl -SL http://melodi.ee.washington.edu/downloads/gmtk/gmtk-1.4.4.tar.gz \
    | tar -xz \
    && cd gmtk-1.4.4 \
    && ./configure CXXFLAGS='-std=c++11'\
    && make \
    && make install \
    && cd .. \
    && rm -rf gmtk-1.4.4

# Update pip
# Install mercurial for hgtools (needed to install for testing)
RUN pip install --upgrade \ 
    pip \
    setuptools \
    mercurial \

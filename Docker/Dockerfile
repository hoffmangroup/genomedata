FROM ubuntu:latest

# Install Ubuntu 18.04.01 LTS (Bionic Beaver)
# Install pkg-config for locationg development files (hdf5)
# Install python 2.7 and pip
# Install python 3.6 and pip3
# Install HDF5
# Install Mercurial for hgtools (needed to install for testing)
RUN apt-get update && apt-get install -y \
    pkg-config \
    python \
    python3 \
    python-pip \
    python3-pip \
    libhdf5-dev \
    hdf5-tools \
    mercurial \
&& rm -rf /var/lib/apt/lists/*


# Update pip
# Install PyFlakes for testing
RUN python2 -m pip install --upgrade \ 
    pip \
    setuptools \
    pyflakes
RUN python3 -m pip install --upgrade \
    pip \
    setuptools \
    pyflakes



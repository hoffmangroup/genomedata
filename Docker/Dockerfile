FROM ubuntu:latest

# Install Ubuntu 18.04.01 LTS (Bionic Beaver)
# Install pkg-config for locationg development files (hdf5)
# Install python 2.7 and pip
# Install python 3.6 and pip3
# Install HDF5
# Install Mercurial for hgtools (needed to install for testing)
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
    build-essential \
    pkg-config \
    python \
    python3 \
    python3-pip \
    libhdf5-dev \
    hdf5-tools \
    git \
    curl \
&& apt-get clean \
&& rm -rf /var/lib/apt/lists/*

# Install pip for Python 2
RUN curl https://bootstrap.pypa.io/get-pip.py --output get-pip.py && python2 get-pip.py

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
FROM ubuntu:latest

# Install latest Ubuntu LTS 
# Install HDF5 and development libraries
# Install pkg-config for locationg development files (hdf5)
# Install Python 3 and development libraries
# Install node for github actions
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
    build-essential \
    pkg-config \
    python3 \
    python3-dev \
    python3-venv \
    python3-pip \
    libhdf5-serial-dev \
    hdf5-tools \
    git \
    nodejs \
&& apt-get clean \
&& rm -rf /var/lib/apt/lists/*

# Install PyFlakes for testing
RUN python3 -m pip install --upgrade \
    build \
    pip \
    setuptools

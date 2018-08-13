#Image for genomedata testing in Python 2 and 3
FROM ubuntu:latest

# Install Ubuntu 18.04.01 LTS (Bionic Beaver)
# Install python 2.7.12 and 3.6.5, and pip 2 and 3

RUN apt-get update && apt-get install -y \
    python2 \
    python3 \
    python2-pip \
    python3-pip \
    libhdf5-dev \
    hdf5-tools 

# Update pip
# Install mercurial for hgtools (needed to install for testing)
RUN python2 -m pip install --upgrade \ 
    pip \
    setuptools \
    mercurial \

RUN python3 -m pip install --upgrade \
    pip \
    setuptools \
    mercurial \



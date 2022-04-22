#!/usr/bin/env python

'''genomedata: tools for accessing large amounts of genomic data

Genomedata is a format for efficient storage of multiple tracks of
numeric data anchored to a genome. The format allows fast random
access to hundreds of gigabytes of data, while retaining a small disk
space footprint. We have also developed utilities to load data into
this format.
'''

from __future__ import absolute_import, division, print_function

# Copyright 2008-2014 Michael M. Hoffman <michael.hoffman@utoronto.ca>

import errno
import os
from shlex import split
import sys
import tokenize

from distutils import sysconfig
from setuptools import Extension, find_packages, setup
from subprocess import CalledProcessError, check_output

DEFAULT_SHELL_ENCODING = "ascii"
LDFLAGS_LIBRARY_SWITCH = "-l"
LDFLAGS_LIBRARY_PATH_SWITCH = "-L"
CFLAGS_INCLUDE_PATH_SWITCH = "-I"

MINIMUM_PYTHON_VERSION_STR = "3.7"

doclines = __doc__.splitlines()
name, short_description = doclines[0].split(": ")
long_description = "\n".join(doclines[2:])

url = "https://pmgenomics.ca/hoffmanlab/proj/%s/" % name.lower()
download_url = "https://pypi.python.org/pypi/%s" % name.lower()

classifiers = ["Natural Language :: English",
               "Development Status :: 5 - Production/Stable",
               "Intended Audience :: Science/Research",
               "License :: OSI Approved :: GNU General Public License v2 "
               "(GPLv2)",
               "Topic :: Scientific/Engineering :: Bio-Informatics",
               "Operating System :: Unix",
               "Programming Language :: Python :: 3"]

entry_points = """
[console_scripts]
genomedata-info = genomedata._info:main
genomedata-query = genomedata._query:main
genomedata-histogram = genomedata._histogram:main
genomedata-load = genomedata.load_genomedata:main
genomedata-load-data = genomedata._load_data:main
genomedata-load-seq = genomedata._load_seq:main
genomedata-load-assembly = genomedata._load_seq:main
genomedata-open-data = genomedata._open_data:main
genomedata-hardmask = genomedata._hardmask:main
genomedata-close-data = genomedata._close_data:main
genomedata-report = genomedata._report:main
genomedata-erase-data = genomedata._erase_data:main
"""

setup_requires = ["setuptools_scm"]  # source control management packaging
# Exclude PyTables 3.4.1 - incorrect binary distribution causes core dumps
# See:
# https://bitbucket.org/hoffmanlab/genomedata/issues/38/pytables-341-causes-a-core-dump-when
# path.py 11 renames 'path' to 'Path'
install_requires = ["numpy", "tables>=3.0,!=3.4.1", "six",
                    "textinput>=0.2.0", "path.py>=11", "pyBigWig"]

# Monkey patches tokenize.detect_encoding() to return a blank string when it
# can't recognize encoding setuptools attempts to process some of the C files
# present, and errors because it can't determine encoding
try:
    _detect_encoding = tokenize.detect_encoding
except AttributeError:
    pass
else:
    def detect_encoding(readline):
        try:
            return _detect_encoding(readline)
        except SyntaxError:
            return "", []
    tokenize.detect_encoding = detect_encoding

source_files = ["src/_c_load_data.c"]
# sz may be needed here if it's statically builtin with an hdf5 distribution?
# Or someone built their own hdf5 version with sz in which case they should
# rely on using LD_LIBRARY_PATH?
libs = ["hdf5", "m", "z"]  # hdf5, math, zlib, (sz?)
library_dirnames = []
include_dirnames = [
    sysconfig.get_config_var("INCLUDEDIR"),  # environment headers
]
c_define_macros = [("H5_NO_DEPRECATED_SYMBOLS", None)]

if "HDF5_DIR" in os.environ:
    hdf5_dir = os.environ["HDF5_DIR"]
    include_dirnames.append(os.path.join(hdf5_dir, "include"))
    library_dirnames.append(os.path.join(hdf5_dir, "lib"))

# Attempt to get HDF5 development directories through pkg-config
try:
    shell_encoding = sys.stdout.encoding
    # Depending on the shell, python version (2), and environment this is not
    # guaranteed to be set. Attempt to fall back to a best-guess if it is not
    # set
    if not shell_encoding:
        shell_encoding = DEFAULT_SHELL_ENCODING

    pkg_config_cflags = split(check_output(
                                ["pkg-config", "--cflags", "hdf5"]
                             ).decode(shell_encoding))
    pkg_config_libs = split(check_output(
                                ["pkg-config", "--libs", "hdf5"]
                           ).decode(shell_encoding))
except OSError as err:
    # OSError ENOENT occurs when pkg-config is not installed
    if err.errno == errno.ENOENT:
        pass
    else:
        raise err
except CalledProcessError:
    # CalledProcessError occurs when hdf5 is not found by pkg-config
    pass
else:
    for path in pkg_config_cflags:
        assert path.find(CFLAGS_INCLUDE_PATH_SWITCH) == 0
        # NB for Python >= 3.9, should swap partition to removeprefix
        include_dirnames.append(path.partition(CFLAGS_INCLUDE_PATH_SWITCH)[2])
    for word in pkg_config_libs:
        if not word.startswith(LDFLAGS_LIBRARY_SWITCH):
            assert word.startswith(LDFLAGS_LIBRARY_PATH_SWITCH)
            library_dirnames.append(
                word.partition(LDFLAGS_LIBRARY_PATH_SWITCH)[2])


load_data_module = Extension(
    '_c_load_data',  # needs to match C file PyInit definition
    sources=source_files,
    include_dirs=include_dirnames,
    libraries=libs,
    library_dirs=library_dirnames,
    define_macros=c_define_macros)


if __name__ == "__main__":
    setup(name=name,
          use_scm_version=True,
          description=short_description,
          author="Michael Hoffman",
          author_email="michael.hoffman@utoronto.ca",
          url=url,
          download_url=download_url,
          classifiers=classifiers,
          long_description=long_description,
          setup_requires=setup_requires,
          install_requires=install_requires,
          zip_safe=False,
          # XXX: this should be based off of __file__ instead
          packages=find_packages("."),  # including "test"
          include_package_data=True,
          entry_points=entry_points,
          # place extension in the base genomedata package
          ext_package="genomedata",
          ext_modules=[load_data_module],
          python_requires=">={}".format(MINIMUM_PYTHON_VERSION_STR)
          )

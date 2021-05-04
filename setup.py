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

import os
import sys
import tokenize

from distutils import sysconfig
from setuptools import Extension, find_packages, setup

if (sys.version_info[0] == 2 and sys.version_info[1] < 7) or \
   (sys.version_info[0] == 3 and sys.version_info[1] < 4):
    print("Genomedata requires Python version 2.7 or 3.4 or later")
    sys.exit(1)

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
               "Programming Language :: Python :: 2.7",
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

setup_requires = ["setuptools_scm"] # source control management packaging
# Exclude PyTables 3.4.1 - incorrect binary distribution causes core dumps
# See:
# https://bitbucket.org/hoffmanlab/genomedata/issues/38/pytables-341-causes-a-core-dump-when
# path.py 11 renames 'path' to 'Path'
install_requires = ["numpy", "tables>=3.0,!=3.4.1", "six",
                    "textinput>=0.2.0", "path.py>=11"]

# Monkey patches tokenize.detect_encoding() to return a blank string when it can't recognize encoding
# setuptools attempts to process some of the C files present, and errors because it can't determine encoding
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

source_files = ["src/genomedata_load_data.c"]
# sz may be needed here if it's statically builtin with an hdf5 distribution? Or someone built their own hdf5 version with sz in which case they should rely on using LD_LIBRARY_PATH?
libs = ["hdf5", "m", "z"] # hdf5, math, zlib, (sz? lossless compression libray for scientific data)
library_dirnames = []
include_dirnames = [
    sysconfig.get_config_var("INCLUDEDIR"), # environment headers
]
c_define_macros = [("H5_NO_DEPRECATED_SYMBOLS", None)]

if "HDF5_DIR" in os.environ:
    hdf5_dir = os.environ["HDF5_DIR"]
    include_dirnames.append(os.path.join(hdf5_dir, "include"))
    library_dirnames.append(os.path.join(hdf5_dir, "lib"))

load_data_module = Extension('_load_data_c_ext', # needs to match C file PyInit definition
                            sources=source_files,
                            include_dirs = include_dirnames,
                            libraries = libs,
                            library_dirs = library_dirnames,
                            define_macros = c_define_macros,
                            extra_compile_args = ["-UNDEBUG"] # NB: Keep assert macros functioning
                            )


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
          ext_package= "genomedata", # place extension in the base genomedata package
          ext_modules=[load_data_module]
          )

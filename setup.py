#!/usr/bin/env python

# Copyright 2008-2022 Michael M. Hoffman <michael.hoffman@utoronto.ca>

import errno
import os
from shlex import split
import sys

from distutils import sysconfig
from setuptools import Extension, setup
from subprocess import CalledProcessError, check_output

DEFAULT_SHELL_ENCODING = "ascii"
LDFLAGS_LIBRARY_SWITCH = "-l"
LDFLAGS_LIBRARY_PATH_SWITCH = "-L"
CFLAGS_INCLUDE_PATH_SWITCH = "-I"

source_files = ["src/_c_load_data.c"]
# sz may be needed here if it's statically builtin with an hdf5 distribution?
# Or someone built their own hdf5 version with sz in which case they should
# rely on using LD_LIBRARY_PATH?
libs = ["hdf5", "m", "z"]  # hdf5, math, zlib, (sz?)
library_dirnames = []
include_dirnames = [
    sysconfig.get_config_var("INCLUDEDIR"),  # environment headers
]

# Build against the stable ABI to support 3.9 onwards
c_define_macros = [("Py_LIMITED_API", "0x030A0000")]

# If possible, use HDF5_DIR environment variable as preferred library source
if "HDF5_DIR" in os.environ:
    hdf5_dir = os.environ["HDF5_DIR"]
    include_dirnames.append(os.path.join(hdf5_dir, "include"))
    library_dirnames.append(os.path.join(hdf5_dir, "lib"))
# Otherwise attempt to get HDF5 development directories through pkg-config
else:
    try:
        shell_encoding = sys.stdout.encoding
        # Depending on the shell, python version (2), and environment this is
        # not guaranteed to be set. Attempt to fall back to a best-guess if it
        # is not set
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
            # NB for Python >= 3.9, should swap partition to remove prefix
            include_dirnames.append(path.partition(
                                    CFLAGS_INCLUDE_PATH_SWITCH)[2])
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
    define_macros=c_define_macros,
    py_limited_api=True)


if __name__ == "__main__":
    # place extension in the base genomedata package
    setup(ext_package="genomedata",
          ext_modules=[load_data_module])
